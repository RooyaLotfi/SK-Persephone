# last updated: 2020/12/30
import datetime  # For date_to_day_of_year
import os  # For reproject_to_file
import pathlib
from typing import Union  # For reproject_to_file

import geopandas as gpd
import numpy as np  # For reproject_tif, reproject_to_file
import numpy.ma as ma  # For reproject_tif ; For mask_data
import rasterio  # For reproject_tif; For reproject_to_file; For mask_data
import scipy.interpolate  # For interpolate
from osgeo import gdal, ogr  # For overlap_detection
# from rasterio.crs import CRS  # For reproject_tif
from rasterio.mask import mask  # For mask_data
from rasterio.warp import \
    calculate_default_transform  # For reproject_to_file; For reproject_tif
from rasterio.warp import Resampling, reproject
from skimage.transform import \
    resize as skimage_resample  # For the resample portion of mask_data

# from osgeo.utils.gdal_edit import * # For whatever reason it does not exist in the GDAL folder when installed through Conda forge. For now I
# copied the python script into a separate .py called gdal_edit and will import it as is
import gdal_edit


def overlap_detection(raster_dir, vector_dir):
    """
    Detects if there is an overlap between a raster and a vector file (e.g. a shapefile)
    Args:
        raster_dir: Directory for the raster file, as a string
        vector_dir: Directory for the vector file, as a string
    Returns:
        bool: Boolean value indicating if there exists an overlap
    """
    raster = gdal.Open(raster_dir)
    vector = ogr.Open(vector_dir)
    if raster is None:
        raise ValueError('overlap_detection: Could not find raster file: ' + str(raster_dir))
    if vector is None:
        raise ValueError('overlap_detection: Could not find vector file: ' + str(vector_dir))

    # Get raster geometry
    transform = raster.GetGeoTransform()
    pixelWidth = transform[1]
    pixelHeight = transform[5]
    cols = raster.RasterXSize
    rows = raster.RasterYSize

    xLeft = transform[0]
    yTop = transform[3]
    xRight = xLeft + cols * pixelWidth
    yBottom = yTop + rows * pixelHeight

    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(xLeft, yTop)
    ring.AddPoint(xLeft, yBottom)
    ring.AddPoint(xRight, yBottom)
    ring.AddPoint(xRight, yTop)
    ring.AddPoint(xLeft, yTop)

    rasterGeometry = ogr.Geometry(ogr.wkbPolygon)
    rasterGeometry.AddGeometry(ring)

    # Get vector geometry
    layer = vector.GetLayer()
    feature = layer.GetFeature(0)
    vectorGeometry = feature.GetGeometryRef()
    overlap = rasterGeometry.Intersect(vectorGeometry)
    return overlap


def interpolate(data, kind='linear', start_date=1, end_date=365, as_dict=False):
    """
    Interpolates the input array, pixel by pixel, across the first dimension that is assumed are the days.
    Assumes that if the data comes as a numpy array, that it is 1-based as in the Day-of-Year calendar, so the value for the first position in the array (0) is day 1.
    If the data comes as a dictionary, it uses the same value as that of the key for each day.
    Limitations: It is only able of interpolating between the first date and the last date of the data available.
    Args:
        data: 2 or 3-dimensional numpy array, the first dimension being the days, so it should look like days x width x height or days x width
        kind: Optional argument, string indicating the kind of interpolation to perform. Based directly on Scipy's interp1d function's argument:
            (‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’, ‘next’, where ‘zero’, ‘slinear’, ‘quadratic’ and ‘cubic’
            refer to a spline interpolation of zeroth, first, second or third order;‘previous’ and ‘next’ simply return the previous or next
            value of the point) or as an integer specifying the order of the spline interpolator to use.
            Default is ‘linear’.
        start_date: Start date in the Julian calendar format, as an integer, for the first value of the interpolation. Default is 1
        end_date: End date in the Julian calendar format, as an integer, for the last value of the interpolation. Default is 365
        as_dict: Boolean that for returning the resulting numpy array as a dictionary instead, with the keys being the days and the values being the resulting numpy arrays for each day
    Returns:
        numpy.ma / dict: 3-dimensional numpy masked array or a dictionary with days for the keys and 2D or 1D arrays for the values, if as_dict is True.
    """
    # Check if the data is a dictionary and convert it to a numpy array

    if isinstance(data, dict):
        # Check that all the arrays in the dictionary are of the same size
        data_keys = list(data.keys())
        for i in range(len(data_keys) - 1):
            if data[data_keys[i]].shape != data[data_keys[i + 1]].shape:
                raise ValueError('interpolate: the shape of the arrays in the dictionary do not match')
        data_shape = ((365,) + data[data_keys[
            -1]].shape)  # Create the shape of the final array (days by width x height or days by width/height)
        data_array = np.ma.array(
            np.full(data_shape, np.nan))  # Create a numpy array filled with np.nan, the size of the data
        for key in data_keys:
            data_array[int(key) - 1] = data[key]  # Add the data of the dictionary into the array
    elif isinstance(data, np.ndarray):
        data_array = np.ma.array(
            data.filled().copy())  # Fill the data if it is a numpy masked array, and copy it to the data_array
    else:
        raise TypeError('interpolate: data needs to be of type dict or numpy array, received: ' + str(type(data)))

    data_array = ma.masked_invalid(data_array)  # Mask invalid data
    ma.set_fill_value(data_array, np.nan)
    indices = np.ndindex(data_array.shape[1:])  # Get the indices
    moved_data = np.moveaxis(data_array, 0, -1)  # Move the days axis to the end of the list
    interpolated_data = np.full(moved_data.shape, np.nan)

    for index in indices:
        interpol_data = moved_data[index].filled().copy()  # Slice the data as a filled copy
        if not all(np.isnan(interpol_data)):  # Check if all the values of the interpolation_data are nan
            day_range = np.arange(start_date, end_date + 1)
            day_range = day_range[~np.isnan(interpol_data)]
            interpolation_y = interpol_data[~np.isnan(interpol_data)]

            y_interp = scipy.interpolate.interp1d(day_range, interpolation_y,
                                                  kind=kind)  # Create the interpolation function

            interpol_dates = np.arange(min(day_range),
                                       max(day_range) + 1)  # Get the maximum and minimum interpolation dates possible
            interpols = y_interp(interpol_dates)  # Get the interpolations
            interpol_data[interpol_dates - 1] = interpols  # Store the interpolations

        interpolated_data[index] = interpol_data  # Save the resulting data
    interpolated_data = np.moveaxis(interpolated_data, -1, 0)  # Reinsert the days axis at the beginning of the list
    interpolated_data = ma.masked_invalid(interpolated_data)  # Mask all the invalid data (i.e. np.nan)
    ma.set_fill_value(interpolated_data, np.nan)  # Set the fill value of the masked array

    if as_dict:
        interpolated_data_dict = {}
        # TODO (Juan): The following part is using day_range, which is the last value found inside of the indices and not the general available data for the functiom
        for i in range(min(day_range), max(day_range) + 1):
            interpolated_data_dict[i] = interpolated_data[i - 1]  # Pass the data into a dictionary
        interpolated_data = interpolated_data_dict  # Replace the numpy array for the dictionary
    return interpolated_data


def date_to_day_of_year(date, format='%Y-%m-%d', return_as_string=True):
    """
    Converts a date string with a specific format into a julian date. If return_as_string is set, the format
    will be with the leading zeroes to make sure the return is always three characters long
    Args:
        date: A string with the date information
        format: the format in which the string is
        return_as_string: Specifies if the function should return a string or an int.
    Returns:
        int/str: The julian date for the specified value, either in int or string. If -1, it means it couldn't convert the date
    """
    julius = -1
    if isinstance(date, str):
        dt = datetime.datetime.strptime(date, format)
        tt = dt.timetuple()
        if return_as_string:
            if tt.tm_yday < 10:
                julius = '00' + str(tt.tm_yday)
            elif tt.tm_yday < 100:
                julius = '0' + str(tt.tm_yday)
            else:
                julius = str(tt.tm_yday)
        else:
            julius = tt.tm_yday
    return julius


def metadata_add_to_raster(raster_path: str, metadata: Union[dict, list, str]) -> int:
    """
    Adds metadata to the raster, be it one value, or a list or dictionary of values
    Args:
        raster_path: path to the raster file
        metadata: values to add as metadata to the raster. The function automatically appends
            the required prefixes to be able to be loaded using GDAL_edit. If a list, it assumes
            that the list contains strings in the form: 'key=value'. If in a dictionary, the
            function assumes that the dictionary keys are the keys to be used for the metadata input
    Returns:
        int: the result of calling the GDAl_edit function. If successful returns 0, otherwise returns -1
    """
    if isinstance(metadata, list):
        formatted_metadata = []
        for value in metadata:
            formatted_metadata.extend(['-mo', value])
    elif isinstance(metadata, str):
        formatted_metadata = ['-mo', metadata]
    elif isinstance(metadata, dict):
        formatted_metadata = []
        for key in metadata.keys():
            formatted_metadata.extend(['-mo', '='.join([str(key), str(metadata[key])])])
    else:
        formatted_metadata = None
    result = -1
    if formatted_metadata is not None:
        gdal_edit_options = ['python_expected_first_value.py']
        gdal_edit_options.extend(formatted_metadata)
        gdal_edit_options.append(raster_path)
        result = gdal_edit.gdal_edit(gdal_edit_options)
    return result


def metadata_delete_from_raster(geotiff_file: str):
    """
    Deletes all the metadata from a raster
    Args:
        geotiff_file: Path to the raster file for which to remove the metadata
    """
    gdal_edit_options = ['argvs_first_value_is_assumed_to_be_script_name', '-unsetmd', geotiff_file]
    gdal_edit.gdal_edit(gdal_edit_options)


def metadata_get_from_raster(raster_path: str) -> dict:
    """
    Gets all the metadata stored in a raster
    Args:
        raster_path: path to the raster
    Returns:
        dict: metadata of the raster, as a dictionary
    """
    dataset = gdal.Open(raster_path)
    meta = dataset.GetMetadata()
    return meta


def reproject_tif(tif_dir, dst_crs):
    """
    Reprojects a Geotiff to a specific CRS and returns it as a numpy masked array.
    Args:
        tif_dir -- Directory for where to find the Geotiff
        dst_crs -- Destination CRS for the reprojection of the Geotiff file
    Returns:
        numpy.ma -- Numpy masked array with the reprojected Geotiff
    """
    # Based on https://gis.stackexchange.com/a/261378
    # Example of a dst_crs = 'EPSG:32614'

    with rasterio.open(tif_dir) as src:
        smap_tif = src.read(1, masked=True)
        src_crs = src.crs
        src_transform = src.transform
        src_bounds = src.bounds
        src_width = src.width
        src_height = src.height

    # Prepare destination array for warping
    dst_transform, dst_width, dst_height = calculate_default_transform(src_crs,
                                                                       dst_crs,
                                                                       src_width,
                                                                       src_height,
                                                                       *src_bounds,
                                                                       )
    newCRS = np.ma.zeros((dst_height, dst_width), dtype=np.float32)
    ma.set_fill_value(newCRS, src.meta['nodata'])

    # reproject array
    reproject(source=smap_tif,
              destination=newCRS,
              src_transform=src_transform,
              src_crs=src_crs,
              dst_transform=dst_transform,
              dst_crs=dst_crs,
              dst_nodata=src.meta['nodata'],
              resampling=Resampling.bilinear
              )

    newCRS = ma.masked_values(newCRS, src.meta['nodata'])

    """ For plotting a comparison beetween both files:
    fig, axs = plt.subplots(1,2, figsize=(20,12))
    ep.plot_bands(smap_tif,
                cmap="viridis",
                ax=axs[0])
    ep.plot_bands(newCRS,
    cmap="viridis",
    ax=axs[1])
    plt.show()
    """
    return newCRS


def reproject_to_file(tif_dir, dst_crs, shape=None, dst_fldr=None):
    """
    Reprojecting a Geotiff and saving it to file, using Rasterio.
    Args:
        tif_dir: Directory for where to find the Geotiff
        dst_crs: Destination CRS for the reprojection of the Geotiff file
        shape: Tuple with the desired height and width (in that order) for the reprojected Geotiff
        dst_fldr: Folder for where to store the resulting Geotiff. If empty, the output folder will be used as the base

    Returns:
        str: Destination filename for where the reprojection file was written to
    """

    # Example of a dst_crs = 'EPSG:32614'
    file_name = os.path.basename(tif_dir)  # Get the filename
    if dst_fldr is None:
        dst_fldr = os.path.join('output', os.path.dirname(tif_dir))
    pathlib.Path(dst_fldr).mkdir(parents=True, exist_ok=True)

    with rasterio.open(tif_dir) as src:
        if shape is None:
            dst_height = None
            dst_width = None
        else:
            dst_height, dst_width = shape
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds,
            dst_width=dst_width, dst_height=dst_height)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        dst_dir = os.path.join(dst_fldr, file_name)
        with rasterio.open(dst_dir, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest)
    return dst_dir


def reproject_to_file_gdal(tif_dir, dst_crs, dst_fldr=None):
    """
    Reprojecting a Geotiff and saving it to file, using gdal.
    Args:
        tif_dir: Directory for where to find the Geotiff
        dst_crs: Destination CRS for the reprojection of the Geotiff file
        dst_fldr: Folder for where to store the resulting Geotiff. If empty, the output folder will be used as the base

    Returns:
        str: Destination filename for where the reprojection file was written to
    """
    if dst_fldr is None:
        dst_crs_fname = dst_crs.replace(':', '_')
        dst_fldr = os.path.join('reprojected_to_' + dst_crs_fname, os.path.dirname(tif_dir))
    pathlib.Path(dst_fldr).mkdir(parents=True, exist_ok=True)

    input_raster = gdal.Open(tif_dir)
    output_raster = dst_fldr

    gdal.Warp(output_raster, input_raster, dstSRS=dst_crs)
    return dst_fldr


def mask_data(tif_dir, shapefile, crop_data=False, all_touched=True, resample=None, fill_value=None):
    """
    Masks the data from the Geotiff file based on an input shapefile.
    Args:
        tif_dir: Directory of the Geotiff file.
        shapefile: Geopandas geodataframe with a 'geometry' column.
        crop_data: If set to true, crops the data to the size/shape of the shapefile.
        all_touched: Include a pixel in the mask if it touches any of the shapes of the shapefile
        resample: Resample the array to the specified tuple
        fill_value: Set the fill value for the masks

    Returns:
        out_image_ma -- The masked array representation of the Geotiff file.
        out_info -- Information of the original Geotiff file.
    """
    with rasterio.open(tif_dir) as src:
        out_image, _ = mask(src, shapefile['geometry'], crop=crop_data, all_touched=all_touched)
        # _ can be replaced as out_transform and be used for mapping pixel coordinates of out_image to another coordinate system.
        out_info = src
    out_image_ma = np.squeeze(ma.array(out_image.astype('float')))

    if resample is not None:
        dst_width, dst_height = resample
        out_image_ma = skimage_resample(out_image_ma, (dst_height, dst_width), order=0)
        out_image_ma = np.ma.masked_array(out_image_ma)
        ma.set_fill_value(out_image_ma, src.meta['nodata'])
        out_image_ma = ma.masked_values(out_image_ma, src.meta['nodata'])

    if out_info.meta['nodata'] is not None:
        if np.issubdtype(out_image_ma.dtype, np.floating):
            out_image_ma = ma.masked_values(out_image_ma, out_info.meta[
                'nodata'])  # Better floating point approximation than np.where
        else:
            out_image_ma = ma.masked_equal(out_image_ma, out_info.meta['nodata'])

    if fill_value is not None:  # Change the fill value to the one specified in the argument
        ma.set_fill_value(out_image_ma, fill_value)
        src.meta[
            'nodata'] = fill_value  # Not relevant for now because the function is using out_info, but letting it here in case there's a change of functionality

    # out_image_ma = ma.masked_values(out_image_ma, src.meta['nodata'])
    return out_image_ma, out_info


def get_tif_info(tif_dir):
    with rasterio.open(tif_dir) as src:
        out_info = src.profile
    return out_info


def check_tif_reprojection(shapefile_dir, tif_dir):
    """
    Check if the shapefile and tif file overlap, if not check if it's because they're in different CRS. If they have different CRS,
    check if there's already a reprojected version of the file, and if not, create it.

    Args:
        shapefile_dir: Shapefile directory location
        tif_dir: Tif directory location
    Returns:
        str: Output directory of the tif file, either if it previously existed or if the function reprojected and created the file
    """
    # TODO: Get a better function name for this, can't think of a better description right now
    if not overlap_detection(tif_dir, shapefile_dir):
        shapefile = gpd.read_file(shapefile_dir)
        # Files don't overlap
        if not shapefile.crs == get_tif_info(tif_dir)[
            'crs']:  # Files don't overlap, check if it's because of a different CRS
            # Reproject the files to be able to use them
            # First, Check if the file has already been created in the output folder, there is an overlap and it is in the correct CRS
            output_dir = (os.path.join('output', tif_dir))
            if not os.path.isfile(output_dir):  # Check if reprojected file exists
                # if not overlap_detection(output_dir, shapefile_dir): # Check if the reprojected file doesn't overlap with the shapefile. NOT WORKING AS INTENDED
                reproject_to_file(tif_dir, shapefile.crs)  # File doesn't exist, create its reprojection
        else:
            raise ValueError("Files have same CRS but don't overlap")
    else:
        output_dir = tif_dir
    return output_dir


def resample(image, width, height, interpolation='NEAREST'):
    import cv2  # For resize
    """
    Resample an image to the specified width and height, using the specified interpolation method.

    Args:
        layer: The image to resample
        width: Desired width
        height: Desired height.
        interpolation: Interpolation method, as a string. The options are the same as of OpenCV:
            INTER_NEAREST – a nearest-neighbor interpolation 
            INTER_LINEAR – a bilinear interpolation
            INTER_AREA – resampling using pixel area relation. It may be a preferred method for image decimation, as it gives moire’-free results. But when the image is zoomed, it is similar to the INTER_NEAREST method. 
            INTER_CUBIC – a bicubic interpolation over 4×4 pixel neighborhood 
            INTER_LANCZOS4 – a Lanczos interpolation over 8×8 pixel neighborhood
    Returns:
        numpy.ma: The resampled image, as a numpy masked array
    """

    if isinstance(interpolation, str):
        if interpolation == 'NEAREST':
            interpol_method = cv2.INTER_NEAREST
        elif interpolation == 'LINEAR':
            interpol_method = cv2.INTER_LINEAR
        elif interpolation == 'AREA':
            interpol_method = cv2.INTER_AREA
        elif interpolation == 'CUBIC':
            interpol_method = cv2.INTER_CUBIC
        elif interpolation == 'LANCZOS4':
            interpol_method = cv2.INTER_LANCZOS4
        else:
            raise ValueError
    else:
        raise TypeError

    img = image
    img_fill = img.fill_value
    img = img.filled()

    dim = (width, height)
    # resize image
    resampled = cv2.resize(img, dim, interpolation=interpol_method)
    resampled = ma.masked_values(resampled, img_fill)
    return resampled


if __name__ == '__main__':
    # Example of the interpolation function and testing the resulting interpolation

    # First let's create a test dictionary with only even numbers, and a full array with what are the expected results from the test dictionary
    interpolate_test_array = np.full((365, 100, 100), np.nan)
    interpolate_test_dict = {}
    for day in range(1, 366):
        dat = np.full((100, 100), day)
        interpolate_test_array[day - 1] = dat
        if day % 2 == 0:
            interpolate_test_dict[day] = dat

    # Remove the values from the first day and last days, as the test data doesn't have them and the interpolation function wouldn't calculate them
    interpolate_test_array[0] = np.full((100, 100), np.nan)
    interpolate_test_array[364] = np.full((100, 100), np.nan)
    interpolate_test_array = ma.masked_invalid(interpolate_test_array)  # Mask all the invalid data (i.e. np.nan)
    ma.set_fill_value(interpolate_test_array, np.nan)  # Set the fill value for the masked array

    data = interpolate(interpolate_test_dict, kind='linear')  # interpolate the data linearly

    # Compare the results
    for i in range(365):
        if not np.array_equal(interpolate_test_array[i].filled(), data[i].filled(), equal_nan=True):
            print('arrays different in line ', str(i), ':', str(interpolate_test_array[i]), 'vs', str(data[i]))
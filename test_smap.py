import os
import geopandas as gpd
import rioxarray as rxr
from Utils.crop_soilgrid_smap import crop
import matplotlib.pyplot as plt
from Utils import crop_soilgrid_smap
import numpy as np
import rioxarray as rxr
from utilities.FileUtils import FileUtils
from rasterio.enums import Resampling
from copy import copy
from utilities.geotiff_utils import interpolate
import utilities.geotiff_utils as geotiff_utils
from scipy import interpolate
from datetime import date


def interpolate_spatial(array, method='cubic'):
    interpolated_array = []
    for band in array:
        x = np.arange(0, band.shape[1])
        y = np.arange(0, band.shape[0])
        # mask invalid values
        band = np.ma.masked_invalid(band)
        xx, yy = np.meshgrid(x, y)
        # get only the valid values
        x1 = xx[~band.mask]
        y1 = yy[~band.mask]
        newarr = band[~band.mask]

        interpolated_band = interpolate.griddata((x1, y1), newarr.ravel(), (xx, yy), method=method)
        interpolated_array.append(interpolated_band)
    return np.array(interpolated_array)


def smap_mask(origin_path, destination_path, RM_IDs, crop_type_id, years, months, days_num, features, resolution):
    CRS = "EPSG:32613"
    buffer = 40000
    f = FileUtils()
    rm_polygons_path = os.path.join("polygon_EPSG_32613/cropped_by_CA_boundaries/", "polygon.shp")
    rm_polygons = gpd.read_file(rm_polygons_path)

    # Read AAFC data
    crop_type_path = os.path.join("AAFC_data", "aci_2019_sk_v1", "masked_cropID__146_153.tif")
    crop_type = rxr.open_rasterio(crop_type_path)

    for rm_id in RM_IDs:
        # Pick a polygon
        rm_polygon = rm_polygons.loc[rm_polygons['id'] == rm_id]
        # Make a directory based on RM id and the crop type
        output_raster_path = os.path.join(destination_path, str(rm_id), "crop_" + str(crop_type_id))
        if not os.path.exists(output_raster_path):
            os.makedirs(output_raster_path)

        # Add buffer to existing polygon
        rm_polygon_buff = rm_polygon.buffer(buffer)

        # Crop the crop mask by polygon
        cropped_crop_type = crop(crop_type, rm_polygon)
        # Save crop mask to a directory
        crop_type_path = os.path.join("data", "SMAP", "3_crop", str(rm_id))
        if not os.path.exists(crop_type_path):
            os.makedirs(crop_type_path)

        crop_type_image = "crop_mask_" + str(rm_id) + '.tif'

        cropped_crop_type.rio.to_raster(os.path.join(crop_type_path, crop_type_image))
        print("Saved crop masks done!", os.path.join(crop_type_path, crop_type_image))

        # Create mask based on crop type
        crop_type_np = cropped_crop_type.data[0]
        crop_type_np[crop_type_np != crop_type_id] = 0
        crop_type_np[crop_type_np == crop_type_id] = 1

        for root, dirs, files in os.walk(origin_path):
            for file in files:
                print("does this file need to be processed ", file)
                arg_0, feature, arg_2, year, month, day = f.parse_smap_file(file)
                if year in years and int(month) in months and int(day) in days_num and feature in features:
                    print("YES! processing file :", file)
                    image_file_path = os.path.join(root, file)
                    image_file = rxr.open_rasterio(image_file_path)
                    # Change CRS
                    image_file = image_file.rio.reproject(dst_crs=CRS)
                    # Cropping with buffered polygon
                    cropped_image_file_buff = crop(image_file, rm_polygon_buff)
                    # Upsampling to 30m
                    cropped_image_file_buff = cropped_image_file_buff.rio.reproject(dst_crs=CRS, resolution=resolution,
                                                                                    resampling=Resampling.bilinear)
                    # Cropping with normal polygons
                    cropped_raster = crop(cropped_image_file_buff, rm_polygon)

                    # A temp raster with size of crop_type
                    temp_raster = copy(cropped_raster)
                    # TODO: Check why the shapes do not match and try to fix it
                    common_x_len = min(np.shape(temp_raster.data)[1], np.shape(cropped_crop_type.data)[1])
                    common_y_len = min(np.shape(temp_raster.data)[2], np.shape(cropped_crop_type.data)[2])

                    for data in temp_raster.data:
                        # data = data*crop_type_np
                        data[0:common_x_len, 0:common_y_len] = data[0:common_x_len, 0:common_y_len] * \
                                                               crop_type_np[0:common_x_len, 0:common_y_len]
                    save_path = os.path.join(output_raster_path, file)
                    temp_raster.rio.to_raster(save_path)
                    print("Save Done! :", save_path)


def get_key(year, month, day):
    f_date = date(year, 1, 1)
    l_date = date(year, month, day)
    delta = l_date - f_date

    return delta.days


def smap_interpolate(origin_path, destination_path, RM_IDs, years, sat_ids, crop_type_id, features, nan_values):
    """
    interpolate function gets a dictionary as: keys = days values = numpy
    this function creates this dictionary
    :param origin_path: path to all rasters
    :param destination_path: path to save data after interpolation
    :param RM_IDs: id's we want to interpolate
    :param years: 3019 ex
    :param sat_ids: s30 l30
    :param crop_type_id: 146 or 154 we need this parameter to find the directory
    :return:
    """
    f = FileUtils()
    for rm_id in RM_IDs:
        for year in years:
            for sat_id in sat_ids:
                hls_path = os.path.join(origin_path, sat_id, year, str(rm_id),
                                        "crop_" + str(crop_type_id))
                print(hls_path)
                # If the directory does not exist or the directory is empty do not proceed
                if not os.path.exists(hls_path) or not os.listdir(hls_path):
                    print(f"No RM with id {rm_id} exists")
                    continue
                save_path = os.path.join(destination_path, sat_id, year, str(rm_id),
                                         "crop_" + str(crop_type_id))
                if not os.path.exists(save_path):
                    os.makedirs(save_path)
                interpolate_test_dict = {}
                for feature in features:
                    for root, dirs, files in os.walk(hls_path):
                        for file in files:
                            if '.tif' in file:
                                image_path = os.path.join(root, file)
                                # print("image path is ", image_path)
                                sat_id, tile, year, day, band_no, arg3 = f.parse_preprocessed_file(file)
                                arg_0, feature_data, arg_2, year, month, day = f.parse_smap_file(file)

                                if feature_data == feature:
                                    raster = rxr.open_rasterio(image_path)
                                    data = raster.data
                                    data[data == nan_values] = np.NaN
                                    key = get_key(year, month, day)
                                    interpolate_test_dict[key] = data

                    interpolated_data = geotiff_utils.interpolate(interpolate_test_dict, kind='linear',
                                                                  as_dict=True)  # interpolate the data linearly
                    raster_temp = raster
                    for key in interpolated_data:
                        print(f"***********************{key}*************************")
                        print("BEFORE ", np.sum(np.isnan(interpolated_data[key].data)))
                        print("shape of interpolated data ", np.shape(interpolated_data[key]))
                        # interpolate spatially
                        interpolated_data[key] = interpolate_spatial(interpolated_data[key], method='linear')
                        print("AFTER ", np.sum(np.isnan(interpolated_data[key])))
                        raster_temp.data = interpolated_data[key]
                        path_save = os.path.join(save_path,
                                                 "sm" + "_" + feature + "_analysis_" + year + str(key) + "_" + '.tif')
                        raster_temp.rio.to_raster(path_save)
                        print("save done! ", path_save)


def main():
    print("*********************************** Implementing cropping to the "
          "extent of RM with buffer**********************************************")
    YEAR = ['2015', '2019']
    RM_IDs = [232, 322]
    resolution = 30
    CROP_TYPE_NO = 153
    CROP_TYPE_NO = 146
    features = ['rootzone', 'surface']
    days_num = np.arange(1, 32).tolist()
    months = np.arange(1, 12).tolist()
    smap_path = os.path.join("data", "smap", "1_original")
    smap_processed_path = os.path.join("data", "smap", "5_resampled_cropped_masked")
    smap_mask(smap_path, smap_processed_path, RM_IDs, CROP_TYPE_NO, YEAR, months, days_num, features, resolution)

    # print("*********************************** Interpolating **********************************************")
    # destination_path = os.path.join("data", "smap", "6_interpolate")
    # # # what is defined as nan for smap
    # nan_values = 3.40282 * (10 ** 38)
    # smap_interpolate(smap_processed_path, destination_path, RM_IDs, YEAR, RM_IDs, CROP_TYPE_NO, features, nan_values)


if __name__ == '__main__':
    main()

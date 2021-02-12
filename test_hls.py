from utilities.Preprocessing import Preprocessing
import os
import rasterio as rio
from matplotlib import pyplot
from utilities.Preprocessing import Preprocessing
import numpy as np
import geopandas as gpd
import rioxarray as rxr
from Utils.crop_soilgrid_smap import crop
from copy import copy
from rasterio.enums import Resampling
import cv2
from utilities.FileUtils import FileUtils
from utilities.geotiff_utils import interpolate
import utilities.geotiff_utils as geotiff_utils
import math
from scipy import interpolate


def hls_QA(origin_path, destination_path, sat_ids, tiles_name, years, days):
    """

    :param origin_path: hls data path before implementing qa layer: in form of : data/hls/year/l30 or S30/ 12UWB or any tile/ ...tif
    :param destination_path: path to save data after implementing qa layer: in form of: data/2_qa_masked/L30 or S30/ tile/ ...tif
    :param sat_ids: Satellite IDs that we want to implement QA layer to: S30 or L30 or both
    :param tiles_name: tiles to implement QA layer: 12UWB
    :param years: years to implement QA layer to: ['2019']
    :param days: string: days to implement QA layer to
    :return: does not return anything just saves data after implementing QA layer to a new directory
    """
    prepare = Preprocessing()
    tiles = [tile[1:] for tile in tiles_name]
    for tile in tiles:
        for year in years:
            for sat_id in sat_ids:
                hls_path = os.path.join(origin_path, tile, year, sat_id)
                # Quality Assessment layer masking
                prepare.apply_qa_mask_hls(hls_path, destination_path, sat_ids, tiles_name, years, days)


def hls_merge(origin_path, destination_path, RM_IDs, sat_ids, years, days):
    """

    :param origin_path: origin path to read files from
    :param destination_path: path to save data
    :param RM_IDs: which RM's we want to work with: We find the tiles that overlap this RM and merge those tiles
    :param sat_ids: S30 or L30 or both
    :param years: 2019 for example
    :param days: which days we want to merge
    :return:
    """
    mgrs_rm_polygons_path = os.path.join("sk-mgrs-rm", "sk-mgrs-rm.shp")
    mgrs_rm_polygon = gpd.read_file(mgrs_rm_polygons_path)
    prepare = Preprocessing()
    for rm_id in RM_IDs:
        mgrs_tiles = mgrs_rm_polygon.loc[mgrs_rm_polygon['id'] == rm_id]['MGRS']
        if mgrs_tiles.empty:
            continue
        tiles_list = ["T" + tile for tile in mgrs_tiles]
        tiles_path = '-'.join(tiles_list)
        for year in years:
            for sat_id in sat_ids:
                for day in days:
                    prepare.merge(origin_path, destination_path, tiles_path, sat_id, mgrs_tiles, year, day)


def hls_mask(origin_path, destination_path, RM_IDs, crop_type_id, sat_ids, years, resolution):
    """
    You have some RM_IDs and you want to mask based on crop type: need to read AAFC data from the path we stored it
    Find the boundary of the polygon associated with RM_IDs:
    1. change the crs of data to fit the crs of RMs
    2. resample the hls data to 30m*30m
    3. crop the aafc data and crop the hls data with the RMs polygons
    4. mask the aafc data with crop_type_id value and keep the ones that have values we want and put other values 0
    5. create a rioxarray object and copy the original raster into this one. Multiply the raster data with the crop type
    mask created in step 3
    :param origin_path: path to find merged data
    :param destination_path: path to save
    :param RM_IDs: id of RMs
    :param crop_type_id: 146: canola 153: spring wheat
    :param sat_ids: s30 or l30
    :param years: 2019 for ex
    :param resolution: 30m we need this because when we mask we need to have the same size in both rasters
    :return:
    """
    CRS = "epsg:32613"
    rm_polygons_path = os.path.join("polygon_EPSG_32613", "cropped_by_CA_boundaries", "polygon.shp")
    rm_polygons = gpd.read_file(rm_polygons_path)

    mgrs_rm_polygons_path = os.path.join("sk-mgrs-rm", "sk-mgrs-rm.shp")
    mgrs_rm_polygon = gpd.read_file(mgrs_rm_polygons_path)

    # Read AAFC data
    crop_type_path = os.path.join("AAFC_data", "aci_2019_sk_v1", "masked_cropID__146_153.tif")
    crop_type = rxr.open_rasterio(crop_type_path)
    for sat_id in sat_ids:
        for rm_id in RM_IDs:
            mgrs_tiles = mgrs_rm_polygon.loc[mgrs_rm_polygon['id'] == rm_id]['MGRS']
            rm_polygon = rm_polygons.loc[rm_polygons["id"] == rm_id]
            if mgrs_tiles.empty or rm_polygon.empty:
                print(f"No RM with id {rm_id} exists")
                continue

            tiles_list = ["T" + tile for tile in mgrs_tiles]
            tiles_path = '-'.join(tiles_list)
            for year in years:

                # Crop AAFC
                # change crs, upsample to 30m * 30m
                crop_type = crop_type.rio.reproject(dst_crs=CRS, resolution=resolution,
                                                    resampling=Resampling.bilinear)

                # crop by polygon
                cropped_crop_type = crop(crop_type, rm_polygon)
                crop_type_path = os.path.join("data", "hls", "4_cropped", sat_id, str(year), str(rm_id))
                if not os.path.exists(crop_type_path):
                    os.makedirs(crop_type_path)
                crop_type_image = "crop_mask_" + str(rm_id) + '.tif'
                cropped_crop_type.rio.to_raster(os.path.join(crop_type_path, crop_type_image))
                print("Save done! ", crop_type_path, " ", crop_type_image)
                # Create mask based on crop type
                crop_type_np = cropped_crop_type.data[0]
                crop_type_np[crop_type_np != crop_type_id] = 0
                crop_type_np[crop_type_np == crop_type_id] = 1
                for sat_id in sat_ids:
                    merged_path = os.path.join(origin_path, sat_id, year, tiles_path)
                    output_raster_path = os.path.join(destination_path, sat_id, year, str(rm_id),
                                                      "crop_" + str(crop_type_id))
                    if not os.path.exists(output_raster_path):
                        os.makedirs(output_raster_path)
                    for root, dirs, files in os.walk(merged_path):
                        for file in files:
                            image_file_path = os.path.join(root, file)
                            data = rxr.open_rasterio(image_file_path)

                            # change crs of hls, upsample to 30m
                            data = data.rio.reproject(dst_crs=CRS, resolution=resolution,
                                                      resampling=Resampling.bilinear)
                            # Crop to the extent of RM without buffer
                            cropped_raster = crop(data, rm_polygon)

                            # A temp raster with size of crop_type
                            temp_raster = copy(cropped_raster)
                            # TODO: Check why the shapes do not match and try to fix it
                            common_x_len = min(np.shape(temp_raster.data)[1], np.shape(cropped_crop_type.data)[1])
                            common_y_len = min(np.shape(temp_raster.data)[2], np.shape(cropped_crop_type.data)[2])

                            for data in temp_raster.data:
                                data[0:common_x_len, 0:common_y_len] = data[0:common_x_len, 0:common_y_len] * \
                                                                       crop_type_np[0:common_x_len, 0:common_y_len]

                            temp_raster.rio.to_raster(os.path.join(output_raster_path, file))
                            print("Save Done! :", output_raster_path, " ", file)


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


def hls_interpolate(origin_path, destination_path, RM_IDs, years, sat_ids, crop_type_id, bands, nan_values):
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
                for band_number in bands:
                    for root, dirs, files in os.walk(hls_path):
                        for file in files:
                            if '.tif' in file:
                                image_path = os.path.join(root, file)
                                # print("image path is ", image_path)
                                sat_id, tile, year, day, band_no, arg3 = f.parse_preprocessed_file(file)
                                if band_no == band_number:
                                    raster = rxr.open_rasterio(image_path)
                                    raster = raster.astype('float64')
                                    data = raster.data
                                    data[data == nan_values] = np.NaN
                                    interpolate_test_dict[int(day)] = data

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
                                                 sat_id + "_" + year + "_" + str(key) + "_" + str(band_number) + '.tif')
                        raster_temp.rio.to_raster(path_save)
                        print("save done! ", path_save)


def get_hls_bands(start, end):
    """

    :param start: starting point of a band for example 1
    :param end: end point
    :return: bands stored in form of '01', ..., '10', '11', '12', '13'
    """
    bands_numbers = np.arange(start, end+1).tolist()
    bands = list(map(str, bands_numbers))
    for i in range(len(bands)):
        if len(bands[i]) == 1:
            bands[i] = '0' + bands[i]
    return bands


def main():
    # Apply QA layer
    print("*********************************** Implementing QA Layer**********************************************")
    qa_masked_path = os.path.join('data', 'hls', '2_qa_masked')
    origin_path = os.path.join('data', 'hls')
    sat_ids = ['L30', 'S30']
    tiles_name = ['T12UWC', 'T12UWB', 'T12UXC']
    years = ['2019']
    days_num = np.arange(149, 181).tolist()
    days = list(map(str, days_num))
    # hls_QA(origin_path, qa_masked_path, sat_ids, tiles_name, years, days)
    print("*********************************** Implementing Merge **********************************************")
    # Apply merge
    dir_region = os.path.join("data", "hls", "3_merged")
    RM_IDs = [13, 2, 322, 232, 292, 1]
    RM_IDs = [322, 322]
    # hls_merge(qa_masked_path, dir_region, RM_IDs, sat_ids, years, days)
    print("*********************************** Implementing Resample"
          ", Crop, Mask **********************************************")
    nan_values = -1000
    CROP_TYPE_NO = 146  # Canola
    CROP_TYPE_NO = 153  # Spring wheat
    resolution = 30
    merged_path = os.path.join("data", "hls", "3_merged")
    filtered_raster_path = os.path.join("data", "hls", "5_resampled_cropped_masked")
    # hls_mask(merged_path, filtered_raster_path, RM_IDs, CROP_TYPE_NO, sat_ids, years, resolution)

    print("*********************************** Interpolating **********************************************")
    processed_path = os.path.join("data", "hls", "5_resampled_cropped_masked")
    save_path = os.path.join("data", "hls", "6_interpolated")

    bands = get_hls_bands(1, 13)
    hls_interpolate(processed_path, save_path, RM_IDs, years, sat_ids, CROP_TYPE_NO, bands, nan_values)


if __name__ == "__main__":
    main()

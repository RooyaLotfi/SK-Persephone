from Utils.merge import merge
from utilities.qa_modis import make_LAI_FPAR_QC_mask
import os
import numpy as np
import rioxarray as rxr
from Utils.crop_soilgrid_smap import crop
from copy import copy
from rasterio.enums import Resampling
import cv2
from utilities.FileUtils import FileUtils
import geopandas as gpd
from test_hls import hls_interpolate


def modis_QA(origin_path, destination_path, sat_ids, tiles, bands, years, days_num, parameter_name, value):
    # TODO: create a dictionary to keep parameter and values as keys and values respectively this way you can manage to
    #  apply QA layer based on mutliple parameters and values
    """

    :param origin_path: raw rasters
    :param destination_path: paths to save
    :param sat_ids: sat IDs to implement QA layer to
    :param tiles: the tiles to implement QA layer to
    :param bands: the bands to implement qa layer to
    :param years: years of interest
    :param days_num: days of intererst
    :param parameter_name: which parameter to work with
    :param value: what values should be kept in that parameter
    :return:
    """
    fileutils = FileUtils()
    for year in years:
        qc_path_save = os.path.join(destination_path, year)
        if not os.path.exists(qc_path_save):
            os.makedirs(qc_path_save)
        for root, dirs, files in os.walk(origin_path):
            for file_mo in files:
                if '.tif' in file_mo:
                    sat_id_mo, tile_mo, year_mo, day_mo, band_no_mo = fileutils.parse_tile_file(file_mo)
                    if sat_id_mo in sat_ids and int(day_mo) in days_num and tile_mo in tiles and band_no_mo in bands:
                        print("band is", band_no_mo)
                        if band_no_mo == 'FparLai_QC':
                            modis_qc_dir = os.path.join(origin_path, year, file_mo)
                            modis_qc = rxr.open_rasterio(modis_qc_dir)
                            modis_qc_array = modis_qc.data[0]
                            filter_SCF_QG_0 = make_LAI_FPAR_QC_mask(modis_qc_array, parameter_name=parameter_name,
                                                                    value=value)
                            for file in files:
                                sat_id, tile, year, day, band_no = fileutils.parse_tile_file(file)
                                if sat_id == sat_id_mo and tile == tile_mo and year == year_mo and day == day_mo and \
                                        band_no != 'FparLai_QC' and band_no in bands:
                                    raster_path = os.path.join(root, file)
                                    raster = rxr.open_rasterio(raster_path)
                                    filtered_raster = copy(raster)
                                    filtered_raster.data[0] = filtered_raster.data[0] * filter_SCF_QG_0
                                    save_path = os.path.join(qc_path_save, file)
                                    filtered_raster.rio.to_raster(save_path)
                                    print("save path ", save_path)


def modis_merge(origin_path, destination_path, sat_ids, tiles, years, days):
    """
    gets sattelites that wants to implement merge to them and the day range and year range
    :param origin_path:
    :param destination_path:
    :param sat_id:
    :param tiles:
    :param years:
    :param days:
    :return:
    """
    region_name = '-'.join(tiles)
    for year in years:
        qa_mask_path_year = os.path.join(origin_path, year)
        for day in days:
            merge(path=qa_mask_path_year, dir_path=destination_path, region_name=region_name, sat_id=sat_ids[0],
                  tiles=tiles,
                  year=year,
                  day=day)


def modis_mask(origin_path, destination_path, sat_ids, RM_IDs, crop_type_id, years, days_num, resolution):
    """
    Gets the RM_IDs and finds the AAFC data and the rasters that fall in a day and year range and crops them to the
     extent of each RM and makes a mask based on crop type (crop_type_id) and filters the pixels that correspond
     to specific crop type
    :param origin_path:
    :param destination_path:
    :param RM_IDs:
    :param crop_type_id:
    :param years:
    :param days_num:
    :param resolution:
    :return:
    """
    CRS = "epsg:32613"
    fileutils = FileUtils()
    buffer = 500000
    rm_polygons_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone",
                                    "polygon_EPSG_32613/cropped_by_CA_boundaries/", "polygon.shp")
    rm_polygons = gpd.read_file(rm_polygons_path)

    for year in years:
        # Read AAFC data
        crop_type_path = os.path.join("AAFC_data", "aci_" + year + "_sk_v1", "masked_cropID__146_153.tif")
        crop_type = rxr.open_rasterio(crop_type_path)

        for rm_id in RM_IDs:
            rm_polygon = rm_polygons.loc[rm_polygons['id'] == rm_id]
            if rm_polygon.empty:
                print(f"No RM with id {rm_id} exists")
                continue

            # crop by polygon
            cropped_crop_type = crop(crop_type, rm_polygon)
            crop_type_path = os.path.join("data", "Modis_tif", "3_crop", year, str(rm_id))

            if not os.path.exists(crop_type_path):
                os.makedirs(crop_type_path)
            crop_type_image = "crop_mask_" + str(rm_id) + '.tif'
            cropped_crop_type.rio.to_raster(os.path.join(crop_type_path, crop_type_image))
            print("Saved crop masks done!", os.path.join(crop_type_path, crop_type_image))
            # Create mask based on crop type
            crop_type_np = cropped_crop_type.data[0]
            crop_type_np[crop_type_np != crop_type_id] = 0
            crop_type_np[crop_type_np == crop_type_id] = 1

            # Buffer around the rm polygon
            rm_polygon_buff = rm_polygon.buffer(buffer)
            if not os.path.exists(origin_path) or not os.listdir(origin_path):
                continue
            for root, dir, files in os.walk(origin_path):
                for file in files:
                    if '.tif' in file:
                        sat_id, tile, year, day, band_no, arg = fileutils.parse_preprocessed_file(file)
                        if year in years and int(day) in days_num:
                            output_raster_path = os.path.join(destination_path, sat_ids[0], year, str(rm_id),
                                                              "crop_" + str(crop_type_id))
                            if not os.path.exists(output_raster_path):
                                os.makedirs(output_raster_path)
                            image_file_path = os.path.join(root, file)
                            image_file = rxr.open_rasterio(image_file_path)
                            # Change CRS
                            image_file = image_file.rio.reproject(dst_crs=CRS)
                            # Cropping with buffered polygon
                            cropped_image_file_buff = crop(image_file, rm_polygon_buff)
                            # Upsampling to 30m
                            cropped_image_file_buff_reprojected = cropped_image_file_buff.rio.reproject(dst_crs=CRS,
                                                                                                        resolution=resolution,
                                                                                                        resampling=Resampling.bilinear)
                            # Cropping with normal polygons
                            cropped_raster = crop(cropped_image_file_buff_reprojected, rm_polygon)

                            # A temp raster with size of crop_type
                            temp_raster = copy(cropped_raster)
                            # TODO: Check why the shapes do not match and try to fix it
                            common_x_len = min(np.shape(temp_raster.data)[1], np.shape(cropped_crop_type.data)[1])
                            common_y_len = min(np.shape(temp_raster.data)[2], np.shape(cropped_crop_type.data)[2])

                            for data in temp_raster.data:
                                # data = data*crop_type_np
                                data[0:common_x_len, 0:common_y_len] = data[0:common_x_len, 0:common_y_len] * \
                                                                       crop_type_np[0:common_x_len, 0:common_y_len]

                            temp_raster.rio.to_raster(os.path.join(output_raster_path, file))
                            print("Save Done! :", file)


# TODO: fix this to keep year in the directory like hls: 6_interpolate/year/rm_id/crop_146/...tif
def main():
    # downloading modis data for year 2019 and days 150-190 same as hls and for tiles that SK is in there
    # test_webscrap_modis.main()
    # converting modis data from hdf to tif file
    # test_hdf2tif_modis.main()
    # Apply QA layer
    print("*********************************** Implementing QA Layer**********************************************")
    parameter_name = 'SCF_QC'
    value = 0
    sat_ids = ['MCD15A3H']
    tiles = ['h10v03', 'h10v04', 'h11v03', 'h11v04', 'h12v03']
    bands = ['Fpar_500m', 'FparLai_QC', 'Lai_500m']
    years = ['2019']
    days_num = np.arange(150, 207).tolist()
    modis_dir = os.path.join("data", "Modis_tif", "1_original")
    qc_path = os.path.join("data", "Modis_tif", "2_qa_masked")
    # modis_QA(modis_dir, qc_path, sat_ids, tiles, bands, years, days_num, parameter_name, value)

    print("*********************************** Implementing Merge **********************************************")
    qa_mask_path = os.path.join("data", "Modis_tif", "2_qa_masked")
    dir_path = os.path.join("data", "Modis_tif", "3_merged")
    days = list(map(str, days_num))
    # modis_merge(qa_mask_path, dir_path, sat_id, tiles, years, days)

    print("*********************************** Implementing Resample"
          ", Crop, Mask **********************************************")
    resolution = 30
    nan_values = -9999
    RM_IDs = [232, 322]
    # CROP_TYPE_NO = 153  # Spring wheat
    CROP_TYPE_NO = 146  # Canola

    modis_merged_path = os.path.join("data", "Modis_tif", "3_merged")
    modis_processed_path = os.path.join("data", "Modis_tif", "5_resampled_cropped_masked")
    # modis_mask(modis_merged_path, modis_processed_path, sat_ids, RM_IDs, CROP_TYPE_NO, years, days_num, resolution)

    print("*********************************** Interpolating **********************************************")
    save_path = os.path.join("data", "Modis_tif", "6_spatial_interpolated")
    bands = ['Fpar_500m', 'Lai_500m']
    hls_interpolate(modis_processed_path, save_path, RM_IDs, years, sat_ids, CROP_TYPE_NO, bands, nan_values)


if __name__ == '__main__':
    main()

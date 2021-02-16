import os
from utilities.Preprocessing import Preprocessing
import numpy as np
import rioxarray as rxr
from Utils.crop_soilgrid_smap import crop
from copy import copy
from rasterio.enums import Resampling
from utilities.FileUtils import FileUtils
import geopandas as gpd
from test_hls import hls_interpolate


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
    buffer = 300000
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
            # crop_type_path = os.path.join("data", "Modis_tif", "3_crop", year, str(rm_id))

            if not os.path.exists(crop_type_path):
                os.makedirs(crop_type_path)
            crop_type_image = "crop_mask_" + str(rm_id) + '.tif'
            # cropped_crop_type.rio.to_raster(os.path.join(crop_type_path, crop_type_image))
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


def main():
    p = Preprocessing()
    print(p)

    print("*********************************** Implementing Hdf2Tif**********************************************")
    hdf_path = os.path.join("data", "Modis_MCD43A4", "1_original")
    tif_path = os.path.join("data", "Modis_MCD43A4", "2_tif")
    years = ["2019"]
    tiles = ['h10v03', 'h10v04', 'h11v03', 'h11v04', 'h12v03']
    days_num = np.arange(151, 191).tolist()
    days = list(map(str, days_num))
    print(days)
    # for year in years:
    #     hdf_path_year = os.path.join(hdf_path, year)
    #     p.hdf2tif(hdf_path_year, tif_path, tiles, years, days)

    print("*********************************** Implementing QA Layer**********************************************")
    qa_masked_path = os.path.join("data", "Modis_MCD43A4", "3_qa")
    tile_path = os.path.join("data", "Modis_MCD43A4", "2_tif")
    sat_ids = ['MCD43A4']
    # p.apply_qa_mask_modis(tile_path, qa_masked_path, sat_ids, tiles, years, days)

    print("*********************************** Implementing Merge **********************************************")
    dir_region = os.path.join("data", "Modis_MCD43A4", "4_merge")
    sat_id = 'MCD43A4'
    region_name = '-'.join(tiles)
    # for year in years:
    #     for day in days:
    #         p.merge(qa_masked_path, dir_region, region_name, sat_id, tiles, year, day)

    print("*********************************** Implementing Resample"
          ", Crop, Mask **********************************************")
    origin_path = os.path.join("data", "Modis_MCD43A4", "4_merge")
    destination_path = os.path.join("data", "Modis_MCD43A4", "5_resampled_cropped_masked")
    RM_IDs = [232, 322]
    crop_type_id = 146
    # crop_type_id = 153
    resolution = 30
    # modis_mask(origin_path, destination_path, sat_ids, RM_IDs, crop_type_id, years, days_num, resolution)

    print("*********************************** Interpolating **********************************************")
    processed_path = os.path.join("data", "Modis_MCD43A4", "5_resampled_cropped_masked")
    save_path = os.path.join("data", "Modis_MCD43A4", "6_interpolated")

    bands_num = np.arange(1, 8).tolist()
    # there is no nan values in
    nan_values = -1000
    bands = list(map(str, bands_num))
    hls_interpolate(processed_path, save_path, RM_IDs, years, sat_ids, crop_type_id, bands, nan_values)


if __name__ == '__main__':
    main()

import os
import geopandas as gpd
import matplotlib.pyplot as plt
import rioxarray as rxr
import numpy as np
from Utils.crop_soilgrid_smap import crop
from Utils.crop_soilgrid_smap import plot_raster_polygon
from Utils.resample import resample
from rasterio.enums import Resampling
from copy import copy

CROP_TYPE_NO = 146  # Canola


# CROP_TYPE_NO = 153  # Sprint wheat


def main():
    # Select a RM
    RM_IDs = [232, 322]
    YEARS = [2019]
    # Pick a buffer and apply it to shapes: RM polygon
    buffer_size = 10000
    resolution = 30
    CRS = "epsg:32613"
    rm_polygons_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone",
                                    "polygon_EPSG_32613/cropped_by_CA_boundaries/", "polygon.shp")
    # Read AAFC data
    crop_type_path = os.path.join("AAFC_data", "aci_2019_sk_v1", "masked_cropID__146_153.tif")
    crop_type = rxr.open_rasterio(crop_type_path)

    for year in YEARS:

        daymet_path = os.path.join("data", "Daymet", "SK-tiles", "mosaic_final", str(year))

        for rm_id in RM_IDs:
            rm_polygons = gpd.read_file(rm_polygons_path)
            rm_polygon = rm_polygons.loc[rm_polygons["id"] == rm_id]
            rm_polygon_buff = rm_polygon.buffer(buffer_size)

            # Make a path to save processed files
            filtered_raster_path = os.path.join("data", "Daymet", "SK-tiles", "mosaic_final_processed", str(year),
                                                str(rm_id), "crop_" + str(CROP_TYPE_NO))
            if not os.path.exists(filtered_raster_path):
                os.makedirs(filtered_raster_path)

            # Crop AAFC
            cropped_crop_type = crop(crop_type, rm_polygon)
            cropped_crop_type.rio.to_raster(
                os.path.join("data", "Daymet", "SK-tiles", "mosaic_final_processed", str(year),
                             str(rm_id), "crop_mask_" + str(rm_id) + '.tif'))

            # Create mask based on crop type
            crop_type_np = cropped_crop_type.data[0]
            crop_type_np[crop_type_np != CROP_TYPE_NO] = 0
            crop_type_np[crop_type_np == CROP_TYPE_NO] = 1

            for root, dirs, files in os.walk(daymet_path):
                for file in files:
                    # Select daymet data
                    image_file_path = os.path.join(root, file)
                    daymet = rxr.open_rasterio(image_file_path)

                    # Crop each daymet data to the extent of RM with buffer
                    cropped_daymet_buff = crop(daymet, rm_polygon_buff)

                    # Upsample to 30m*30m
                    cropped_daymet_buff = cropped_daymet_buff.rio.reproject(CRS, resolution=resolution,
                                                                            resampling=Resampling.bilinear)
                    # Crop to the extent of RM without buffer
                    cropped_daymet = crop(cropped_daymet_buff, rm_polygon)

                    # Mask based on crop type
                    # TODO: Remove this later when the correct data was uploaded in the server
                    filtered_raster = copy(cropped_daymet)
                    common_x_len = min(np.shape(filtered_raster.data)[1], np.shape(cropped_crop_type.data)[1])
                    common_y_len = min(np.shape(filtered_raster.data)[2], np.shape(cropped_crop_type.data)[2])

                    # Crop all the bands
                    for data in filtered_raster.data:
                        data[0:common_x_len, 0:common_y_len] = data[0:common_x_len, 0:common_y_len] * crop_type_np[
                                                                                                      0:common_x_len,
                                                                                                      0:common_y_len]
                    filtered_raster.rio.to_raster(os.path.join(filtered_raster_path, file))
                    print("Save done!", os.path.join(filtered_raster_path, file))


if __name__ == '__main__':
    main()

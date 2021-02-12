import os
import geopandas as gpd
import rioxarray as rxr
import numpy as np
import matplotlib.pyplot as plt

YEAR = {"2019"}
TILE = ["T12UWC"]
RM_ID = [232, 322]
# CROP_TYPE = 146  # Canola


CROP_TYPE = 146  # Sprint wheat
target_resolution = 30


# TODO : change this script to a module and send it to utils folder

def main():
    print("hiiiiii")
    rm_polygons_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone",
                                    "polygon_EPSG_32613/cropped_by_CA_boundaries/", "polygon.shp")
    mgrs_rm_polygons_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone/sk-mgrs-rm", "sk-mgrs-rm.shp")

    rm_polygons = gpd.read_file(rm_polygons_path)
    mgrs_rm_polygon = gpd.read_file(mgrs_rm_polygons_path)

    # Pick either one of these paths:
    soilgrid_path = os.path.join("data", "Soil_Grid")
    # SMAP_path = os.path.join("data", "SMAP")
    PATH = soilgrid_path

    for rm_id in RM_ID:
        for year in YEAR:
            mgrs_tiles = mgrs_rm_polygon.loc[mgrs_rm_polygon['id'] == rm_id]['MGRS']
            if mgrs_tiles.empty:
                continue
            rm_polygon = rm_polygons.loc[rm_polygons['id'] == rm_id]
            tiles_list = ["T" + tile for tile in mgrs_tiles]
            tiles_path = '_'.join(tiles_list)
            rasters_path = os.path.join(PATH, year, tiles_path, str(rm_id))
            masked_rasters_path = os.path.join(rasters_path, "masked_crop_" + str(CROP_TYPE))
            if not os.path.exists(masked_rasters_path):
                os.makedirs(masked_rasters_path)

            crop_type_path = os.path.join(PATH, year, tiles_path, str(rm_id), "crop_mask_" + str(rm_id) + ".tif")
            print(crop_type_path)
            crop_type = rxr.open_rasterio(crop_type_path)
            crop_type_np = crop_type.data[0]
            crop_type_np[crop_type_np != CROP_TYPE] = 0
            crop_type_np[crop_type_np == CROP_TYPE] = 1
            print(np.unique(crop_type_np))

            for root, dirs, image_files in sorted(os.walk(rasters_path)):
                for image_file in image_files:
                    if '.tif' in image_file and not ("crop" in image_file):
                        image_path = os.path.join(rasters_path, image_file)
                        raster = rxr.open_rasterio(image_path)
                        raster_np = raster.data[0]
                        # the sizes might be different in 1 or 2 pixels
                        # RESIZE RASTER TO CROP_TYPE
                        print("shape of raster file", np.shape(raster_np))
                        print("shape of crop type file", np.shape(crop_type[0]))
                        raster_np = raster.data[0]
                        common_x_len = min(np.shape(raster_np)[0], np.shape(crop_type_np)[0])
                        common_y_len = min(np.shape(raster_np)[1], np.shape(crop_type_np)[1])
                        cropped_raster = raster
                        # the sizes might be different in 1 or 2 pixels
                        cropped_raster.data[0] = raster_np * crop_type_np
                        raster_masked_path = os.path.join(masked_rasters_path, image_file)
                        print(raster_masked_path)
                        cropped_raster.rio.to_raster(raster_masked_path)

                        # call mask


if __name__ == '__main__':
    main()

#
# if '.tif' in image_file and not ("crop" in image_file):
#     image_path = os.path.join(rasters_path, image_file)
#     raster = rxr.open_rasterio(image_path)
#     raster_np = raster.data[0]
#     common_x_len = min(np.shape(raster_np)[0], np.shape(crop_type_np)[0])
#     common_y_len = min(np.shape(raster_np)[1], np.shape(crop_type_np)[1])
#     cropped_raster = raster
#     # the sizes might be different in 1 or 2 pixels
#     cropped_raster.data[0] = raster.data[0][0:common_x_len, 0:common_y_len] * crop_type_np[
#                                                                               0:common_x_len,
#                                                                               0:common_y_len]
#     raster_masked_path = os.path.join(rasters_masked_path, image_file)
#     print(raster_masked_path)
#     cropped_raster.rio.to_raster(raster_masked_path)

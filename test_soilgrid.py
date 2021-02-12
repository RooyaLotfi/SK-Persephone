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


def main():
    print("*********************************** Implementing cropping to the "
          "extent of RM with buffer**********************************************")
    CRS = "EPSG:32613"
    buffer = 10000
    RM_IDs = [232, 322]
    resolution = 30
    CROP_TYPE_NO = 153
    CROP_TYPE_NO = 146
    features = ['bdod', 'cec', 'cfvo', 'clay', 'phh2o', 'sand', 'silt', 'soc']
    layers = ['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']
    # features = ['soc', 'silt']
    # layers =['0-5', '5-15', '15-30', '30-60', '60-100', '100-200']
    smap_path = os.path.join("data", "soilgrid", "1_original")
    smap_processed_path = os.path.join("data", "soilgrid", "5_resampled_cropped_masked")
    f = FileUtils()
    rm_polygons_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone",
                                    "polygon_EPSG_32613/cropped_by_CA_boundaries/", "polygon.shp")
    rm_polygons = gpd.read_file(rm_polygons_path)

    # Read AAFC data
    crop_type_path = os.path.join("AAFC_data", "aci_2019_sk_v1", "masked_cropID__146_153.tif")
    crop_type = rxr.open_rasterio(crop_type_path)

    for rm_id in RM_IDs:
        rm_polygon = rm_polygons.loc[rm_polygons['id'] == rm_id]
        output_raster_path = os.path.join(smap_processed_path, str(rm_id), "crop_" + str(CROP_TYPE_NO))
        if not os.path.exists(output_raster_path):
            os.makedirs(output_raster_path)
        rm_polygon_buff = rm_polygon.buffer(buffer)

        # crop by polygon
        cropped_crop_type = crop(crop_type, rm_polygon)
        crop_type_path = os.path.join("data", "soilgrid", "3_crop", str(rm_id))
        if not os.path.exists(crop_type_path):
            os.makedirs(crop_type_path)
        crop_type_image = "crop_mask_" + str(rm_id) + '.tif'
        cropped_crop_type.rio.to_raster(os.path.join(crop_type_path, crop_type_image))
        print("Saved crop masks done!", os.path.join(crop_type_path, crop_type_image))
        # Create mask based on crop type
        crop_type_np = cropped_crop_type.data[0]
        crop_type_np[crop_type_np != CROP_TYPE_NO] = 0
        crop_type_np[crop_type_np == CROP_TYPE_NO] = 1

        for root, dir, files in os.walk(smap_path):
            for file in files:
                feature, layer, arg3, province = f.parse_soilgrid_file(file)
                if feature in features and layer in layers:
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

                    temp_raster.rio.to_raster(os.path.join(output_raster_path, file))
                    print("Save Done! :", file)


if __name__ == '__main__':
    main()

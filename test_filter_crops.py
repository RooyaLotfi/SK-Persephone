import os
import geopandas as gpd
import rioxarray as rxr
import numpy as np
import matplotlib.pyplot as plt

SATTELLITE = ["S30", "L30"]
YEAR = ["2019"]
TILE = ["T12UWC"]
RM_ID = [322]
# CROP_TYPE = 146  # Canola


CROP_TYPE = 146  # Sprint wheat
# TODO : change this script to a module and send it to utils folder
# TODO: Write a function called get tiles which gets an ID and returns all the tiles in it
def main():
    for satellite in SATTELLITE:
        for year in YEAR:
            for tile in TILE:
                for rm_id in RM_ID:
                    # TODO: fix tiles part to get it based on ID
                    crop_type_path = os.path.join("data", "2_qa_masked", satellite, year, tile, str(rm_id), "cropped",
                                                  "crop_mask_" + str(rm_id) + ".tif")
                    crop_type = rxr.open_rasterio(crop_type_path)
                    crop_type_np = crop_type.data[0]
                    crop_type_np[crop_type_np != CROP_TYPE] = 0
                    crop_type_np[crop_type_np == CROP_TYPE] = 1
                    print(np.unique(crop_type_np))
                    rasters_path = os.path.join("data", "2_qa_masked", satellite, year, tile, str(rm_id), "cropped")
                    rasters_masked_path = os.path.join("data", "2_qa_masked", satellite, year, tile, str(rm_id),
                                                       "hls_masked_crop_" + str(CROP_TYPE))
                    if not os.path.exists(rasters_masked_path):
                        os.makedirs(rasters_masked_path)
                    for root, dirs, image_files in sorted(os.walk(rasters_path)):
                        for image_file in image_files:
                            # print(image_file)
                            if '.tif' in image_file and not ("crop" in image_file):
                                image_path = os.path.join(rasters_path, image_file)
                                raster = rxr.open_rasterio(image_path)
                                raster_np = raster.data[0]
                                # TODO: change bellow line to copy array by reference
                                cropped_raster = raster
                                # the sizes might be different in 1 or 2 pixels
                                cropped_raster.data[0] = raster.data[0] * crop_type_np
                                raster_masked_path = os.path.join(rasters_masked_path, image_file)
                                print(raster_masked_path)
                                cropped_raster.rio.to_raster(raster_masked_path)
                    # print(rasters_path)

    # crop_type_path = os.path.join("data", "2_qa_masked", "S30", "2019", "T12UWB", "232", "crop_mask_232.tif")
    # raster_path = os.path.join("data", "2_qa_masked", "S30", "2019", "T12UWB", "232", "S30_T12UWB_2019_151_01.tif")
    # raster_masked_path = os.path.join("data", "2_qa_masked", "S30", "2019", "T12UWB", "232", "masked",
    #                                   "S30_T12UWB_2019_151_01.tif")
    # crop_type = rxr.open_rasterio(crop_type_path)
    # raster = rxr.open_rasterio(raster_path)
    # raster_crop = raster
    # print(crop_type.rio)
    # crop_type_np = crop_type.data[0]
    # raster_np = raster_crop.data[0]
    # print(raster_crop.rio.resolution())
    # print(crop_type.rio.resolution())
    # print(np.shape(crop_type_np))
    # print(np.shape(raster_np))
    # common_x_len = min(np.shape(raster_np)[0], np.shape(crop_type_np)[0])
    # common_y_len = min(np.shape(raster_np)[1], np.shape(crop_type_np)[1])
    # print(" common len x ", common_x_len, " common len y", common_y_len)
    #
    # crop_type_np = crop_type_np[0:common_x_len, 0:common_y_len]
    # raster_np = raster_np[0:common_x_len, 0:common_y_len]
    # print(np.shape(crop_type_np))
    # print(np.shape(raster_np))
    #
    # raster_cropped = raster_crop
    # raster_cropped.data[0] = raster_crop.data[0][0:common_x_len, 0:common_y_len] * crop_type_np
    #
    # print(raster_cropped)
    # raster_cropped.rio.to_raster(raster_masked_path)
    # fig, ax = plt.subplots(figsize=(3, 5))
    # raster_cropped.plot()
    # plt.show()
    # #
    # # crop_146 = np.copy(crop_type_np)
    # # crop_146[crop_146 != 146] = 0
    # # crop_146[crop_146 == 146] = 1
    # # crop_153 = np.copy(crop_type_np)
    # # crop_153[crop_153 != 153] = 0
    # # crop_153[crop_153 == 153] = 1
    # # print(np.unique(crop_146))
    # # print(np.unique(crop_153))
    # #
    # # print(crop_146*raster_np)
    # # print(crop_153*raster_np)
    # # raster_146 = crop_146*raster_np
    # # raster_153 = crop_153*raster_np
    # # x = raster_146[:, 0]
    # # print(np.shape(x))
    # # y = raster_146[:, 1]
    # # print(np.shape(y))
    # #
    # # plt.imshow(raster_146)
    # # plt.show()
    # #
    # # plt.imshow(raster_153)
    # # plt.show()
    # # plt.imshow(crop_153)
    # #
    # # plt.show()
    # # plt.imshow(crop_146)
    # # plt.show()


if __name__ == '__main__':
    main()

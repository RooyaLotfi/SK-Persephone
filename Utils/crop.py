import re
import rioxarray
import os
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import mapping
# from utilities.FileUtils import FileUtils
import numpy as np

ID_num = np.arange(1, 622).tolist()

RM_IDs = [322, 261, 292]
RM_IDs = ID_num
# Default values for satellite, tile, extension of files to download, and year
SATELLITE = {'L30', 'S30'}
YEAR = {'2019'}


def crop(raster, polygon):
    if raster.rio.crs != polygon.crs:
        print("Raster CRS: ", raster.rio.crs)
        raster = raster.rio.reproject(polygon.crs)
        print("Is changed to: ", polygon.crs)
    cropped_raster = raster.rio.clip(polygon.geometry.apply(mapping))
    # Write the data to a new geotiff file
    return cropped_raster


def plot_raster_polygon(raster, polygon):
    f, ax = plt.subplots(figsize=(10, 5))
    raster[0].plot.imshow(ax=ax)

    polygon.plot(ax=ax, alpha=.5)
    ax.set(title="Raster Layer with Shapefile Overlayed")

    ax.set_axis_off()
    plt.show()


def crop_RM(qa_mask_path, rm_polygons_path, mgrs_rm_polygons_path, crop_type_path):
    rm_polygons = gpd.read_file(rm_polygons_path)
    mgrs_rm_polygon = gpd.read_file(mgrs_rm_polygons_path)
    crop_type = rioxarray.open_rasterio(crop_type_path)

    for rm_id in RM_IDs:
        mgrs_tiles = mgrs_rm_polygon.loc[mgrs_rm_polygon['id'] == rm_id]['MGRS']
        if mgrs_tiles.empty:
            continue
        # print("ID: ", rm_id, " is in MGRS tiles: \n", mgrs_tiles.to_list())
        rm_polygon = rm_polygons.loc[rm_polygons['id'] == rm_id]
        for satellite in SATELLITE:
            for year in YEAR:
                tiles_list = ["T" + tile for tile in mgrs_tiles]
                tiles_path = '_'.join(tiles_list)
                rasters_path = os.path.join(qa_mask_path,
                                            satellite, year, tiles_path)
                if not os.path.exists(rasters_path):
                    os.makedirs(rasters_path)
                print("RASTER ", rasters_path)
                id_path = os.path.join(rasters_path, str(rm_id))
                if not os.path.exists(id_path):
                    os.makedirs(id_path)

                cropped_crop_type_path = os.path.join(id_path, "crop_mask_" + str(rm_id) + ".tif")
                print("AAFC ", cropped_crop_type_path)
                cropped_crop_type = crop(crop_type, rm_polygon)
                cropped_crop_type.rio.to_raster(cropped_crop_type_path)
                for root, dirs, image_files in sorted(os.walk(rasters_path)):
                    for image_file in image_files:
                        if '.tif' in image_file:
                            image_path = os.path.join(rasters_path, image_file)
                            raster = rioxarray.open_rasterio(image_path)
                            cropped_raster = crop(raster, rm_polygon)
                            cropped_raster.rio.to_raster(os.path.join(id_path, image_file))
                            print("CROPPED IMAGE", os.path.join(id_path, image_file))
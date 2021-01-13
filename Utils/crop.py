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


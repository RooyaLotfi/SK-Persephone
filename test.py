# last updated: 2020/12/30
import datetime  # For date_to_day_of_year
import os  # For reproject_to_file
import pathlib
from typing import Union  # For reproject_to_file
import matplotlib.pyplot as plt
import geopandas as gpd
import numpy as np  # For reproject_tif, reproject_to_file
import numpy.ma as ma  # For reproject_tif ; For mask_data
import rasterio  # For reproject_tif; For reproject_to_file; For mask_data
import scipy.interpolate  # For interpolate
from osgeo import gdal, ogr  # For overlap_detection
# from rasterio.crs import CRS  # For reproject_tif
from rasterio.mask import mask  # For mask_data
from rasterio.warp import \
    calculate_default_transform  # For reproject_to_file; For reproject_tif
from rasterio.warp import Resampling, reproject
from skimage.transform import \
    resize as skimage_resample  # For the resample portion of mask_data
import earthpy as et
import earthpy.plot as ep
import earthpy.spatial as es
import rasterio as rio
from rasterio.plot import show


SATELLITE = 'L30'
# First three characters of the tiles
TILE = '12UXA'
YEAR = '2019'
CRS = "epsg:32613"


def mask_data(tif_dir, tif_out, shapefile, crop_data=False, all_touched=True, resample=None, fill_value=None):
    """
    Masks the data from the Geotiff file based on an input shapefile.
    Args:
        tif_dir: Directory of the Geotiff file.
        shapefile: Geopandas geodataframe with a 'geometry' column.
        crop_data: If set to true, crops the data to the size/shape of the shapefile.
        all_touched: Include a pixel in the mask if it touches any of the shapes of the shapefile
        resample: Resample the array to the specified tuple
        fill_value: Set the fill value for the masks

    Returns:
        out_image_ma -- The masked array representation of the Geotiff file.
        out_info -- Information of the 1_original Geotiff file.
    """
    with rasterio.open(tif_dir) as src:
        out_image, out_transform = mask(src, shapefile['geometry'], crop=crop_data, all_touched=all_touched)
        # _ can be replaced as out_transform and be used for mapping pixel coordinates of out_image to another coordinate system.
        out_meta = src.meta
    out_image_ma = np.squeeze(ma.array(out_image.astype('float')))

    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})

    with rasterio.open(tif_out, "w", **out_meta) as dest:
        dest.write(out_image)

    #return out_image_ma, out_transform

rm_tile_path = os.path.join("sk-mgrs-rm", "sk-mgrs-rm.shp")
rm_tiles = gpd.read_file(rm_tile_path)
print("location of item 0 : ", rm_tiles.iloc[[1]])
rm_tiles = rm_tiles.loc[rm_tiles['MGRS'] == TILE]
print(rm_tiles.MGRS)
rm_tiles = rm_tiles.iloc[[0]]
rm_tiles.to_file("sample.shp")

print(rm_tiles)
tile_raster = os.path.join("data", "hls_tmp_tif", TILE, YEAR, SATELLITE, "HLS.L30.T12UXA.2019001.v1.4.09.tif.epsg.32613.tif")
tile_raster_out = os.path.join("data", "hls_tmp_tif", TILE, YEAR, SATELLITE, "cropped", "HLS.L30.T12UXA.2019001.v1.4.092.tif")

mask_data(tile_raster, tile_raster_out, rm_tiles)

with rio.open(tile_raster) as src:
    lidar_chm_im = src.read(masked=True)[0]
    extent = rio.plot.plotting_extent(src)
    soap_profile = src.profile
print("_______________", src.crs)


fig, ax = plt.subplots(figsize=(10, 10))
ep.plot_bands(lidar_chm_im,
              cmap='terrain',
              extent=extent,
              ax=ax,
              cbar=False)
rm_tiles.plot(ax=ax, alpha=.6, color='g')
plt.show()

crop_type_path = os.path.join("AAFC_data", "aci_2019_sk_v1", "masked_cropID__146_153.tif")
crop_type_path_cropped = os.path.join("data", "hls_tmp_tif", TILE, YEAR, SATELLITE, "cropped", "masked_cropID__146_153_masked.tif")
#crop_mask = rio.open(crop_type_path)

mask_data(crop_type_path, crop_type_path_cropped, rm_tiles)

with rio.open(crop_type_path) as src:
    lidar_chm_im = src.read(masked=True)[0]
    extent = rio.plot.plotting_extent(src)
    soap_profile = src.profile
print("_______________", src.crs)

fig, ax = plt.subplots(figsize=(10, 10))
ep.plot_bands(lidar_chm_im,
              cmap='terrain',
              extent=extent,
              ax=ax,
              cbar=False)
rm_tiles.plot(ax=ax, alpha=.6, color='g')
plt.show()


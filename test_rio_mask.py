import os
import fiona
import rasterio
import rasterio.mask
import rasterio as rio
# Default values for satellite, tile, extension of files to download, and year
SATELLITE = 'L30'
# First three characters of the tiles
TILE = '12UXA'
YEAR = '2019'
CRS = "epsg:32613"



# last updated: 2020/12/30
import datetime  # For date_to_day_of_year
import os  # For reproject_to_file
import pathlib
from typing import Union  # For reproject_to_file

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

# from osgeo.utils.gdal_edit import * # For whatever reason it does not exist in the GDAL folder when installed through Conda forge. For now I
# copied the python script into a separate .py called gdal_edit and will import it as is


def mask_data(tif_dir, shapefile, crop_data=False, all_touched=True, resample=None, fill_value=None):
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
        out_image, _ = mask(src, shapefile['geometry'], crop=crop_data, all_touched=all_touched)
        # _ can be replaced as out_transform and be used for mapping pixel coordinates of out_image to another coordinate system.
        out_info = src
    out_image_ma = np.squeeze(ma.array(out_image.astype('float')))

    if resample is not None:
        dst_width, dst_height = resample
        out_image_ma = skimage_resample(out_image_ma, (dst_height, dst_width), order=0)
        out_image_ma = np.ma.masked_array(out_image_ma)
        ma.set_fill_value(out_image_ma, src.meta['nodata'])
        out_image_ma = ma.masked_values(out_image_ma, src.meta['nodata'])

    if out_info.meta['nodata'] is not None:
        if np.issubdtype(out_image_ma.dtype, np.floating):
            out_image_ma = ma.masked_values(out_image_ma, out_info.meta[
                'nodata'])  # Better floating point approximation than np.where
        else:
            out_image_ma = ma.masked_equal(out_image_ma, out_info.meta['nodata'])

    if fill_value is not None:  # Change the fill value to the one specified in the argument
        ma.set_fill_value(out_image_ma, fill_value)
        src.meta[
            'nodata'] = fill_value  # Not relevant for now because the function is using out_info, but letting it here in case there's a change of functionality

    # out_image_ma = ma.masked_values(out_image_ma, src.meta['nodata'])
    return out_image_ma, out_info

rm_tile_path = os.path.join("sk-mgrs-rm", "sk-mgrs-rm.shp")
rm_tiles = gpd.read_file(rm_tile_path)

tile_raster = os.path.join("data", "hls_tmp_tif", TILE, YEAR, SATELLITE, "HLS.L30.T12UXA.2019001.v1.4.01.tif")
tile_raster_out = os.path.join("data", "hls_tmp_tif", TILE, YEAR, SATELLITE, "cropped", "HLS.L30.T12UXA.2019001.v1.4.01.tif")

out_image_ma, out_info = mask_data(tile_raster, rm_tiles)

print(out_image_ma)
print(out_info)
#!/usr/local/bin/python3

import os
import numpy as np
import matplotlib.pyplot as plt
import subprocess, glob
from osgeo import gdal, gdalnumeric, ogr
import pandas as pd
import csv
import ast
import pathlib
import timeit
import rasterio
import sys
import shutil
from utilities.FileUtils import FileUtils
from shapely.geometry import box
import geopandas as gpd
from fiona.crs import from_epsg
import json
import pycrs
from rasterio.mask import mask
import re
from rasterio.features import shapes
import fiona
import datetime
from shapely.geometry import Polygon, Point
from geojson import Feature, FeatureCollection, dump
from scipy import stats
import numpy.ma as ma


class Preprocessing:
    def __init__(self):
        pass

    def hdf2tif(self, hdf_path, tif_path, tiles, years, days):
        # Extracts tif subdatasets from HDF files stored in hdf_path directory
        # Extracted tif files are stored in tif_path directory
        if not os.path.exists(tif_path):
            os.mkdir(tif_path)
            print("Directory ", tif_path, "created ")
        else:
            print("Directory ", tif_path, "already exists")

        file_lists = []
        file_util = FileUtils()
        for root, dirs, files in sorted(os.walk(hdf_path)):
            files.sort(key=file_util.natural_keys)
            file_lists.append(files)

        for files in file_lists:
            for file in files:
                if ".hdf" in file:
                    sat_id, tile, year, day, ext = file_util.parse_tile_file(file)

                    if tile in tiles and year in years and day in days:
                        path = f"{hdf_path}/{file}"
                        hdf_file = gdal.Open(path, gdal.GA_ReadOnly)
                        subdatasets = hdf_file.GetSubDatasets()

                        for i in range(0, len(subdatasets)):
                            subdataset_name = subdatasets[i][0]

                            if 'Quality' in subdataset_name:
                                band_name = 'Q' + subdataset_name[-1]
                                data_type = gdal.GDT_Byte
                            else:
                                band_name = subdataset_name[-1]
                                data_type = gdal.GDT_Int16

                            band = gdal.Open(subdataset_name, gdal.GA_ReadOnly)
                            band_image = band.ReadAsArray()

                            path_out = f"{tif_path}/{file[:-4]}.{band_name}.tif"
                            tif_file = gdal.GetDriverByName('GTiff').Create(path_out, band.RasterXSize,
                                                                            band.RasterYSize, 1, data_type)

                            tif_file.SetGeoTransform(band.GetGeoTransform())
                            tif_file.SetProjection(band.GetProjection())
                            tif_file.GetRasterBand(1).WriteArray(band_image)
                            print("Save Done: %s" % path_out)
                            tif_file.FlushCache()
                            del tif_file

    def reproject(self, tile_path, reprojected_path, tiles, years, days):
        # Reprojects tiles stored in tile_path directory
        # Results are stored in reprojected_path directory
        if not os.path.exists(reprojected_path):
            os.mkdir(reprojected_path)
            print("Directory ", reprojected_path, "created ")
        else:
            print("Directory ", reprojected_path, "already exists")

        file_lists = []
        file_util = FileUtils()
        for root, dirs, files in sorted(os.walk(tile_path)):
            files.sort(key=file_util.natural_keys)
            file_lists.append(files)

        for files in file_lists:
            for file in files:
                if ".tif" in file:
                    sat_id, tile, year, day, band_no = file_util.parse_tile_file(file)

                    if tile in tiles and year in years and day in days:
                        path_in = f"{tile_path}/{file}"
                        image = gdal.Open(path_in, gdal.GA_ReadOnly)
                        proj_in = image.GetProjection()
                        path_out = f"{reprojected_path}/{file}"
                        command = f"gdalwarp -s_srs {proj_in} -t_srs '+proj=utm +zone=14 +datum=WGS84 +units=m +no_defs' -of GTiff -overwrite {path_in} {path_out}"
                        output = os.system(command)
                        print(output)

    def clip_by_extent(self, qa_masked_path, clipped_path, extent_layer_path, clipped_tile_name, tiles, years, days):
        # Clips tiles stored in qa_masked_path directory using a layer extent
        # Results are stored in clipped_path directory
        # extent_layer_path: Path to extent layer
        # clipped_tile_name: String representing clipped tile's location
        if not os.path.exists(clipped_path):
            os.mkdir(clipped_path)
            print("Directory ", clipped_path, "created ")
        else:
            print("Directory ", clipped_path, "already exists")

        extent_layer = gdal.Open(extent_layer_path, gdal.GA_ReadOnly)
        geo_trans = extent_layer.GetGeoTransform()

        minx = geo_trans[0]
        maxy = geo_trans[3]
        maxx = minx + geo_trans[1] * extent_layer.RasterXSize
        miny = maxy + geo_trans[5] * extent_layer.RasterYSize
        bounding_box = [minx, maxy, maxx, miny]

        file_lists = []
        file_util = FileUtils()
        for root, dirs, files in sorted(os.walk(qa_masked_path)):
            files.sort(key=file_util.natural_keys)
            file_lists.append(files)

        for files in file_lists:
            for file in files:
                if ".tif" in file:
                    sat_id, tile, year, day, band_no = file_util.parse_qamasked_file(file)

                    if tile in tiles and year in years and day in days:
                        dir_sat = clipped_path + '/' + sat_id
                        if not os.path.exists(dir_sat):
                            os.mkdir(dir_sat)
                            print("Directory ", dir_sat, " Created ")

                        dir_year = dir_sat + '/' + year
                        if not os.path.exists(dir_year):
                            os.mkdir(dir_year)
                            print("Directory ", dir_year, " Created ")

                        dir_tile = dir_year + '/' + clipped_tile_name
                        if not os.path.exists(dir_tile):
                            os.mkdir(dir_tile)
                            print("Directory ", dir_tile, " Created ")

                        path = f"{qa_masked_path}/{sat_id}/{year}/{tile}/{file}"
                        image = gdal.Open(path, gdal.GA_ReadOnly)

                        file_name_out = file.replace(tile, clipped_tile_name)
                        path_out = f"{clipped_path}/{sat_id}/{year}/{clipped_tile_name}/{file_name_out}"

                        image = gdal.Translate(path_out, image, projWin=bounding_box)
                        image = None
                        print("Save Done: %s" % path_out)

    def resample(self, qa_masked_path, resampled_path, target_resolution, resampling_method, sat_ids, tiles, years,
                 days):
        # Resamples tiles stored in qa_masked_path directory into target_resolution using a resampling_method
        # Results are stored in resampled_path directory
        if not os.path.exists(resampled_path):
            os.mkdir(resampled_path)
            print("Directory ", resampled_path, "created ")
        else:
            print("Directory ", resampled_path, "already exists")

        file_lists = []
        file_util = FileUtils()
        for root, dirs, files in sorted(os.walk(qa_masked_path)):
            files.sort(key=file_util.natural_keys)
            file_lists.append(files)

        for files in file_lists:
            for file in files:
                if ".tif" in file:
                    sat_id, tile, year, day, band_no = file_util.parse_qamasked_file(file)

                    if sat_id in sat_ids and tile in tiles and year in years and day in days:
                        dir_sat = resampled_path + '/' + sat_id
                        if not os.path.exists(dir_sat):
                            os.mkdir(dir_sat)
                            print("Directory ", dir_sat, " Created ")

                        dir_year = dir_sat + '/' + year
                        if not os.path.exists(dir_year):
                            os.mkdir(dir_year)
                            print("Directory ", dir_year, " Created ")

                        dir_tile = dir_year + '/' + tile
                        if not os.path.exists(dir_tile):
                            os.mkdir(dir_tile)
                            print("Directory ", dir_tile, " Created ")

                        path_in = f"{qa_masked_path}/{sat_id}/{year}/{tile}/{file}"
                        path_out = f"{resampled_path}/{sat_id}/{year}/{tile}/{file}"

                        command = f"gdalwarp -s_srs EPSG:32614 -t_srs EPSG:32614 -tr {target_resolution}.0 {target_resolution}.0 -r {resampling_method} -tap -of GTiff {path_in} {path_out}"
                        output = os.system(command)
                        print(output)

    def resize(self, qa_masked_path, resized_path, target_size, tiles, years, days):
        # Resizes tiles stored in qa_masked_path directory into target_size x target_size
        # Results are stored in resized_path directory
        if not os.path.exists(resized_path):
            os.mkdir(resized_path)
            print("Directory ", resized_path, "created ")
        else:
            print("Directory ", resized_path, "already exists")

        file_lists = []
        file_util = FileUtils()
        for root, dirs, files in sorted(os.walk(qa_masked_path)):
            files.sort(key=file_util.natural_keys)
            file_lists.append(files)

        for files in file_lists:
            for file in files:
                if ".tif" in file:
                    sat_id, tile, year, day, band_no = file_util.parse_qamasked_file(file)

                    if tile in tiles and year in years and day in days:
                        dir_sat = resized_path + '/' + sat_id
                        if not os.path.exists(dir_sat):
                            os.mkdir(dir_sat)
                            print("Directory ", dir_sat, " Created ")

                        dir_year = dir_sat + '/' + year
                        if not os.path.exists(dir_year):
                            os.mkdir(dir_year)
                            print("Directory ", dir_year, " Created ")

                        dir_tile = dir_year + '/' + tile
                        if not os.path.exists(dir_tile):
                            os.mkdir(dir_tile)
                            print("Directory ", dir_tile, " Created ")

                        path_in = f"{qa_masked_path}/{sat_id}/{year}/{tile}/{file}"
                        image = gdal.Open(path_in, gdal.GA_ReadOnly)
                        image_val = image.ReadAsArray()
                        image_val = image_val[:target_size, :target_size]

                        path_out = f"{resized_path}/{sat_id}/{year}/{tile}/{file}"
                        out_file = gdal.GetDriverByName('GTiff').Create(path_out, target_size, target_size, 1,
                                                                        image.GetRasterBand(1).DataType)

                        out_file.SetGeoTransform(image.GetGeoTransform())
                        out_file.SetProjection(image.GetProjection())
                        out_file.GetRasterBand(1).WriteArray(image_val)
                        print("Save Done: %s" % path_out)
                        out_file.FlushCache()
                        del out_file

    def cirrus(self, pixel):
        # check for cirrus
        pix = int(pixel / (2 ** 0))
        bit = pix - (int(pix / 2) * 2)
        return bit

    def cloud(self, pixel):
        # check for clouds
        pix = int(pixel / (2 ** 1))
        bit = pix - (int(pix / 2) * 2)
        return bit

    def adjclouds(self, pixel):
        # check for adjacent clouds
        pix = int(pixel / (2 ** 2))
        bit = pix - (int(pix / 2) * 2)
        return bit

    def clouds_shadow(self, pixel):
        # check for cloud shadows
        pix = int(pixel / (2 ** 3))
        bit = pix - (int(pix / 2) * 2)
        return bit

    def water(self, pixel):
        # check for water
        pix = int(pixel / (2 ** 5))
        bit = pix - (int(pix / 2) * 2)
        return bit

    def create_qa_mask_hls(self, qa_layer, return_filled=True):
        # Creates QA mask for clouds, adjacent clouds, cloud shadows, water, and cirrus from HLS QA layer
        qa_layer_ma = ma.array(qa_layer.astype('float'))
        for pixel in range(0, 256):
            if pixel in qa_layer:
                cloud_found = self.cloud(pixel)
                adjcloud_found = self.adjclouds(pixel)
                clouds_shadow_found = self.clouds_shadow(pixel)
                water_found = self.water(pixel)
                cirrus_found = self.cirrus(pixel)

                if (cloud_found or adjcloud_found or clouds_shadow_found or water_found or cirrus_found) == 1:
                    qa_layer_ma = ma.masked_equal(qa_layer_ma, pixel)
                else:
                    qa_layer_ma[qa_layer == pixel] = 1.0
        ma.set_fill_value(qa_layer_ma, np.nan)
        if return_filled:
            qa_layer_ma = qa_layer_ma.filled()
        return qa_layer_ma

    def create_qa_mask_modis(self, qa_layer):
        # Creates QA mask for fill-values from MODIS QA layer
        qa_layer[qa_layer == 1] = 0
        qa_layer[qa_layer == 255] = 1
        qa_layer = 1 - qa_layer
        return qa_layer

    def apply_qa_mask_hls(self, tile_path, qa_masked_path, sat_ids, tiles, years, days):
        # Applies QA masks on each HLS band
        if not os.path.exists(qa_masked_path):
            os.mkdir(qa_masked_path)
            print("Directory ", qa_masked_path, "created ")
        else:
            print("Directory ", qa_masked_path, "already exists")

        file_lists = []
        file_util = FileUtils()
        for root, dirs, files in sorted(os.walk(tile_path)):
            files.sort(key=file_util.natural_keys)
            file_lists.append(files)

        for files in file_lists:
            for image_path in files:
                if ".tif" in image_path:
                    sat_id, tile, year, day, band_no = file_util.parse_tile_file(image_path)

                    if sat_id in sat_ids and tile in tiles and year in years and day in days:
                        print(sat_id, tile, year, day, band_no)

                        dir_sat = qa_masked_path + '/' + sat_id
                        if not os.path.exists(dir_sat):
                            os.mkdir(dir_sat)
                            print("Directory ", dir_sat, " Created ")

                        dir_year = dir_sat + '/' + str(year)
                        if not os.path.exists(dir_year):
                            os.mkdir(dir_year)
                            print("Directory ", dir_year, " Created ")

                        dir_tile = dir_year + '/' + tile
                        if not os.path.exists(dir_tile):
                            os.mkdir(dir_tile)
                            print("Directory ", dir_tile, " Created ")

                        if sat_id == 'S30':
                            qa_layer_no = '14'
                        elif sat_id == 'L30':
                            qa_layer_no = '11'

                        if band_no != qa_layer_no:
                            qa_layer_path = "{}/HLS.{}.{}.{}{}.v1.4.hdf_sds_{}.tif".format(tile_path, sat_id, tile,
                                                                                           year, day, qa_layer_no)
                            qa_layer = np.array(gdal.Open(qa_layer_path).ReadAsArray(), dtype='uint8')
                            qa_mask = self.create_qa_mask_hls(qa_layer)

                            image = rasterio.open(tile_path + '/' + image_path)
                            tif_file = image.read(1)
                            qa_masked = tif_file * qa_mask

                            rastout = "{}/{}/{}/{}/{}_{}_{}_{}_{}.tif".format(qa_masked_path, sat_id, year, tile,
                                                                              sat_id, tile, year, day, band_no)
                            kwargs = image.meta
                            kwargs.update(dtype=rasterio.float32)
                            with rasterio.open(rastout, 'w', **kwargs) as dst:
                                dst.write(qa_masked.astype(rasterio.float32), 1)
                                print("Save Done: %s" % rastout)

    def apply_qa_mask_modis(self, tile_path, qa_masked_path, sat_ids, tiles, years, days):
        # Applies QA Masks on MODIS bands
        if not os.path.exists(qa_masked_path):
            os.mkdir(qa_masked_path)
            print("Directory ", qa_masked_path, "created ")
        else:
            print("Directory ", qa_masked_path, "already exists")

        file_lists = []
        file_util = FileUtils()
        for root, dirs, files in sorted(os.walk(tile_path)):
            files.sort(key=file_util.natural_keys)
            file_lists.append(files)

        for files in file_lists:
            for file in files:
                if ".tif" in file:
                    sat_id, tile, year, day, band_no = file_util.parse_tile_file(file)

                    if sat_id in sat_ids and tile in tiles and year in years and day in days:
                        print(sat_id, tile, year, day, band_no)

                        dir_sat = qa_masked_path + '/' + sat_id
                        if not os.path.exists(dir_sat):
                            os.mkdir(dir_sat)
                            print("Directory ", dir_sat, " Created ")

                        dir_year = dir_sat + '/' + year
                        if not os.path.exists(dir_year):
                            os.mkdir(dir_year)
                            print("Directory ", dir_year, " Created ")

                        dir_tile = dir_year + '/' + tile
                        if not os.path.exists(dir_tile):
                            os.mkdir(dir_tile)
                            print("Directory ", dir_tile, " Created ")

                        if sat_id == "MCD43A4":
                            qa_layers = ['Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7']

                        if band_no not in qa_layers:
                            qa_file_name = file.replace(f".{band_no}.", f".Q{band_no}.")
                            qa_layer_path = f"{tile_path}/{qa_file_name}"
                            qa_layer = np.array(gdal.Open(qa_layer_path).ReadAsArray(), dtype='uint8')
                            qa_mask = self.create_qa_mask_modis(qa_layer)

                            image = rasterio.open(f"{tile_path}/{file}")
                            tif_file = image.read(1)
                            qa_masked = tif_file * qa_mask

                            rastout = f"{qa_masked_path}/{sat_id}/{year}/{tile}/{sat_id}_{tile}_{year}_{day}_{band_no}.tif"
                            kwargs = image.meta
                            kwargs.update(dtype=rasterio.int16)
                            with rasterio.open(rastout, 'w', **kwargs) as dst:
                                dst.write(qa_masked.astype(rasterio.int16), 1)
                                print("Save Done: %s" % rastout)

    def apply_crop_masks(self, qa_masked_path, crop_mask_path, preprocessed_path, tiles, years, days, classes):
        # Applies crop masks on the QA masked output
        if not os.path.exists(preprocessed_path):
            os.mkdir(preprocessed_path)
            print("Directory ", preprocessed_path, "created ")
        else:
            print("Directory ", preprocessed_path, "already exists")

        file_lists = []
        file_util = FileUtils()
        for root, dirs, files in sorted(os.walk(qa_masked_path)):
            files.sort(key=file_util.natural_keys)
            file_lists.append(files)

        for files in file_lists:
            for image_path in files:
                if ".tif" in image_path:
                    sat_id, tile, year, day, band_no = file_util.parse_qamasked_file(image_path)

                    if tile in tiles and year in years and day in days:
                        print(sat_id, tile, year, day, band_no)

                        dir_year = preprocessed_path + '/' + str(year)
                        if not os.path.exists(dir_year):
                            os.mkdir(dir_year)
                            print("Directory ", dir_year, " Created ")

                        dir_tile = dir_year + '/' + tile
                        if not os.path.exists(dir_tile):
                            os.mkdir(dir_tile)
                            print("Directory ", dir_tile, " Created ")

                        dir_day = dir_tile + '/' + day
                        if not os.path.exists(dir_day):
                            os.mkdir(dir_day)
                            print("Directory ", dir_day, " Created ")

                        image = rasterio.open("{}/{}/{}/{}".format(qa_masked_path, year, tile, image_path))
                        tif_file = image.read(1)

                        for crop_type in classes:

                            dir_crop = dir_day + '/' + crop_type
                            if not os.path.exists(dir_crop):
                                os.mkdir(dir_crop)
                                print("Directory ", dir_crop, " Created ")

                            mask_path = "{}/{}/{}/land_cover_{}_{}_{}.tif".format(crop_mask_path, year, tile, tile,
                                                                                  year, crop_type)
                            mask = rasterio.open(mask_path)
                            mask_tif = mask.read(1)
                            crop_mask = np.copy(mask_tif)
                            crop = tif_file * crop_mask
                            rastout = "{}/{}/{}/{}/{}/{}_{}_{}_{}_{}_{}.tif".format(preprocessed_path, year, tile, day,
                                                                                    crop_type, sat_id, tile, year, day,
                                                                                    band_no, crop_type)

                            kwargs = image.meta
                            kwargs.update(dtype=rasterio.float32)
                            with rasterio.open(rastout, 'w', **kwargs) as dst:
                                dst.write(crop.astype(rasterio.float32), 1)
                                print("Save Done: %s" % rastout)

    def create_mvc(self, qa_masked_path, sat_id, tile, year, first_day, last_day):
        # Creates Maximum-Value Composite from days list for particular sat_id/tile/year combination
        days = range(int(first_day), int(last_day) + 1)
        days_str = first_day + '-' + last_day

        dir_sat = qa_masked_path + '/' + sat_id
        if not os.path.exists(dir_sat):
            os.mkdir(dir_sat)
            print("Directory ", dir_sat, " Created ")

        dir_year = dir_sat + '/' + year
        if not os.path.exists(dir_year):
            os.mkdir(dir_year)
            print("Directory ", dir_year, " Created ")

        dir_tile = dir_year + '/' + tile
        if not os.path.exists(dir_tile):
            os.mkdir(dir_tile)
            print("Directory ", dir_tile, " Created ")

        file_lists = []
        file_util = FileUtils()
        for root, dirs, files in sorted(os.walk(qa_masked_path)):
            files.sort(key=file_util.natural_keys)
            file_lists.append(files)

        for files in file_lists:
            bands_dict = {str(i): [] for i in range(1, 14)}
            for file in files:
                if ".tif" in file:
                    sat_id_read, tile_read, year_read, day, band_no = file_util.parse_qamasked_file(file)

                    if re.match(r'^\d{3}$', day):
                        if sat_id_read == sat_id and tile_read == tile and year_read == year and int(day) in days:
                            path = f"{qa_masked_path}/{sat_id}/{year}/{tile}/{file}"
                            bands_dict[band_no].append(path)

            for band in bands_dict:
                images = []
                geotrans = []
                projections = []
                for path in bands_dict[band]:
                    file = gdal.Open(path, gdal.GA_ReadOnly)
                    images.append(file.ReadAsArray())
                    geotrans.append(file.GetGeoTransform())
                    projections.append(file.GetProjection())
                    print(path, " read")

                if images:
                    mvc = images[0]
                    for i in range(1, len(images)):
                        mvc = np.maximum(mvc, images[i])

                    same_geotrans = False
                    same_projection = False
                    if geotrans:
                        same_geotrans = all(transform == geotrans[0] for transform in geotrans)
                    if projections:
                        same_projection = all(projection == projections[0] for projection in projections)

                    if same_geotrans and same_projection:
                        geotransform = geotrans[0]
                        projection = projections[0]

                    if not same_geotrans:
                        print("Images do not share same GeoTransform")
                    if not same_projection:
                        print("Images do not share same Projection")

                    format = "GTiff"
                    driver = gdal.GetDriverByName(format)
                    outDataRaster_path = f"{dir_tile}/{sat_id}_{tile}_{year}_{days_str}_{band}.tif"
                    outDataRaster = driver.Create(outDataRaster_path, mvc.shape[1], mvc.shape[0], 1, gdal.GDT_Float32)
                    outDataRaster.SetGeoTransform(geotransform)
                    outDataRaster.SetProjection(projection)
                    outDataRaster.GetRasterBand(1).WriteArray(mvc)
                    print("Save Done: %s" % outDataRaster_path)
                    outDataRaster.FlushCache()
                    del outDataRaster

    def merge(self, qa_masked_path, region_name, sat_id, tiles, year, day):
        # Merges tiles for a specific sat_id/year/day
        dir_sat = qa_masked_path + '/' + sat_id
        if not os.path.exists(dir_sat):
            os.mkdir(dir_sat)
            print("Directory ", dir_sat, " Created ")

        dir_year = dir_sat + '/' + year
        if not os.path.exists(dir_year):
            os.mkdir(dir_year)
            print("Directory ", dir_year, " Created ")

        dir_region = dir_year + '/' + region_name
        if not os.path.exists(dir_region):
            os.mkdir(dir_region)
            print("Directory ", dir_region, " Created ")

        file_lists = []
        file_util = FileUtils()
        for root, dirs, files in sorted(os.walk(qa_masked_path)):
            files.sort(key=file_util.natural_keys)
            file_lists.append(files)

        if sat_id == "S30":
            bands_dict = {str(i): [] for i in range(1, 14)}
            data_type = "Float32"
        elif sat_id == "L30":
            bands_dict = {str(i): [] for i in range(1, 11)}
            data_type = "Float32"
        elif sat_id == "MCD43A4":
            bands_dict = {str(i): [] for i in range(1, 8)}
            data_type = "Int16"

        for files in file_lists:
            for file in files:
                if ".tif" in file:
                    sat_id_read, tile, year_read, days_read, band_no = file_util.parse_qamasked_file(file)

                    if sat_id_read == sat_id and tile in tiles and year_read == year and days_read == day:
                        path_in = f"{qa_masked_path}/{sat_id}/{year}/{tile}/{file}"
                        bands_dict[str(int(band_no))].append(path_in)

        for band in bands_dict:
            if bands_dict[band]:
                paths = " ".join(bands_dict[band])
                print(paths, " read")
                command = f"gdal_merge.py -o {dir_region}/{sat_id}_{region_name}_{year}_{day}_{band}.tif -of gtiff -ot {data_type} -n 0 -tap {paths}"
                output = os.system(command)
                print(output)

    def polygonize(self, rasters_path, vectors_path, tiles, years):
        # Polygonizes rasters in rasters_path directory
        # Vector output is stored in vectors_path directory
        if not os.path.exists(vectors_path):
            os.mkdir(vectors_path)
            print("Directory ", vectors_path, "created ")
        else:
            print("Directory ", vectors_path, "already exists")

        file_lists = []
        file_util = FileUtils()
        for root, dirs, files in sorted(os.walk(rasters_path)):
            files.sort(key=file_util.natural_keys)
            file_lists.append(files)

        for files in file_lists:
            for file in files:
                if ".tif" in file:
                    tile, year, ext = file_util.parse_mask_file(file)

                    if tile in tiles and year in years:
                        dir_year = vectors_path + '/' + year
                        if not os.path.exists(dir_year):
                            os.mkdir(dir_year)
                            print("Directory ", dir_year, " Created ")

                        path_in = f"{rasters_path}/{year}/{file}"
                        path_out = f"{vectors_path}/{year}/land_cover_{tile}_{year}.geojson"

                        mask = None
                        src = rasterio.open(path_in)
                        image = src.read(1)
                        results = ({'properties': {'ID_PARCEL': i, 'CODE_GROUP': int(v)}, 'geometry': s} for i, (s, v)
                                   in enumerate(shapes(image, mask=mask, transform=src.transform)))

                        with fiona.open(path_out, 'w', driver="GeoJSON", crs=src.crs,
                                        schema={'properties': [('ID_PARCEL', 'int'), ('CODE_GROUP', 'int')],
                                                'geometry': 'Polygon'}) as dst:
                            dst.writerecords(results)
                        print("Save Done: %s" % path_out)

    def stack_rasters(self, path_in, path_out, tile, year, day):
        # Stacks rasters in path_in directory
        # Stores new rasters in path_out directory
        if not os.path.exists(path_out):
            os.mkdir(path_out)
            print("Directory ", path_out, "created ")
        else:
            print("Directory ", path_out, "already exists")

        file_lists = []
        file_util = FileUtils()
        for root, dirs, files in sorted(os.walk(path_in)):
            files.sort(key=file_util.natural_keys)
            file_lists.append(files)

        bands = []
        for files in file_lists:
            for file in files:
                if ".tif" in file:
                    sat_id, tile_read, year_read, day_read, band_no = file_util.parse_qamasked_file(file)

                    if tile_read == tile and year_read == year and day_read == day:
                        dir_year = path_out + '/' + year
                        if not os.path.exists(dir_year):
                            os.mkdir(dir_year)
                            print("Directory ", dir_year, " Created ")

                        dir_tile = dir_year + '/' + tile
                        if not os.path.exists(dir_tile):
                            os.mkdir(dir_tile)
                            print("Directory ", dir_tile, " Created ")

                        path = f"{path_in}/{year}/{tile}/{file}"
                        path = f"{path_in}/{file}"
                        band = gdal.Open(path, gdal.GA_ReadOnly)
                        bands.append(band.ReadAsArray())
                        print(path, " read")

        if bands:
            bands = np.stack(bands, axis=0)
            geotransform = band.GetGeoTransform()
            projection = band.GetProjection()

            date = datetime.datetime.strptime(f"{year}{day}", '%Y%j').date()
            date = date.strftime('%Y%m%d')

            format = "GTiff"
            driver = gdal.GetDriverByName(format)
            outDataRaster_path = f"{path_out}/{year}/{tile}/{date}.tif"
            outDataRaster = driver.Create(outDataRaster_path, bands.shape[2], bands.shape[1], bands.shape[0],
                                          gdal.GDT_Float32)
            outDataRaster.SetGeoTransform(geotransform)
            outDataRaster.SetProjection(projection)
            for i, band in enumerate(bands, 1):
                outDataRaster.GetRasterBand(i).WriteArray(band)
            print("Save Done: %s" % outDataRaster_path)
            outDataRaster.FlushCache()
            del outDataRaster

    def reformat_polygon_file(self, vector_path_in, vector_path_out, raster_path):
        vector_in = gpd.read_file(vector_path_in)

        raster = gdal.Open(raster_path)
        ulx, xres, xskew, uly, yskew, yres = raster.GetGeoTransform()
        lrx = ulx + (raster.RasterXSize * xres)
        lry = uly + (raster.RasterYSize * yres)
        coords = [(ulx, lry), (ulx, uly), (lrx, uly), (lrx, lry), (ulx, lry)]
        raster_poly = Polygon(coords)
        raster = rasterio.open(raster_path)

        features = []
        parcel_id = 0
        for polygon in vector_in['geometry']:
            if polygon.centroid.within(raster_poly):
                crop, transform = rasterio.mask.mask(raster, [polygon], crop=True)
                label = stats.mode(crop.flatten())[0][0]
                features.append(
                    Feature(geometry=polygon, properties={"ID_PARCEL": parcel_id, "CODE_GROUP": int(label)}))
                parcel_id = parcel_id + 1

        crs = {"type": "name", "properties": {"name": "EPSG:32614"}}
        feature_collection = FeatureCollection(features, crs=crs)
        with open(vector_path_out, 'w') as file:
            dump(feature_collection, file)
        print("Save Done: ", vector_path_out)

    def create_crop_mask_aafc(self, path_to_aafc_crop_map, crop_id_list):

        """
        Create mask layers out of AAFC dataset for the provided crop id and store the masked layer as a new tif file
        :param: path_to_aafc_crop_map: pathway to the AAFC tif file
        :param: crop_id_list:list of crop ids to be masked based on the information provided for crop ids
        e.g. [145] or [145, 146]
        """

        # creates crop mask filter and store the tiff file containing only the requested crop id

        # reads AAFC raster file
        src = rasterio.open(path_to_aafc_crop_map)
        src_array = src.read(1)

        src_array_filter_list = []
        for crop_id in crop_id_list:
            src_array_filter = src_array == crop_id
            src_array_filter_list.append(src_array_filter)

        src_array[:] = 0
        for i in range(0, len(crop_id_list)):
            src_array[src_array_filter_list[i]] = crop_id_list[i]

        # write the masked AAFC raster as a new tiff file
        with rasterio.Env():
            # Write an array as a raster band to a new 8-bit file. For
            # the new file's profile, we start with the profile of the source
            profile = src.profile

            # And then change the band count to 1, set the
            # dtype to uint8, and specify LZW compression.
            profile.update(
                dtype=rasterio.uint8,
                count=1,
                compress='lzw')

            path_to_masked_aafc = ''
            for i in path_to_aafc_crop_map.split('/')[:-1]: path_to_masked_aafc += '%s/' % i

            output_name = 'masked_cropID_'
            for crop_id in crop_id_list:
                output_name += '_%s' % crop_id

            path_to_masked_aafc += output_name + '.tif'

            with rasterio.open(path_to_masked_aafc, 'w', **profile) as dst:
                dst.write(src_array.astype(rasterio.uint8), 1)
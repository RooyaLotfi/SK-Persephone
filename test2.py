import os
from rasterio.warp import reproject, Resampling
from rasterio.warp import calculate_default_transform
import matplotlib.pyplot
import rasterio.mask
import rasterio as rio
from rasterio.plot import show
import earthpy.plot as ep
from osgeo import gdal, gdalnumeric, ogr
import rasterio
import geopandas as gpd
import fiona
from shapely.geometry import Polygon, Point
#from geojson import Feature, FeatureCollection, dump
from scipy import stats
import numpy as np
import numpy.ma as ma  # For reproject_tif ; For mask_data


def reproject_et(inpath, outpath, new_crs):
    '''
    I have to change this in a way that it would not record it in a path and just gets a tiff and returns a tiff
    with a new crs

    This function reprojects the raster to a new crs
    :param inpath: pathway to the input tif file
    :param outpath: pathway to the output tif file
    :param new_crs: The new crs
    '''
    dst_crs = new_crs
    with rio.open(inpath) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })
        # fix these parts because the path might not exist
        # if not(os.path.exists(outpath)):
        #     os.mkdir(outpath)
        with rio.open(outpath, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(source=rio.band(src, i),
                          destination=rio.band(dst, i),
                          src_transform=src.transform,
                          src_crs=src.crs,
                          dst_transform=transform,
                          dst_crs=dst_crs,
                          resampling=Resampling.nearest)


def change_tiff_crs(path_in, AAFC_raster, dst_crs):
    """

    :param path_in: path to raster file: AAFC_raster
    :param AAFC_raster: the raster file itself
    :param dst_crs: final crs that we want to transform
    :return: nothing. the AAFC_raster file changes completely (it writes on it with a new crs)
    """
    with rasterio.open(path_in, 'w', **kwargs) as dst:
        for i in range(1, AAFC_raster.count + 1):
            reproject(
                source=rasterio.band(AAFC_raster, i),
                destination=rasterio.band(dst, i),
                src_transform=AAFC_raster.transform,
                src_crs=AAFC_raster.crs,
                dst_transform=transform,
                dst_crs=dst_crs,
                resampling=Resampling.nearest)


def create_crop_mask_aafc(path_to_aafc_crop_map, crop_id_list):
    """
    Create mask layers out of AAFC dataset for the provided crop id and store the masked layer as a new tif file
    :param: path_to_aafc_crop_map: pathway to the AAFC tif file
    :param: crop_id_list:list of crop ids to be masked based on the information provided for crop ids
    e.g. [145] or [145, 146]
    """

    # creates crop mask filter and store the tiff file containing only the requested crop id

    # reads AAFC raster file
    src = rio.open(path_to_aafc_crop_map)
    src_array = src.read(1)

    src_array_filter_list = []
    for crop_id in crop_id_list:
        src_array_filter = src_array == crop_id
        src_array_filter_list.append(src_array_filter)

    src_array[:] = 0
    for i in range(0, len(crop_id_list)):
        src_array[src_array_filter_list[i]] = crop_id_list[i]

    # write the masked AAFC raster as a new tiff file
    with rio.Env():
        # Write an array as a raster band to a new 8-bit file. For
        # the new file's profile, we start with the profile of the source
        profile = src.profile

        # And then change the band count to 1, set the
        # dtype to uint8, and specify LZW compression.
        profile.update(
            dtype=rio.uint8,
            count=1,
            compress='lzw')

        path_to_masked_aafc = ''
        for i in path_to_aafc_crop_map.split('/')[:-1]: path_to_masked_aafc += '%s/' % i

        output_name = 'masked_cropID_'
        for crop_id in crop_id_list:
            output_name += '_%s' % crop_id

        path_to_masked_aafc += output_name + '.tif'

        with rio.open(path_to_masked_aafc, 'w', **profile) as dst:
            dst.write(src_array.astype(rio.uint8), 1)
        return path_to_masked_aafc


def reformat_polygon_file(vector_path_in, vector_path_out, raster_path):
    vector_in = gpd.read_file(vector_path_in)
    raster = rasterio.open(raster_path)

    features = []
    parcel_id = 0
    for polygon in vector_in['geometry']:
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


if __name__ == '__main__':
    path_in = os.path.join("AAFC_data", "aci_2019_sk_v1", "aci_2019_sk.tif")
    CRS = "EPSG:32613"
    _output_crs = CRS.replace(":", "_")

    AAFC_raster_path = os.path.join("AAFC_data", "aci_2019_sk_v1", _output_crs+"_aci_2019_sk.tif")

    reproject_et(path_in, AAFC_raster_path, CRS)
    AAFC_raster = rio.open(AAFC_raster_path)
    print(AAFC_raster.crs)
    aafc_crop_mask_path = create_crop_mask_aafc(AAFC_raster_path, [146, 153])

    AAFC_crop_mask = rio.open(aafc_crop_mask_path)


    # #
    # # output_path = os.path.join("polygon", "boundaries.shp")
    # # polygons = gpd.read_file(output_path)
    #
    # # print("crs three is ", polygons.crs)
    # #
    # # reproject_et(path, "new_AAFCC_CRS_EPSG:32614.tif", "EPSG:32614")
    # #
    # # src2 = rio.open("new_AAFCC_CRS_EPSG:32614.tif")
    # # print(src2.crs)
    #
    # # # matplotlib.pyplot.imshow(src2.read(1), cmap='pink')
    # # # matplotlib.pyplot.show()
    # #
    # # shapes = polygons.geometry
    # # # num = 0
    # # # for shape in shapes:
    # # #     print(num)
    # # #     num = num + 1
    # # #     print(shape)
    # #
    # # with fiona.open(output_path, "r") as shapefile:
    # #     shapes = [feature["geometry"] for feature in shapefile]
    # #
    # # with rasterio.open("aci2.tif") as src:
    # #     out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
    # #     out_meta = src.meta
    # #
    # # out_meta.update({"driver": "GTiff",
    # #                  "height": out_image.shape[1],
    # #                  "width": out_image.shape[2],
    # #                  "transform": out_transform})
    # #
    # # with rasterio.open("cropped_masked.tif", "w", **out_meta) as dest:
    # #     dest.write(out_image)
    # #
    # # src3 = rio.open("cropped_masked.tif")
    # #
    # # polygon_0 = polygons.loc[[0], 'geometry']
    # # polygon_1 = polygons.loc[[1], 'geometry']
    # # polygon_20 = polygons.loc[[20], 'geometry']
    # # polygon_200 = polygons.loc[[200], 'geometry']
    # # polygon_150 = polygons.loc[[150], 'geometry']
    # #
    # # extent = rio.plot.plotting_extent(src3)
    # #
    # # #
    # # # matplotlib.pyplot.imshow(src3.read(1))
    # # #
    # # # # plot the data
    # # # fig, ax = matplotlib.pyplot.subplots(figsize=(6, 6))
    # # # polygon_0.plot(ax=ax)
    # # # ax.set_title("Shapefile Imported into Python \nCrop Extent for Soaproot Saddle Field Site",
    # # #              fontsize=16)
    # # # #ax.set_axis_off();
    # # # matplotlib.pyplot.show()
    # # #
    # # # print(src3.crs)
    # # # print(polygons.crs)
    # # #
    # # fig, ax = matplotlib.pyplot.subplots(figsize=(10, 10))
    # # #plt.plot(src3.read(1))
    # # ep.plot_bands(src3.read(1),
    # #               cmap='terrain',
    # #               extent=extent,
    # #               ax=ax,
    # #               cbar=False)
    # # polygon_0.plot(ax=ax, alpha=.6, color='g')
    # # polygon_1.plot(ax=ax, alpha=.6, color='b')
    # # polygon_20.plot(ax=ax, alpha=.6, color='r')
    # # polygon_200.plot(ax=ax, alpha=.6, color='g')
    # # polygon_150.plot(ax=ax, alpha=.6, color='g')
    # # polygons.plot(ax=ax, alpha=.6, color ='g')
    # #
    # # matplotlib.pyplot.show()
    # # #
    # # # """**********************************"""
    # #
    # # fig, ax = matplotlib.pyplot.subplots(figsize=(10, 10))
    # # # plt.plot(src3.read(1))
    # # ep.plot_bands(src2.read(1),
    # #               cmap='terrain',
    # #               extent=extent,
    # #               ax=ax,
    # #               cbar=False)
    # # polygon_0.plot(ax=ax, alpha=.6, color='g')
    # # polygon_1.plot(ax=ax, alpha=.6, color='b')
    # # polygon_20.plot(ax=ax, alpha=.6, color='r')
    # # polygon_200.plot(ax=ax, alpha=.6, color='g')
    # # polygon_150.plot(ax=ax, alpha=.6, color='g')
    # # polygons.plot(ax=ax, alpha=.6, color='g')
    # #
    # # matplotlib.pyplot.show()
    # # # # print(polygons.count)
    # # # # for i in range(polygons.shape[0]):
    # # # #     print(i)
    # #
    # #
    # # # """for plotting a small mask :D """
    # # with rasterio.open("aci2.tif") as src:
    # #     out_image, out_transform = rasterio.mask.mask(src, polygon_150, crop=True)
    # #     out_meta = src.meta
    # #
    # # out_meta.update({"driver": "GTiff",
    # #                  "height": out_image.shape[1],
    # #                  "width": out_image.shape[2],
    # #                  "transform": out_transform})
    # # with rasterio.open("small_masked.tif", "w", **out_meta) as dest:
    # #     dest.write(out_image)
    # # #
    # #
    # # fig, ax = matplotlib.pyplot.subplots(2, figsize=(10, 10))
    # # src2 = rio.open("small_masked.tif")
    # # #plt.plot(src3.read(1))
    # # ep.plot_bands(src2.read(1),
    # #               cmap='terrain',
    # #               extent=extent,
    # #               ax=ax[0],
    # #               cbar=False)
    # # polygon_150.plot(ax=ax[1], alpha=.6, color ='r')
    # # matplotlib.pyplot.show()
    # #
    # #
    # #

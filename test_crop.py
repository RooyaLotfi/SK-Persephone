import geopandas as gpd
import os
import pathlib
import gdal
import rasterio as rio
import matplotlib.pyplot as plt
import earthpy.plot as ep
from rasterio.plot import show
import rioxarray as rxr
from shapely.geometry import mapping
import earthpy.spatial as es
import numpy as np


def crop_raster(input_raster, input_shapefile, output_raster):
    # Create convex polygon
    df = gpd.read_file(input_shapefile)
    polygon_geom = df.unary_union.convex_hull #TODO (Juan): This results in, well, a convex polygon, and not as precise as can be
    crs = 'epsg:32613'
    polygon = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[polygon_geom])
    # print(polygon.geometry)
    # Save polygon file
    polygon_dir = os.path.join('output','cropping_shapefiles', os.path.basename(input_shapefile))
    pathlib.Path(polygon_dir).parents[0].mkdir(parents = True, exist_ok = True)
    polygon.to_file(filename=polygon_dir, driver="ESRI Shapefile")
    warp_options = gdal.WarpOptions(format='GTiff',
        cropToCutline=True,
        cutlineDSName=polygon_dir,
        dstAlpha=True,
        dstNodata=0,
        srcSRS='epsg:32613',
        dstSRS='epsg:32613'
    )
    try:
        out = gdal.Warp(output_raster, input_raster, options=warp_options)
        out.FlushCache()
        out = None
        del out
    except:
        print('crop_raster: Could not convert:', input_raster)
    os.remove(polygon_dir)
    return output_raster


# Which RM? pick it
# Which QA path?
# Find out which tiles consist of it
# Merge those tiles or not
# Crop them and multiply it to
#

rm_boundary_path = os.path.join("polygon_EPSG_32613", "cropped_by_CA_boundaries", "polygon.shp")
rm_boundary = gpd.read_file(rm_boundary_path)

sk_mgrs_rm_path = os.path.join("sk-mgrs-rm", "sk-mgrs-rm.shp")
mgrs_rm = gpd.read_file(sk_mgrs_rm_path)


RM_ID = 232

# Select shape file for Raster ID RM_ID
RM_232 = mgrs_rm.loc[mgrs_rm['id'] == RM_ID]
print(RM_232.crs)
raster_path = os.path.join("data", "2_qa_masked", "S30", "2019", "T12UWA_T12UWB_T12UWC", "S30_T12UWA_T12UWB_T12UWC_2019_154_1.tif")

lidar_dem = rxr.open_rasterio(raster_path,
                              masked=True).squeeze()
# Check the CRS
print("1_original one ", lidar_dem.rio.crs)
lidar_dem_wgs84 = lidar_dem.rio.reproject(RM_232.crs)
print("transformed one ",lidar_dem_wgs84.rio.crs)


# Plot your newly converted data
f, ax = plt.subplots(figsize=(10, 4))

lidar_dem_wgs84.plot.imshow(ax=ax,
                            cmap='Greys')
RM_232.plot(ax=ax, alpha=.6, color='g')
ax.set(title="Plot Showing Roads Overlayed on Elevation Data")
ax.set_axis_off()
plt.show()


lidar_clipped = lidar_dem_wgs84.rio.clip(RM_232.geometry.apply(mapping))

f, ax = plt.subplots(figsize=(10, 4))
lidar_clipped.plot(ax=ax)
RM_232.plot(ax=ax, alpha=.6, color='g')

ax.set(title="Raster Layer Cropped to Geodataframe Extent")
ax.set_axis_off()
plt.show()


# with rio.open(raster_path) as src:
#     lidar_chm_im = src.read(masked=True)[0]
#     extent = rio.plot.plotting_extent(src)
#     soap_profile = src.profile
#
# RM_232 = RM_232.to_crs(src.crs)
# print(src.crs)
# print(RM_232.crs)
# fig, ax = plt.subplots(figsize=(10, 10))
# ep.plot_bands(lidar_chm_im,
#               cmap='terrain',
#               extent=extent,
#               ax=ax,
#               cbar=False)
# RM_232.plot(ax=ax, alpha=.6, color='g')
# plt.show()
crop_extent_soap = RM_232.to_crs(lidar_dem.rio.crs)

with rio.open(raster_path) as src:
    lidar_chm_crop, soap_lidar_meta = es.crop_image(src, crop_extent_soap)

    # Update the metadata to have the new shape (x and y and affine information)
soap_lidar_meta.update({"driver": "GTiff",
                        "height": lidar_chm_crop.shape[0],
                        "width": lidar_chm_crop.shape[1],
                        "transform": soap_lidar_meta["transform"]})

# generate an extent for the newly cropped object for plotting
cr_ext = rio.transform.array_bounds(soap_lidar_meta['height'],
                                    soap_lidar_meta['width'],
                                    soap_lidar_meta['transform'])

bound_order = [0, 2, 1, 3]
cr_extent = [cr_ext[b] for b in bound_order]
cr_extent, crop_extent_soap.total_bounds

# mask the nodata and plot the newly cropped raster layer
lidar_chm_crop_ma = np.ma.masked_equal(lidar_chm_crop[0], -9999.0)
ep.plot_bands(lidar_chm_crop_ma, cmap='terrain', cbar=False)
# Save to disk so you can use the file later.

path_out = "data/soap_lidar_chm_crop.tif"
with rio.open(path_out, 'w', **soap_lidar_meta) as ff:
    ff.write(lidar_chm_crop_ma, 1)
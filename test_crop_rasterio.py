import rioxarray
import os
import matplotlib.pyplot as plt
import geopandas as gpd
import earthpy.plot as ep
from shapely.geometry import mapping


raster_path = os.path.join("data", "2_qa_masked", "S30", "2019", "T12UWA_T12UWB_T12UWC",
                           "S30_T12UWA_T12UWB_T12UWC_2019_154_1.tif")
CRS = "epsg:32613"
xarray = rioxarray.open_rasterio(raster_path)
xarray = xarray.rio.reproject(CRS)
print(xarray.rio.crs)
print(xarray.data[0])

f, ax = plt.subplots(figsize=(10, 10))
xarray[0].plot.imshow()
ax.set(title="S30_T12UWA_T12UWB_T12UWC_2019_154_1")

ax.set_axis_off()
plt.show()
sk_mgrs_rm_path = os.path.join("sk-mgrs-rm", "sk-mgrs-rm.shp")
mgrs_rm = gpd.read_file(sk_mgrs_rm_path)
RM_ID = 232
# Select shape file for Raster ID RM_ID
RM_232 = mgrs_rm.loc[mgrs_rm['id'] == RM_ID]

f, ax = plt.subplots(figsize=(10, 5))
xarray[0].plot.imshow(ax=ax)

RM_232.plot(ax=ax,
                 alpha=.5)
ax.set(title="Raster Layer with Shapefile Overlayed")

ax.set_axis_off()
plt.show()

lidar_clipped = xarray.rio.clip(RM_232.geometry.apply(mapping))

f, ax = plt.subplots(figsize=(10, 10))
lidar_clipped.plot(ax=ax)
ax.set(title="Raster Layer Cropped to Geodataframe Extent")
ax.set_axis_off()
plt.show()

# Write the data to a new geotiff file
lidar_clipped.rio.to_raster("rasteriocropped.tif")

# print(RM_232)

# with rio.open(raster_path) as src:
#     # avalin band ro mikhone
#     band1 = src.read(masked=True)[0]
#     # tool va arze joghrafiayi
#     extent = rio.plot.plotting_extent(src)
#     soap_profile = src.profile
#     # print(extent)
#     # print(band1)
#     # print(soap_profile)
#
# ep.plot_bands(band1,
#               extent=extent,
#               title="S30 T12UWA + T12UWB + T12UWC 2019 154 1",
#               )
#

# #
# print('crop extent crs: ', RM_232.crs)
# # print('lidar crs: ', soap_profile['crs'])
#
# # plot the data
# fig, ax = plt.subplots(figsize = (6, 6))
# RM_232.plot(ax=ax)
# ax.set_title("Shapefile Imported into Python \nCrop Extent for RM",
#              fontsize = 16)
# ax.set_axis_off();
# plt.show()
#
#
# fig, ax = plt.subplots(figsize=(10, 10))
# ep.plot_bands(band1,
#               cmap='terrain',
#               extent=extent,
#               ax=ax,
#               cbar=False)
# RM_232.plot(ax=ax, alpha=.6, color='g');
# plt.show()
import os
import rasterio as rio
import geopandas as gpd
import matplotlib.pyplot as plt
from rasterio.warp import calculate_default_transform, reproject, Resampling
import fiona
import rasterio
import rasterio.mask


# Default values for satellite, tile, extension of files to download, and year
SATELLITE = 'L30'
# First three characters of the tiles
TILE = '12UXA'
YEAR = '2019'
CRS = "epsg:32613"
CROP_ID = 146

dst_crs = CRS
hls_file_src = "HLS.L30.T12UXA.2019001.v1.4.01.tif"
tile_raster = os.path.join("data", "hls_tmp_tif", TILE, YEAR, SATELLITE, "HLS.L30.T12UXA.2019001.v1.4.09.tif")
hls_file_dst = os.path.join("data", "hls_tmp_tif", TILE, YEAR, SATELLITE,
                            "HLS.L30.T12UXA.2019001.v1.4.09.tif." + dst_crs + ".tif")

with rio.open(tile_raster) as src:
    transform, width, height = calculate_default_transform(
        src.crs, dst_crs, src.width, src.height, *src.bounds)
    kwargs = src.meta.copy()
    kwargs.update({
        'crs': dst_crs,
        'transform': transform,
        'width': width,
        'height': height
    })

    with rio.open(hls_file_dst, 'w', **kwargs) as dst:
        for i in range(1, src.count + 1):
            reproject(
                source=rio.band(src, i),
                destination=rio.band(dst, i),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=transform,
                dst_crs=dst_crs,
                resampling=Resampling.nearest)

crop_type_path = os.path.join("AAFC_data", "aci_2019_sk_v1", "masked_cropID__146_153.tif")
crop_mask = rio.open(crop_type_path)

rm_tile_path = os.path.join("sk-mgrs-rm", "sk-mgrs-rm.shp")
rm_vor_tile = gpd.read_file(rm_tile_path)

tile_raster = rio.open(hls_file_dst)
out_meta = tile_raster.meta
# tile_raster.to_crs(CRS)

# print(crop_mask.read(1))
# print(rm_vor_tile.columns)
# print(rm_vor_tile.loc[rm_vor_tile['MGRS'] == TILE])

filter_rm_tile = rm_vor_tile.loc[rm_vor_tile['MGRS'] == TILE]
# filter_rm_tile.plot()
# plt.show()
# print("checking crs : \nrm_tile: ", rm_vor_tile.crs)
# print("crop_type: ", crop_mask.crs)
# print("tile_raster: ", tile_raster.crs)

polygons = filter_rm_tile["geometry"]
print(filter_rm_tile)
#print(filter_rm_tile.columns)
for row in filter_rm_tile.iterrows():
    #polygon = row['geometry']
    #print(row['geometry'])
    print(row)
    shapes = [feature["geometry"] for feature in row]
    out_image, out_transform = rasterio.mask.mask(tile_raster, shapes, crop=True)
    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})

    cropped_file_dst = os.path.join("data", "hls_tmp_tif", TILE, YEAR, SATELLITE, "cropped"
                                "HLS.L30.T12UXA.2019001.v1.4.01." + index + ".tif")
    with rasterio.open(cropped_file_dst, "w", **out_meta) as dest:
        dest.write(out_image)

    #x, y = polygon.exterior.xy
#     # plt.plot(x, y)
#     # plt.show()
#     print(index)
#     out_image, out_transform = rio.mask.mask(tile_raster, polygon, crop=True)
#     out_meta = tile_raster.meta
#     out_meta.update({"driver": "GTiff",
#                      "height": out_image.shape[1],
#                      "width": out_image.shape[2],
#                      "transform": out_transform})
#     name = polygon.id
#     with rasterio.open("", "w", **out_meta) as dest:
#         dest.write(out_image)

#
# with rasterio.open("tests/data/RGB.byte.tif") as src:
#     out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
#     out_meta = src.meta
import os
import geopandas as gpd
import rioxarray as rxr
from Utils.crop_soilgrid_smap import crop
import matplotlib.pyplot as plt

sk_boundary_path = os.path.join("Canada-Boundaries", "Saskatchewan", "sk_boundary.shp")
sk_boundary = gpd.read_file(sk_boundary_path)
CRS = "EPSG:32613"
sk_boundary = sk_boundary.to_crs(CRS)
smap_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone/data/smap_originalResolution/2019")
image_file = os.path.join("sm_rootzone_analysis_20150401.tif")
buffer = 40000
sk_boundary_buffer = sk_boundary.buffer(buffer)


cropped_smap_path = os.path.join(smap_path, "cropped_sk_boundary_buffer_"+str(buffer))
if not os.path.exists(cropped_smap_path):
    os.makedirs(cropped_smap_path)

for root, dirs, image_files in sorted(os.walk(smap_path)):
    for image_file in image_files:
        if '.tif' in image_file:
            print(image_file)
            smap = rxr.open_rasterio(os.path.join(smap_path, image_file))
            cropped_smap = crop(smap, sk_boundary_buffer)
            cropped_smap.rio.to_raster(os.path.join(cropped_smap_path, image_file))

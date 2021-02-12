import os
import geopandas as gpd
from Utils.resample import resample


#
target_resolution = 30
resampling_method = 'bilinear'


# Pick either one of these paths:
smap_original_resolution_path = os.path.join(
    "/Users/roya.lotfi/AIRM_Project/SK-Persephone/data/smap_originalResolution_selected/2019")
SMAP_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone/data/smap_originalResolution_selected/test_resample/")
rasters_path = smap_original_resolution_path
resampled_smap_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone/data/smap_originalResolution_selected/test_resample/")
if not os.path.exists(resampled_smap_path):
    os.makedirs(resampled_smap_path)
for root, dirs, image_files in sorted(os.walk(smap_original_resolution_path)):
    for image_file in image_files:
        if '.tif' in image_file:
            image_path = os.path.join(root, image_file)
            resampled_raster_path = os.path.join(resampled_smap_path, image_file)
            print(resampled_raster_path)
            resample(image_path, resampled_raster_path, target_resolution=target_resolution, resampling_method=resampling_method)

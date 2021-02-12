import os
import geopandas as gpd
from Utils.resample import resample


target_resolution = 30
resampling_method = 'bilinear'


def main():
    soilgrid_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone/data/Soli_Grid_originalResolution")
    resampled_rasters_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone/data/Soil_Grid/2019")
    if not os.path.exists(resampled_rasters_path):
        os.makedirs(resampled_rasters_path)
    for root, dirs, image_files in sorted(os.walk(soilgrid_path)):
        for image_file in image_files:
            if '.tif' in image_file:
                image_path = os.path.join(root, image_file)
                resampled_raster_path = os.path.join(resampled_rasters_path, image_file)
                resample(image_path, resampled_raster_path, target_resolution, resampling_method)


if __name__ == '__main__':
    main()

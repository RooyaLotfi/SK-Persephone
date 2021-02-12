import os
import geopandas as gpd
from Utils.resample import resample

#
target_resolution = 30
resampling_method = 'bilinear'


def main():
    YEAR = ['2019']

    for year in YEAR:
        # Pick either one of these paths:
        modis_original_resolution_path = os.path.join(
            "/Users/roya.lotfi/AIRM_Project/SK-Persephone/data/Modis_tif/4_cropped_sk", year)
        resampled_modis_path = os.path.join(
            "/Users/roya.lotfi/AIRM_Project/SK-Persephone/data/Modis_tif/5_resampled", year)
        if not os.path.exists(resampled_modis_path):
            os.makedirs(resampled_modis_path)
        for root, dirs, image_files in sorted(os.walk(modis_original_resolution_path)):
            for image_file in image_files:
                if '.tif' in image_file:
                    image_path = os.path.join(root, image_file)
                    resampled_raster_path = os.path.join(resampled_modis_path, image_file)
                    print(resampled_raster_path)
                    resample(image_path, resampled_raster_path, target_resolution=target_resolution,
                             resampling_method=resampling_method)


if __name__ == '__main__':
    main()

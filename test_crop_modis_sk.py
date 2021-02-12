import os
import geopandas as gpd
import rioxarray as rxr
from Utils.crop_soilgrid_smap import crop
import matplotlib.pyplot as plt


def main():
    sk_boundary_path = os.path.join("Canada-Boundaries", "Saskatchewan", "sk_boundary.shp")
    sk_boundary = gpd.read_file(sk_boundary_path)
    CRS = "EPSG:32613"
    sk_boundary = sk_boundary.to_crs(CRS)
    YEAR = ['2019']
    for year in YEAR:
        modis_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone/data/Modis_tif/3_merged", year)
        buffer = 40000
        sk_boundary_buffer = sk_boundary.buffer(buffer)

        cropped_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone/data/Modis_tif/", "4_cropped_sk", year)
        if not os.path.exists(cropped_path):
            os.makedirs(cropped_path)

        for root, dirs, image_files in sorted(os.walk(modis_path)):
            for image_file in image_files:
                if '.tif' in image_file:
                    smap = rxr.open_rasterio(os.path.join(modis_path, image_file))
                    cropped_modis = crop(smap, sk_boundary_buffer)
                    print("cropped path is ", cropped_path)
                    cropped_modis.rio.to_raster(os.path.join(cropped_path, image_file))


if __name__ == '__main__':
    main()

import os
import rioxarray
import geopandas as gpd
import numpy as np
from Utils.crop_soilgrid_smap import crop, plot_raster_polygon
from Utils import crop_soilgrid_smap

CROP_TYPE = 146
#CROP_TYPE = 153


def main():
    curr_path = os.getcwd()
    soilgrid_path = os.path.join(curr_path, "data", "Soil_Grid")
    print(curr_path)

    # path to RM boundaries
    rm_polygons_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone",
                                    "polygon_EPSG_32613/cropped_by_CA_boundaries/", "polygon.shp")
    # Path to MGRS RM boundaries which has an attribute that tells the RM is in which MGRS tiles
    mgrs_rm_polygons_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone/sk-mgrs-rm", "sk-mgrs-rm.shp")
    # We want to crop both the Soilgrid data and the crop type raster file because we need to mask raster files and
    # keep pixels that the crop grows in it
    crop_type_path = os.path.join("AAFC_data", "aci_2019_sk_v1", "masked_cropID__146_153.tif")

    # which IDs we want to work with
    crop_soilgrid_smap.RM_IDs = [322, 232, 18]
    # Which year we want to work with
    crop_soilgrid_smap.YEAR = {'2019'}
    crop_soilgrid_smap.crop_RM(soilgrid_path, rm_polygons_path, mgrs_rm_polygons_path, crop_type_path)


if __name__ == '__main__':
    main()

import os
import geopandas as gpd
from Utils.resample import resample

SATTELITE = {"L30", "S30"}
YEAR = {"2019"}
TILE = {"T12UWB", "T12UWC"}
RM_IDs = [322, 232]
resolution = 30


# TODO: replace the TILE with the tiles that concatenates tiles based on MGRS IDs
def main():
    rm_polygons_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone",
                                    "polygon_EPSG_32613/cropped_by_CA_boundaries/", "polygon.shp")
    mgrs_rm_polygons_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone/sk-mgrs-rm", "sk-mgrs-rm.shp")
    rm_polygons = gpd.read_file(rm_polygons_path)
    mgrs_rm_polygon = gpd.read_file(mgrs_rm_polygons_path)

    hls_path = os.path.join("data", "2_qa_masked")
    for rm_id in RM_IDs:
        mgrs_tiles = mgrs_rm_polygon.loc[mgrs_rm_polygon['id'] == rm_id]['MGRS']
        if mgrs_tiles.empty:
            continue
        tiles_list = ["T" + tile for tile in mgrs_tiles]
        tiles_path = '_'.join(tiles_list)

        for sattelite in SATTELITE:
            for year in YEAR:
                hls_path_original_resolution = os.path.join(hls_path, sattelite, year, tiles_path,
                                                            "original_res")
                hls_path_resampled = os.path.join(hls_path, sattelite, year, tiles_path, "res_" + resolution)
                print("1_original path", hls_path_original_resolution)
                if not os.path.exists(hls_path_original_resolution):
                    print("the path did not exist")
                for root, dir, files in os.walk(hls_path_original_resolution):
                    for file in files:
                        if '.tif' in file:
                            path_in = os.path.join(hls_path_original_resolution, file)
                            path_out = os.path.join(hls_path_resampled, file)
                            resample(path_in, path_out, out_crs='EPSG:32613',
                                     target_resolution=resolution, resampling_method='bilinear')

                # print("1_original path", hls_path_original_resolution)
                # hls_path_resampled = os.path.join(hls_path, sattelite, year, tiles_path)
                # print("5_resampled path", hls_path_resampled)


if __name__ == '__main__':
    main()

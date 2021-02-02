import os
from utilities.FileUtils import FileUtils


def merge(path, dir_path, region_name, sat_id, tiles, year, day):
    # Merges tiles for a specific sat_id/year/day

    dir_path = os.path.join(dir_path, year)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    file_lists = []
    file_util = FileUtils()
    for root, dirs, files in sorted(os.walk(path)):
        files.sort(key=file_util.natural_keys)
        file_lists.append(files)
    bands_dict = []
    if sat_id == "S30":
        bands_dict = {str(i): [] for i in range(1, 14)}
        data_type = "Float32"
    elif sat_id == "L30":
        bands_dict = {str(i): [] for i in range(1, 11)}
        data_type = "Float32"
    elif sat_id == "MCD43A4":
        bands_dict = {str(i): [] for i in range(1, 8)}
        data_type = "Int16"
    elif sat_id == "MCD15A3H":
        bands_dict = {"Fpar_500m": [], "FparExtra_QC": [], "FparStdDev_500m": [], "Lai_500m": [], "LaiStdDev_500m": []}
        data_type = "Int16"

    for files in file_lists:
        for file in files:
            if ".tif" in file:
                sat_id_read, tile, year_read, days_read, band_no = file_util.parse_tile_file(file)
                if sat_id_read == sat_id and tile in tiles and year_read == year and days_read == day:
                    path_in = f"{path}/{file}"
                    bands_dict[str(band_no)].append(path_in)
    for band in bands_dict:
        if bands_dict[band]:
            paths = " ".join(bands_dict[band])
            print(paths, " read")
            command = f"gdal_merge.py -o {dir_path}/{sat_id}_{region_name}_{year}_{day}_{band}.tif -of gtiff -ot {data_type} -n 0 -tap {paths}"
            output = os.system(command)
            print(output)

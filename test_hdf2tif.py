import os
from osgeo import gdal, gdalnumeric, ogr

import re


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    return [atoi(c) for c in re.split(r'(\d+)', text)]


def parse_tile_file(image):
    file_name = image.replace('.', ' ').replace('_', ' ').split()

    if file_name[0] == 'HLS':
        sat_id = file_name[1]
        tile = file_name[2]
        time = file_name[3]

        for i in range(0, len(time), 4):
            year = time[:4]
            day = time[i:i + 4]
        # band_no = file_name[8]

    elif file_name[0] == 'MCD43A4':
        sat_id = file_name[0]
        year = file_name[1][1:5]
        day = file_name[1][5:]
        tile = file_name[2]
        # band_no = file_name[5]

    return sat_id, tile, year, day


def hdf2tif(hdf_path, tif_path, tiles, years, days):
    # Extracts tif subdatasets from HDF files stored in hdf_path directory
    # Extracted tif files are stored in tif_path directory
    if not os.path.exists(tif_path):
        print(tif_path)
        os.mkdir(tif_path)
        print("Directory ", tif_path, "created ")
    else:
        print("Directory ", tif_path, "already exists")

    file_lists = []
    for root, dirs, files in sorted(os.walk(hdf_path)):
        files.sort(key=natural_keys)
        file_lists.append(files)

    for files in file_lists:
        for file in files:
            if file.endswith(".hdf"):
                #sat_id, tile, year, day = parse_tile_file(file)

                # if tile in tiles and year in years and day in days:
                if True:
                    path = f"{hdf_path}/{file}"
                    print("Path is ________", path)
                    hdf_file = gdal.Open(path, gdal.GA_ReadOnly)
                    subdatasets = hdf_file.GetSubDatasets()

                    for i in range(0, len(subdatasets)):
                        subdataset_name = subdatasets[i][0]

                        if 'Grid:band' in subdataset_name:
                            band_name = subdataset_name.partition('Grid:band')[2]
                            # uncomment the line bellow if your data set is modis
                            # band_name = subdataset_name.partition('MOD_Grid_MCD15A3H:')[2]
                            data_type = gdal.GDT_Byte
                        else:
                            band_name = subdataset_name.partition('Grid:')[2]
                            # uncomment the line bellow if your data set is modis
                            # band_name = subdataset_name.partition('MOD_Grid_MCD15A3H:')[2]
                            data_type = gdal.GDT_Int16

                        band = gdal.Open(subdataset_name, gdal.GA_ReadOnly)
                        band_image = band.ReadAsArray()

                        path_out = f"{tif_path}/{file[:-4]}.{band_name}.tif"
                        tif_file = gdal.GetDriverByName('GTiff').Create(path_out, band.RasterXSize, band.RasterYSize, 1,
                                                                        data_type)

                        tif_file.SetGeoTransform(band.GetGeoTransform())
                        tif_file.SetProjection(band.GetProjection())
                        tif_file.GetRasterBand(1).WriteArray(band_image)
                        print("Save Done: %s" % path_out)
                        tif_file.FlushCache()
                        del tif_file


"""For HLS data """
# TILE = '12UYA'
# TILE = '12UXB'
# TILE = '12UYB'
# tiles = [TILE]
# years = ['2019']
# days = ['001', '008']
# hdf_path = os.path.join("data", 'hls_tmp', TILE, '2019', 'L30')
# tif_path = os.path.join("data", 'hls_tmp_tif', TILE, '2019', 'L30')
# hdf2tif(hdf_path=hdf_path, tif_path=tif_path, tiles=tiles, years=years,
#         days=days)
# #
# hdf_path = os.path.join("data", 'hls_tmp', TILE, '2019', 'S30')
# tif_path = os.path.join("data", 'hls_tmp_tif', TILE, '2019', 'S30')
# hdf2tif(hdf_path=hdf_path, tif_path=tif_path, tiles=tiles, years=years,
#         days=days)

"""For modis data"""
modis_path = os.path.join("data", 'Modis')
modis_path_tif = os.path.join("data", "Modis_tif")
hdf2tif(hdf_path=modis_path, tif_path=modis_path_tif, tiles=1, years=1, days=1)

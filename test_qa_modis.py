"""QC layer """
import os
from utilities.FileUtils import FileUtils
from utilities.qa_modis import make_LAI_FPAR_QC_mask
import rioxarray as rxr
from copy import copy

YEAR = {"2019"}


def main():
    fileutils = FileUtils()
    modis_dir = os.path.join("data", "Modis_tif", "1_original")

    for year in YEAR:
        rasters_path = os.path.join("data", "Modis_tif", "1_original", year)
        qc_path = os.path.join("data", "Modis_tif", "2_qa_masked", year)
        if not os.path.exists(qc_path):
            os.makedirs(qc_path)
        for root, dirs, files in os.walk(rasters_path):
            for file_mo in files:
                if '.tif' in file_mo:
                    sat_id_mo, tile_mo, year_mo, day_mo, band_no_mo = fileutils.parse_tile_file(file_mo)
                    if band_no_mo == 'FparLai_QC':
                        modis_qc_dir = os.path.join(modis_dir, year, file_mo)
                        print("modis dir ", modis_qc_dir)
                        modis_qc = rxr.open_rasterio(modis_qc_dir)
                        modis_qc_array = modis_qc.data[0]
                        filter_SCF_QG_0 = make_LAI_FPAR_QC_mask(modis_qc_array, parameter_name='SCF_QC', value=0)
                        for file in files:
                            sat_id, tile, year, day, band_no = fileutils.parse_tile_file(file)
                            if sat_id == sat_id_mo and tile == tile_mo and year == year_mo and day == day_mo and \
                                    band_no != 'FparLai_QC':
                                raster_path = os.path.join(root, file)
                                raster = rxr.open_rasterio(raster_path)
                                filtered_raster = copy(raster)
                                filtered_raster.data[0] = filtered_raster.data[0] * filter_SCF_QG_0
                                save_path = os.path.join(qc_path, file)
                                filtered_raster.rio.to_raster(save_path)
                                print("save path ", save_path)


if __name__ == '__main__':
    main()

from utilities.Preprocessing import Preprocessing
import os
import rasterio as rio
from matplotlib import pyplot
from utilities.Preprocessing import Preprocessing
import numpy as np


def main():
    qa_masked_path = os.path.join('data', 'hls', '2_qa_masked')
    sat_ids = ['L30', 'S30']
    tiles_name = ['T12UWC']
    tiles = ['12UWC']
    years = ['2019']
    days_num = np.arange(149, 181).tolist()
    days = list(map(str, days_num))
    prepare = Preprocessing()

    for tile in tiles:
        for year in years:
            for sat_id in sat_ids:
                hls_path = os.path.join('data', 'hls', '1_original', tile, year, sat_id)
                print(hls_path)
                # Quality Assessment layer masking
                # prepare.apply_qa_mask_hls(hls_path, qa_masked_path, sat_ids, tiles_name, years, days)


if __name__ == "__main__":
    main()

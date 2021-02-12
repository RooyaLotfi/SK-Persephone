import os
import numpy as np
from utilities.Preprocessing import Preprocessing

# TODO: Fix merge so that it would get all the tiles based on ID
if __name__ == '__main__':
    qa_masked_path = os.path.join('data', '2_qa_masked')
    sat_ids = ['L30']
    tiles = ['T12UWA', 'T12UWB', 'T12UWC']
    years = ['2019']
    days_num = np.arange(154, 157).tolist()
    days = list(map(str, days_num))
    prepare = Preprocessing()
    region_name = '_'.join(tiles)
    print(region_name)
    for sat_id in sat_ids:
        for year in years:
            for day in days:
                print(prepare.merge(qa_masked_path, region_name, sat_id, tiles, year, day))

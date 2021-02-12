import numpy as np
import os
from Utils.merge import merge

qa_mask_path = os.path.join("data", "Modis_tif", "2_qa_masked", "2019")
dir_path = os.path.join("data", "Modis_tif", "3_merged")
sat_ids = 'MCD15A3H'
tiles = ['h10v03', 'h10v04', 'h11v03', 'h11v04', 'h12v03']
years = ['2019']
days_num = np.arange(150, 207).tolist()
days = list(map(str, days_num))
region_name = '-'.join(tiles)
for year in years:
    for day in days:
        # print(f"{year},{day}")
        merge(path=qa_mask_path, dir_path=dir_path, region_name=region_name, sat_id=sat_ids, tiles=tiles, year=year,
              day=day)

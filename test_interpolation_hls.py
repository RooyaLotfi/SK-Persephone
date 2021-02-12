import os
import rioxarray as rxr
import numpy as np
from utilities.FileUtils import FileUtils
from utilities.geotiff_utils import interpolate
import matplotlib.pyplot as plt
import numpy.ma as ma
import math


def main():
    hls_path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone/data/qa_masked/S30/2019/T12UWB/232",
                            "hls_masked_crop_146")
    f = FileUtils()
    interpolate_test_dict = {}
    band_number = 1
    for root, dirs, files in os.walk(hls_path):
        for file in files:
            if '.tif' in file:
                image_path = os.path.join(root, file)
                sat_id, tile, year, day, band_no, crop_type = f.parse_preprocessed_file(file)
                if int(band_no) == band_number:
                    raster = rxr.open_rasterio(image_path)
                    data = raster.data
                    # print("in valid stuff ",ma.masked_invalid(data))
                    data[data == -1000] = np.NaN
                    interpolate_test_dict[int(day)] = data

    interpolated_data = interpolate(interpolate_test_dict, kind='linear', as_dict=True)  # interpolate the data linearly

    path = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone/data/qa_masked/S30/2019/T12UWB/232",
                        "hls_masked_crop_146_interpolated")
    if not os.path.exists(path):
        os.makedirs(path)
    print()
    print("number of nans for 151 ", np.sum(np.isnan(interpolated_data[151])))
    print("number of nans for 155 ", np.sum(np.isnan(interpolated_data[155])))
    print("number of nans for 161 ", np.sum(np.isnan(interpolated_data[161])))
    print("number of nans for 162 ", np.sum(np.isnan(interpolated_data[162])))
    print("number of nans for 170 ", np.sum(np.isnan(interpolated_data[170])))
    print("number of nans for 179 ", np.sum(np.isnan(interpolated_data[179])))

    # For saving in the file
    raster_temp = raster
    # print(interpolated_data)
    for key in interpolated_data:
        print("key is ", key)
        print("value is ",interpolated_data[key])
        raster_temp.data = interpolated_data[key]
        path_save = os.path.join(path, str(key)+'.tif')
        raster_temp.rio.to_raster(path_save)
    #     #print("value is ", value)
    #     print("******************************************************* saved ", path_save)

    # plt.imshow(interpolated_data[166][0][1150:1250, 370:450])
    # print("some of nan values of interpolated data ", np.sum(interpolated_data[166][0][1150:1250, 370:450] == -1000))
    # plt.show()


if __name__ == '__main__':
    main()

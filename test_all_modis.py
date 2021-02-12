import test_webscrap_modis
import test_hdf2tif_modis
import os
import test_qa_modis
import numpy as np
import test_merge_modis
import test_crop_modis_sk
import test_resample_modis
import test_crop_modis
import test_filter_crops_modis

# downloading modis data for year 2019 and days 150-190 same as hls and for tiles that SK is in there
test_webscrap_modis.main()
# converting modis data from hdf to tif file
test_hdf2tif_modis.main()
# Apply QA layer
test_qa_modis.main()
test_merge_modis.main()
test_crop_modis_sk.main()
test_resample_modis.main()
test_crop_modis.main()
# test_filter_crops_modis.main()

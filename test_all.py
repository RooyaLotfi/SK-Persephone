import test_webscrap
import test_voronoi
import test_merge
import test_QA
import test_crop_module
import test_filter_crops
import test_resample

# test_voronoi.main()
test_webcscrap.main()
test_QA.main()
test_merge.main()
test_resample.main()

test_crop_module.main()
test_filter_crops.main()

# weather : daymet - 1km^2 (upsample)
# calculate vegitation index (hls: 30m)
# physical parameters (lai, chlorophyl)
# modis products (up sample) 500m^2,
# soil moisture from smap
# soil properties can get it from soil grid (250m^2)
# crop
# crop AAFC
# mask : hls, daymet, physical parameters
# use data cube for region of interest
# stacking

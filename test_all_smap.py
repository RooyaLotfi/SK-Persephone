import test_crop_SMAP
import test_resample_soilgrid
import test_crop_skboundary
import test_resample_smap
import test_filter_crops_smap

# The smap files were very large to upsample so I first cropped it to the extent of saskatchewan boudaries and I
# defined a buffer around the sk boundary file because If we crop it and then upsample it some pixels will be missing
# so we defined a buffer to keep all the pixels
test_crop_skboundary.main()
test_resample_smap.main()
test_crop_SMAP.main()
test_filter_crops_smap.main()

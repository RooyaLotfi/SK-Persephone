import os


def resample(path_in, in_crs, path_out, out_crs='EPSG:32613', target_resolution=30, resampling_method='bilinear'):
    command = f"gdalwarp -s_srs {in_crs} -t_srs {out_crs} -tr {target_resolution}.0 {target_resolution}.0 -r {resampling_method} -tap -of GTiff {path_in} {path_out}"
    output = os.system(command)
    print(output)

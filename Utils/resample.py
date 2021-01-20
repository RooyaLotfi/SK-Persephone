import os

def resample(path_in, path_out, target_resolution, resampling_method):
    command = f"gdalwarp -s_srs EPSG:4326 -t_srs EPSG:32613 -tr {target_resolution}.0 {target_resolution}.0 -r {resampling_method} -tap -of GTiff {path_in} {path_out}"
    output = os.system(command)
    print(output)
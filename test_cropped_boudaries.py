import os
import geopandas as gpd
import matplotlib.pyplot as plt

path = os.path.join("polygon_EPSG_32613", "cropped_by_CA_boundaries", "polygon.shp")
gdframe = gpd.read_file(path)

print(gdframe)
gdframe.plot()
plt.show()

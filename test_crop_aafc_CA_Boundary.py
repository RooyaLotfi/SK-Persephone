import os
import geopandas as gpd
from preprocessing.voronoi import write_geodf
import rasterio as rio
import matplotlib.pyplot as plt
import rioxarray as rxr

CRS = "epsg:32613"

# Import all of your data at the top of your notebook to keep things organized.
sk_bouondary_path = os.path.join("Canada-Boundaries", "Saskatchewan", "sk_boundary.shp")
sk_boundary = gpd.read_file(sk_bouondary_path)

sk_boundary = sk_boundary.to_crs(CRS)

sk_AAFC_path = os.path.join("AAFC_data", "aci_2019_sk_v1", "aci_2019_sk.tif")

sk_AAFC = rio.open(sk_AAFC_path)

f, ax = plt.subplots(figsize=(10, 5))
plt.plot(sk_AAFC.read(1))
ax.set(title="AAFC")
plt.show()


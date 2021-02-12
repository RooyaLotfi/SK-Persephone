"""
Cropping Voronoi boundaries by the SK boundary that is defined in Canada-Boundaries file
Visualizing the SK boundary, Voronoi RM Polygons, Cropped Voronoi
"""

import os
import geopandas as gpd
from preprocessing.voronoi import write_geodf
import matplotlib.pyplot as plt

sk_RM_path = os.path.join("polygon_EPSG_32613", "polygon.shp")
sk_RM_voronoi = gpd.read_file(sk_RM_path)

CRS = "epsg:32613"

# Import all of your data at the top of your notebook to keep things organized.
sk_bouondary_path = os.path.join("Canada-Boundaries", "Saskatchewan", "sk_boundary.shp")
sk_boundary = gpd.read_file(sk_bouondary_path)

sk_boundary = sk_boundary.to_crs(CRS)
# Clip the data using GeoPandas clip
points_clip = gpd.clip(sk_RM_voronoi, sk_boundary)

# Plot the data
fig, ax = plt.subplots(figsize=(12, 8))
sk_bound = sk_boundary.plot(ax=ax, color='mistyrose', label="Saskatchewan Boundary")
vor_bound = sk_RM_voronoi.plot(ax=ax, color='skyblue', label="Voronoi Boundaries")
overlap_bound = points_clip.plot(ax=ax, color='royalblue', label="Overlap")

plt.suptitle('Overlapping region between the voronoi regions and Saskatchewan boundary', fontsize=20, color='dimgray')

plt.xlabel("Longitude", fontsize=15, color='lightslategray')
plt.ylabel("Latitude", fontsize=15, color='lightslategray')
plt.savefig("figures/cropped_voronoi_regions.png")
plt.show()

write_geodf(points_clip, "polygon_EPSG_32613/cropped_by_CA_boundaries", "polygon.shp")

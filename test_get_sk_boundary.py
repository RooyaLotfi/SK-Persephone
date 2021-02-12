"""
Extracting SK boundaries from all boundaries in Canada-Boundaries file (whole canada province boundaries)
"""

import os
import geopandas as gpd
from preprocessing.voronoi import write_geodf

SK_RM_path = os.path.join("polygon_EPSG_32613", "polygon.shp")
SK_RM = gpd.read_file(SK_RM_path)

# Import all of your data at the top of your notebook to keep things organized.
country_boundary_ca_path = os.path.join("Canada-Boundaries", "gpr_000b11a_e.shp")
country_boundary_ca = gpd.read_file(country_boundary_ca_path)

gpd_sk = country_boundary_ca.loc[country_boundary_ca["PRENAME"]=="Saskatchewan"]
path = os.path.join("Canada-Boundaries", "Saskatchewan")
write_geodf(gpd_sk, path, "sk_boundary.shp")

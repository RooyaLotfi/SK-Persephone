"""
Module for creating diagram Voronoi for a geodataframe that has centroid datapoints

"""

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
from scipy.spatial import distance
from scipy.spatial import Voronoi
import pandas as pd
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

eps = sys.float_info.epsilon


def _in_box(coords, bounding_box):
    """
    Checks if all points in coords are within a rectangular area. The rectangle axis are
    defined by bounding_box.
    :param coords: Numpy array of coordinates in 2d space
    :param bounding_box: Numpy array of size 4: [minimum x, maximum x, minimum y, maximum y]
    :return: A numpy array of True and False. True/False: the point is/isn't inside rectangle
    """
    return np.logical_and(np.logical_and(bounding_box[0] <= coords[:, 0],
                                         coords[:, 0] <= bounding_box[1]),
                          np.logical_and(bounding_box[2] <= coords[:, 1],
                                         coords[:, 1] <= bounding_box[3]))


def _centroid_region(vertices):
    """
    Finding centroid of a polygon. If vertices are: v1=[a,b,c] v2=[a2,b2,c2] v3=[a3,b3,c3]
    We need centroids to connect voronoi regions to voronoi centroids
    center is: [(a1+a2+a3)/3,(b1+b2+b3)/3,(c1+c2+c3)/3]
    :param vertices: polygon vertices
    :return: center of the polygon
    """
    # Polygon's signed area
    area = 0
    # Centroid's x
    c_x = 0
    # Centroid's y
    c_y = 0
    for i in range(0, len(vertices) - 1):
        s = (vertices[i, 0] * vertices[i + 1, 1] - vertices[i + 1, 0] * vertices[i, 1])
        area = area + s
        c_x = c_x + (vertices[i, 0] + vertices[i + 1, 0]) * s
        c_y = c_y + (vertices[i, 1] + vertices[i + 1, 1]) * s
    area = 0.5 * area
    c_x = (1.0 / (6.0 * area)) * c_x
    c_y = (1.0 / (6.0 * area)) * c_y
    return np.array([[c_x, c_y]])


def _find_idx_point(vertices, point):
    """
    Detects the vertex index in which the point is inside it
    :param vertices: a set of polygon vertices
    :param point: a point inside one of the vertices
    :return: the index of the polygon that the vertex is inside it's area
    """
    centroid = _centroid_region(vertices)
    # second, if we have more than one point inside polygon we need to remove one of them
    # noinspection DuplicatedCode
    dist = distance.cdist(point, centroid, 'euclidean')
    center_index = np.argmin(dist)
    # first the point needs to be inside rectangle boundary that is defined by minimum
    # and maximum over x,y of all vertices
    bounding_box = [np.min(vertices[:, 0]), np.max(vertices[:, 0]),
                    np.min(vertices[:, 1]), np.max(vertices[:, 1])]

    is_in_polygon = np.logical_and(np.logical_and(bounding_box[0] <= point[center_index, 0],
                                                  point[center_index, 0] <= bounding_box[1]),
                                   np.logical_and(bounding_box[2] <= point[center_index, 1],
                                                  point[center_index, 1] <= bounding_box[3]))
    if is_in_polygon:
        return center_index
    print("The point is not inside any of the vertices!")
    return 0


def _join_center2boundary(centroids, vor, cols, is_shapely_bndry):
    """
    Gets a geodataframe (centroids) and it's corresponding voronoi object (vor) and
    the name of columns for creating dataframe. Finds which centroid (centroids) belongs
    to which voronoi vertices (vor) and makes a column of each centroid and it's
    corresponding region.
    :param centroids: a geodataframe of centroid
    :param vor: voronoi diagram that is gotten from coordinates in the centroids geodataframe
    :param cols: the column names for [voronoi center, voronoi boundary vertices]
    :param is_shapely_bndry: True: convert numpy voronoi boundary to shapely.
    False: keep the boundary as numpy array
    :return: dataframe with columns [voronoi center, voronoi boundary vertices]
    """
    boundary_df = pd.DataFrame(columns=cols)
    for region in vor.filtered_regions:
        vertices = vor.vertices[region + [region[0]], :]
        idx_vor_center = _find_idx_point(vertices, centroids)
        centroid = centroids[idx_vor_center]
        if is_shapely_bndry:
            boundary_df = boundary_df.append({cols[0]: Polygon(vertices),
                                              cols[1]: np.array(centroid)},
                                             ignore_index=True)
        else:
            boundary_df = boundary_df.append({cols[0]: vertices, cols[1]: centroid},
                                             ignore_index=True)
    return boundary_df


def _find_id(dframe, point):
    """
    The geodataframes have a column named id. This function finds the id of of a
    given point in the dataframe example: array([1,2,3]). dataframe["id":{12},
    "geometry":{[1,2,3]}]. id = 12
    :param dframe: geopandas dataframe with columns ["id", "geometry"]
    :param point: a numpy array point
    :return: the "id" corresponding to "geometry" = point
    """
    # Geometry column in geopandas dataframe is in form of Shapely object so we need to convert
    # numpy to shapely to compare if the values are equal
    loc = dframe.loc[dframe["geometry"] == Point(point)]
    center_id = pd.to_numeric(loc["id"].iloc[0])
    return center_id


def _add_id(bdry_cent_df, dframe):
    """
    Iterates through each row of bdry_cent_df and extracts value of "vor_centroid" and calls
     function _find_id to find "id" value corresponding to "vor_centroid" and stores them in a list
    :param bdry_cent_df: a dataframe of ["vor_centroid" or voronoi centroid, voronoi boundary]
    :param dframe: Geodataframe with ["id", "geometry" or voronoi centroid]
    :return: New dataframe with appended id
    """
    id_list = []
    # extract list of IDs corresponding to each vor_centroid in boundary_df
    for index, row in bdry_cent_df.iterrows():
        id_list.append(_find_id(dframe, row["vor_centroid"]))
    # add list of all IDs for each row in bdry_cent_df
    bdry_cent_df['id'] = id_list
    return bdry_cent_df


def create_geodf(dframe, cols, col_num):
    """
    gets a dataframe and returns a geodataframe
    :param dframe: the dataframe that we want to transform it to geodataframe
    :param cols: the columns in dframe that you wish to keep
    :param col_num: which column should be geometry column (shapely object)
    :return: geodataframe with cols and geometry column
    """
    new_dframe = dframe[cols]
    get_col = cols[col_num]
    new_dframe.set_geometry(get_col)
    new_dframe.rename(columns={get_col: "geometry"}, inplace=True)
    new_df_gpd = gpd.GeoDataFrame(new_dframe, crs=dframe.crs, geometry='geometry')
    return new_df_gpd


def get_voronoi_df(centroid_df, factor):
    """
    Gets a dataframe in form of geopandas and does the following steps:
     1. extract centroid ("geometry") column
     2. calculates diagram voronoi and gets vor object that has voronoi regions and vertices
     3. defines a new dataframe with 2 columns ('vor_boundary', 'vor_centroid') to match voronoi
     boundaries to voronoi centroids
     4. the ID's that were lost while making voronoi diagram will be added to the dataframe that
      was created in step 3
     5. now we merge two dataframes by column "id"
    :param centroid_df: Geodataframe of centroids
    :param factor: bounding factor to create diagram voronoi
    :return: centroid_df with an added column:
    ['vor_boundary' or vertices of voronoi boundary, 'vor_centroid' or centroids of voronoi diagram]
    """
    # Extract centroids
    centroids = extract_coords(centroid_df)
    # bounding factor = how we should mirror the points
    vor = voronoi(centroids, factor)
    # storing vertices and it's centroid in a pandas dataframe:
    columns = ('vor_boundary', 'vor_centroid')
    centroid_boundary_df = _join_center2boundary(centroids, vor, columns, True)
    centroid_boundary_df = _add_id(centroid_boundary_df, centroid_df)
    # merge sk_centers dataframe with boundary dataframe based on equal values in geometry
    merged_df = pd.merge(centroid_df, centroid_boundary_df, on="id")
    return merged_df


def write_geodf(gdframe, path):
    """
    gets a geodataframe and stores in in a shapefile in the path directory
    :param gdframe: geodataframe to be saved in a shape file
    :param path: path to save the geodataframe
    """
    save_path = os.path.join(path, "boundaries.shp")
    if os.path.exists(path):
        gdframe.to_file(save_path)
    else:
        os.mkdir(path)
        gdframe.to_file(save_path)


def extract_coords(dframe):
    """
    Extracts geometry column from pandas dataframe and transfrom it from a shapely file
     to a numpy array
    :param dframe: Geopandas dataframe
    :return: a numpy array of "geometry" column
    """
    coords = np.empty((0, 2), dtype=float)
    for index, row in dframe.iterrows():
        coords = np.append(coords, np.array(list(row['geometry'].coords)), axis=0)
    return coords


def voronoi(coords, factor):
    """
    Calculate diagram voronoi for a set of datapoints (coords) with a bounding factor
    :param coords: Voronoi diagram will be applied to coords points
    :param factor: Factor by which we want to limit our diagram
    :return: An object of Voronoi type with properties attributes like regions,ets.
    """
    bounding_box = np.array(
        [min(coords[:, 0]) - factor, max(coords[:, 0]) + factor,
         min(coords[:, 1]) - factor,
         max(coords[:, 1]) + factor])  # [x_min, x_max, y_min, y_max] in the map

    # Select coords inside the bounding box
    i = _in_box(coords, bounding_box)
    # Mirror points
    points_center = coords[i, :]
    points_left = np.copy(points_center)
    points_left[:, 0] = bounding_box[0] - (points_left[:, 0] - bounding_box[0])
    points_right = np.copy(points_center)
    points_right[:, 0] = bounding_box[1] + (bounding_box[1] - points_right[:, 0])
    points_down = np.copy(points_center)
    points_down[:, 1] = bounding_box[2] - (points_down[:, 1] - bounding_box[2])
    points_up = np.copy(points_center)
    points_up[:, 1] = bounding_box[3] + (bounding_box[3] - points_up[:, 1])
    points = np.append(points_center,
                       np.append(np.append(points_left,
                                           points_right,
                                           axis=0),
                                 np.append(points_down,
                                           points_up,
                                           axis=0),
                                 axis=0),
                       axis=0)

    plt.show()
    # Compute Voronoi
    vor = Voronoi(points)
    # Filter regions
    regions = []
    for region in vor.regions:
        flag = True
        for index in region:
            if index == -1:
                flag = False
                break
            x = vor.vertices[index, 0]
            y = vor.vertices[index, 1]
            if not (bounding_box[0] - eps <= x <= bounding_box[1] + eps and bounding_box[2]
                    - eps <= y <= bounding_box[3] + eps):
                flag = False
                break
        if region != [] and flag:
            regions.append(region)
    vor.filtered_points = points_center
    vor.filtered_regions = regions
    return vor


def plot_voronoi(centroid, vor):
    """
    Uses the diagram voronoi object from scipy.spatial library to plot each voronoi region.
    :param centroid: Center of Diagram Voronoi
    :param vor: Is Diagram Voronoi object corresponding to centroid points
    """
    # Plot initial centers that we read from file
    plt.subplots(figsize=(10, 10))
    centroid_points, = plt.plot(centroid[:, 0], centroid[:, 1], 'o', color='salmon')

    # Plot ridges
    for region in vor.filtered_regions:
        vertices = vor.vertices[region + [region[0]], :]
        region_edge, = plt.plot(vertices[:, 0], vertices[:, 1], '-', color='lightgrey')

    # Plot ridges points
    for region in vor.filtered_regions:
        vertices = vor.vertices[region, :]
        vertice_points, = plt.plot(vertices[:, 0], vertices[:, 1], '*', color='steelblue')

    plt.legend([vertice_points, region_edge, centroid_points],
               ["Voronoi Vertices ", "Voronoi Boundaries", "Centroids"])
    plt.plot()
    plt.suptitle('Diagram Voronoi', fontsize=20, color='dimgray')
    plt.xlabel("Latitude", fontsize=15, color='lightslategray')
    plt.ylabel("Longitude", fontsize=15, color='lightslategray')
    path = os.path.join("figures")
    save_path = os.path.join("figures", "voronoi.png")
    if os.path.exists(path):
        plt.savefig(save_path)
    else:
        os.mkdir(path)
        plt.savefig(save_path)
    plt.show()

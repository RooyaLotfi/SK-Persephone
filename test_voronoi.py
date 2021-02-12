from preprocessing.voronoi import *
import geopandas as gpd
import matplotlib.pyplot as plt


def main():
    centroids_path = os.path.join("RM_centroids", "original", "centroid.shp")
    sk_centers_df = gpd.read_file(centroids_path)

    # I changed this CRS because Brock mentioned I need to do so
    CRS = "EPSG:32613"
    BOUNDING_FACTOR = 20000
    # changing CRS
    sk_centers_df = sk_centers_df.to_crs(CRS)

    # A dataframe with boundary, centroid and id
    df_vor = get_voronoi_df(sk_centers_df, BOUNDING_FACTOR)
    _output_crs = CRS.replace(":", "_")
    polygon_output_path = os.path.join("RM_centroids", "polygon_" + _output_crs)
    centroids_output_path = os.path.join("RM_centroids", "centroid_" + _output_crs)

    geo_col_num = 3
    columns = ['id', 'xcoord', 'ycoord', 'vor_boundary']

    new_df = create_geodf(df_vor, columns, geo_col_num)
    write_geodf(new_df, polygon_output_path, "polygon.shp")
    write_geodf(sk_centers_df, centroids_output_path, "centroid.shp")

    # Extract RM_centroids
    center_coord = extract_coords(sk_centers_df)
    voronoi_diagram = voronoi(center_coord, BOUNDING_FACTOR)
    plot_voronoi(center_coord, voronoi_diagram)

    # df_polygon = gpd.read_file(polygon_output_path)
    # df_centroid = gpd.read_file(centroids_output_path)
    # df_polygon.plot()
    # df_centroid.plot()
    # plt.show()


if __name__ == '__main__':
    main()

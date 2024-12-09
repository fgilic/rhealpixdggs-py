import shapely
import geopandas as gpd
import pandas as pd
from dask.distributed import Client
import numpy as np
from itertools import chain

if __name__ == "__main__":
    client = Client(processes=True, n_workers=16)

    # cleaning from duplicates (duplicated points on boundaries between strips)
    # gdf = gpd.read_parquet("cell_nuclei_built-up_raw.parquet")
    # gdf_clean = gdf.drop_duplicates("suid")
    #
    # print(len(gdf_clean))
    #
    # gdf_clean.to_parquet("cell_nuclei_built-up_no_duplicates.parquet")

    points = gpd.read_parquet("cell_nuclei_built-up_no_duplicates.parquet")
    polygons = gpd.read_file("reporting_units.geojson", crs="EPSG:4326")
    points_list = points.to_dict("records")
    point_groups = np.array_split(points_list, 128)

    def func(point_group):
        points_to_return = []
        for point in point_group:
            for index, polygon in polygons.iterrows():
                country_name = polygon["name_en"]
                if shapely.contains(polygon["geometry"], point["geometry"]):
                    point["country"] = country_name
                    break
                else:
                    point["country"] = None
            points_to_return.append(point)
        return points_to_return

    futures = client.map(func, point_groups)
    results = client.gather(futures)
    results = list(chain.from_iterable(results))
    gpd.GeoDataFrame(results, crs="EPSG:4326").to_parquet(
        "cell_nuclei_overlaid_reporting_units.parquet"
    )

    #     points_within = gpd.GeoDataFrame(points_within)
    #     points_within.to_parquet(f"nucleus_built-up_{country}.parquet", crs="EPSG:4326")
    #     print(f"Total number of nuclei for {country}: {len(points_within)}")
    #     print(f"Number of nuclei for 2010 for {country}: {points_count_2010}")
    #     print(f"Number of nuclei for 2020 for {country}: {points_count_2020}\n")
    #     polygon["nuclei_count_built-up_2010"] = points_count_2010
    #     polygon["nuclei_count_built-up_2020"] = points_count_2020
    #     polygons_output.append(polygon)
    # gpd.GeoDataFrame(polygons_output).to_parquet("reporting_units_count_nuclei.parquet", crs="EPSG:4326")
    # print(f"Total output number of written nuclei: {total_output_number_of_nuclei}")
    # print(f"Total input number of nuclei: {len(points)}")

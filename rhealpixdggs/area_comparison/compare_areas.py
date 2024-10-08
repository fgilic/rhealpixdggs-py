import pandas as pd
import geopandas as gpd

# input: - two parquet files with data calculated using code from tag 0.5.4 and using new code
#        - one FlatGeobuf file with cell geometry (low density vertices on polygon edges - for visualization)
# output: FltGeobuf with cell geometry (low density) and attribute data for each cell (area obtained by new and old
# code, theoretical area, differences


level = 0

# theoretical area from tag 0.5.4 is wrong (issue #18)
grid_old = pd.read_parquet(f"grid_{level}_0-5-4.parquet")[
    [
        "Cell_suid",
        "Cell_region",
        "Cell_shape",
        "Cell_nucleus",
        "Calculated_area_of_cell",
    ]
]
grid_old = grid_old.rename(
    columns={
        "Cell_nucleus": "Cell_nucleus_old",
        "Calculated_area_of_cell": "Calculated_area_old",
    }
)
grid_new = pd.read_parquet(f"grid_{level}.parquet")[
    ["Cell_suid", "Theoretical_area_of_cell", "Cell_nucleus", "Calculated_area_of_cell"]
]
grid_new = grid_new.rename(
    columns={
        "Cell_nucleus": "Cell_nucleus_new",
        "Calculated_area_of_cell": "Calculated_area_new",
    }
)
grid_geometry_lq = gpd.read_file(f"grid_{level}_LQ.fgb")[["Cell_suid", "geometry"]]

grid_combined = pd.merge(grid_old, grid_new, on="Cell_suid")
grid_combined["Area_difference_old"] = (
    grid_combined["Theoretical_area_of_cell"] - grid_combined["Calculated_area_old"]
).abs()
grid_combined["Area_difference_new"] = (
    grid_combined["Theoretical_area_of_cell"] - grid_combined["Calculated_area_new"]
).abs()

print(grid_combined["Area_difference_old"].describe())
print(grid_combined["Area_difference_new"].describe())

# grid_combined["nucleus_geom"] = grid_combined["Cell_nucleus_new"].apply(lambda coords : Point(ast.literal_eval(coords)))
grid_combined = pd.merge(grid_combined, grid_geometry_lq, on="Cell_suid")

grid_gdf = gpd.GeoDataFrame(data=grid_combined, geometry="geometry", crs="EPSG:4326")
grid_gdf.to_file(f"grid_{level}_combined.fgb", driver="FlatGeobuf")

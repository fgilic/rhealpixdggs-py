import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import geopandas as gpd
import pandas as pd
import pyproj
import matplotlib as mpl
from rhealpixdggs.qpix_dggs import QPixDGGS
from rhealpixdggs.dggs import RHEALPixDGGS
from rhealpixdggs import ellipsoids
import shapely
import mapclassify

n_side = 3
cells_data = pd.read_csv("..\\compactness_res-5_nside-3.csv")
cells_data_054 = pd.read_csv("../areaRHPcells_res-5_nside-3_054.csv")
plot = "rhp" # "rhp" or "qpix"

cells_data_concat = pd.concat([cells_data["rCell_suid"], abs(cells_data["Theoretical_area_of_rcell"]-cells_data["Calculated_area_of_rcell"]), abs(cells_data_054["Theoretical_area_of_rcell"]-cells_data_054["Calculated_area_of_rcell"])], axis=1)
cells_data_concat = cells_data_concat.rename(columns={0:"err_old",1:"err_new"})

errors_concat = pd.concat(
    [
        cells_data_concat.iloc[:,1],
        cells_data_concat.iloc[:,2],
    ]
)
vmin = errors_concat.min()
vmax = errors_concat.max()

cells_data_concat = cells_data_concat.to_dict(orient="records")



a_wgs84 = pyproj.list.get_ellps_map()["WGS84"]["a"]
f_wgs84 = 1 / pyproj.list.get_ellps_map()["WGS84"]["rf"]
wgs84_crs = pyproj.CRS.from_epsg(4326)
wgs84_ellipsoid = ellipsoids.Ellipsoid(a=a_wgs84, f=f_wgs84)
rdggs = RHEALPixDGGS(
    ellipsoid=wgs84_ellipsoid, north_square=0, south_square=0, N_side=n_side
)



cell_data_with_geometry = []
for cell_data in cells_data_concat:
    cell_suid_str = cell_data["rCell_suid"]
    cell_suid = []
    for i in cell_suid_str:
        try:
            cell_suid.append(int(i))
        except:
            cell_suid.append(i)
    cell_suid = tuple(cell_suid)
    if plot == "rhp":
        cell = rdggs.cell(cell_suid)
        # exclude cells that intersect meridian with lambda = +-180°
        if cell.intersects_meridian(lam=-180.0) or cell.intersects_meridian(lam=-180.0):
            continue
        cell_boundary_points = cell.boundary(plane=False, interior=True)
        cell_geom = shapely.Polygon(cell_boundary_points)



    cell_data["geometry"] = cell_geom
    cell_data_with_geometry.append(cell_data)



fig = plt.figure(figsize=(20, 10), frameon=False)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
# ax.axis("off")
ax.set_global()
ax.coastlines(linewidth=1.75)
# ax.stock_img()

grid_gdf = gpd.read_file("grid_30deg.geojson")
grid_gdf = grid_gdf.to_crs(ccrs.Robinson())

# world_gdf = gpd.read_file("coastlines_simplified.geojson")
# world_gdf = world_gdf.to_crs(ccrs.Robinson())
cells_gdf = gpd.GeoDataFrame(cell_data_with_geometry, crs=wgs84_crs)
cells_gdf = cells_gdf.to_crs(ccrs.Robinson())


errors_concat = pd.concat(
    [
        cells_gdf.iloc[:,1],
        cells_gdf.iloc[:,2],
    ]
)
vmin = errors_concat.min()
vmax = errors_concat.max()


# mask = (cells_gdf["rCell_c_sp"] < 0.92) & (cells_gdf["qCell_c_sp"] < 0.92)
# cells_gdf_masked = cells_gdf[mask]


# rhp_qpix_csp_concat = pd.concat(
#     [
#         cells_gdf_masked["qCell_c_sp"],
#         cells_gdf_masked["rCell_c_sp"],
#     ])
#
# mean_compactness = rhp_qpix_csp_concat.mean()
# cells_gdf["rCell_c_sp_diff"] = abs(cells_gdf["rCell_c_sp"] - mean_compactness)
# cells_gdf["qCell_c_sp_diff"] = abs(cells_gdf["qCell_c_sp"] - mean_compactness)

# mean_rhp_compactness = cells_gdf_masked["qCell_c_sp"].mean()
# mean_qpix_compactness = cells_gdf_masked["rCell_c_sp"].mean()
# cells_gdf["rCell_c_sp_diff"] = abs(cells_gdf["rCell_c_sp"] - mean_rhp_compactness)
# cells_gdf["qCell_c_sp_diff"] = abs(cells_gdf["qCell_c_sp"] - mean_qpix_compactness)

# # cells_gdf_masked = cells_gdf[mask]
# rhp_qpix_csp_concat_diff = pd.concat(
#     [
#         cells_gdf_masked["qCell_c_sp_diff"],
#         cells_gdf_masked["rCell_c_sp_diff"],
#     ])
# vmin = rhp_qpix_csp_concat_diff.min()
# vmax = rhp_qpix_csp_concat_diff.max()


# if plot == "rhp":
#     mask = cells_gdf["rCell_c_sp"] < 0.92
#     cells_gdf = cells_gdf[mask]
# elif plot == "qpix":
#     mask = cells_gdf["qCell_c_sp"] < 0.92
#     cells_gdf = cells_gdf[mask]

grid_gdf.plot(ax=ax, edgecolor="black", linewidth=2.5, zorder=3)

# values in (value, color) defined by bin edges with Quantiles, 6 bins
# mapclassify.classify(rhp_qpix_csp_concat, scheme="Quantiles", k=6)
# and normalized from 0.0 to 1.0
cmap_spectral = mpl.colors.LinearSegmentedColormap.from_list(
    N=500000,
    name="spectral",
    colors=[
        (0, "#006837"),
        (0.0000486195, "#4bb05c"),
        (0.0001068702, "#b7e075"),
        (0.0007345004, "#fffebe"),
        (0.0042498136, "#fdbf6f"),
        (0.0087190610, "#ea5739"),
        (1, "#a50026"),
    ],
)


cells_gdf.plot(
    ax=ax,
    linewidth=0.4,
    column="err_new",
    zorder=1,
    legend=True,
    edgecolor="face",
    cmap=cmap_spectral,
    vmin=vmin,
    vmax=vmax,
)
print("RHP_nside3_res5")
print(vmin)
print(vmax)

plt.savefig("area_error_new.png",  dpi=900)

# # START only plot compactness for nside2, resolution6 and nside3, resolution4
# cmap_blueish = mpl.colors.LinearSegmentedColormap.from_list(
#     N=500000,
#     name="blueish",
#     colors=[
#         (0, "#0096A0"),
#         (0.17, "#26A6AF"),
#         (0.32, "#67C1C9"),
#         (0.50, "#7CCAD1"),
#         (0.67, "#A3DBE1"),
#         (0.83, "#CEEEF2"),
#         (1, "#F0FCFF"),
#     ],
# )
#
# cmap_pink = mpl.colors.LinearSegmentedColormap.from_list(
#     N=500000,
#     name="pink",
#     colors=[
#         (0, "#A00083"),
#         (0.17, "#AC2193"),
#         (0.32, "#BB48A7"),
#         (0.50, "#C769B7"),
#         (0.67, "#D388C6"),
#         (0.83, "#DDA2D3"),
#         (1, "#F7E9F6"),
#     ],
# )
#
# # vmin = 0.73 # for QPix
# # vmax = 0.80 # for QPix
#
# vmin = 0.58 # for rHP
# vmax = 0.79 # for rHP
#
# cells_gdf.plot(
#     ax=ax,
#     linewidth=0.4,
#     column="rCell_c_sp",
#     zorder=1,
#     legend=True,
#     edgecolor="face",
#     cmap=cmap_pink.reversed(),
#     vmin=vmin,
#     vmax=vmax,
# )
# # plt.savefig("compactness_nside2_res6.svg")
# plt.savefig("RHP_compactness_nside3_res5.png", dpi=900)
# # END

plt.show()


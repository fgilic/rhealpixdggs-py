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

n_side = 3
cells_data = pd.read_csv("..\\compactness_res-0_nside-3.csv").to_dict(orient='records')
plot = "rhp" # "rhp" or "qpix"

a_wgs84 = pyproj.list.get_ellps_map()["WGS84"]["a"]
f_wgs84 = 1 / pyproj.list.get_ellps_map()["WGS84"]["rf"]
wgs84_crs = pyproj.CRS.from_epsg(4326)
wgs84_ellipsoid = ellipsoids.Ellipsoid(a=a_wgs84, f=f_wgs84)
rdggs = RHEALPixDGGS(
    ellipsoid=wgs84_ellipsoid, north_square=0, south_square=0, N_side=n_side
)
qdggs = QPixDGGS(
    ellipsoid=wgs84_ellipsoid, north_square=0, south_square=0, N_side=n_side
)
cell_data_with_geometry = []

def check_antimeridian(cell):
    # only for QPix DGGS, and only when north_south and south_square are (0,0)
    boundary = cell.boundary(plane=True)
    cell_geom = shapely.Polygon(boundary)
    x_n, y_n = qdggs.cell(("N",)).nucleus(plane=True)
    x_vo, y_vo = qdggs.cell(("O",)).ul_vertex(plane=True)
    x_vs, y_vs = qdggs.cell(("S",)).ul_vertex(plane=True)
    x_s, y_s = qdggs.cell(("S",)).nucleus()
    meridian_west_geom = shapely.LineString([(x_n, y_n), (x_vo, y_vo), (x_vs, y_vs), (x_s, y_s)])

    x_vr, y_vr = qdggs.cell(("R",)).ul_vertex(plane=True)
    x_vr = x_vr + qdggs.cell(("R",)).width(plane=True)
    meridian_east_geom = shapely.LineString([(x_vr, y_vr), (x_vr, y_vr - qdggs.cell(("R",)).width(plane=True))])




    if cell_geom.intersects(meridian_west_geom) or cell_geom.intersects(meridian_east_geom):
        return True
    else:
        return False


for cell_data in cells_data:
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
        # if cell.intersects_meridian(lam=-180.0) or cell.intersects_meridian(lam=-180.0):
        #     continue
        cell_boundary_points = cell.boundary(n=20, plane=False, interior=True)
        cell_geom = shapely.Polygon(cell_boundary_points)
    elif plot == "qpix":
        cell = qdggs.cell(cell_suid)
        # if check_antimeridian(cell):
        #     continue
        cell_boundary_points = cell.boundary(n=20, plane=False, interior=True)
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
cells_gdf.to_file("rhp_res0.gpkg")

grid_gdf.plot(ax=ax, edgecolor="black", linewidth=2.5, zorder=3)


cells_gdf.plot(
    ax=ax,
    linewidth=0.4,
    zorder=1,
    legend=True,
    edgecolor="black",

)

# plt.savefig("RHP_compactness_diff_nside3_res5.png",  dpi=900)
plt.show()

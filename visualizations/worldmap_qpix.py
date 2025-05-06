# Module for generating visualizations for rHEALPix DGGS

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import geopandas as gpd
import pyproj
import pandas as pd
import shapely
import cartopy
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from rhealpixdggs.projection_wrapper import Projection
from itertools import pairwise
import numpy as np
from rhealpixdggs.utils import auth_lat
from rhealpixdggs.qpix_dggs import QPixDGGS
from matplotlib import font_manager
import math

fig = plt.figure(figsize=(20, 10), frameon=False)
ax = fig.add_subplot(1, 1, 1)
# ax.set_global()

north_square = 0
south_square = 0
lon_0 = 0
lat_0 = 0  # TODO
ellipsoid = WGS84_ELLIPSOID
ellipsoid.lon_0 = lon_0
ellipsoid.lat_0 = lat_0

rosca_plonka_projection = Projection(
    ellipsoid=ellipsoid,
    proj="rosca_plonka",
    north_square=north_square,
    south_square=south_square,
)
phi_0_auth = np.arcsin(2.0 / 3) * 180 / np.pi  # phi_0 on sphere
phi_0 = auth_lat(phi_0_auth, e=WGS84_ELLIPSOID.e, inverse=True)  # phi_0 on ellipsoid

# coastline from cartopy (Natural Earth)
# coastline_geometry = list(cartopy.feature.NaturalEarthFeature(category="physical", name="coastline", scale="110m").geometries())

# coastline from cartopy exported and generalized/simplified in QGIS
coastline_geometry = gpd.read_file("coastline_simplified.geojson").geometry
coastline_gdf = gpd.GeoDataFrame({"geometry": coastline_geometry})



def region_from_point(p):
    # determine resolution 0 cell from (lam, phi)
    # expects ellipsoidal coordinates in degrees
    lam, phi_ellips = p
    lam = lam * math.pi / 180.0
    phi_ellips = phi_ellips * math.pi / 180.0
    phi = auth_lat(phi_ellips, e=WGS84_ELLIPSOID.e, inverse=False, radians=True)
    lam_r = lam - math.pi / 4
    x_r = math.cos(phi) * math.cos(lam_r)
    y_r = math.cos(phi) * math.sin(lam_r)
    z_r = math.sin(phi)

    if x_r == z_r or -x_r == z_r or y_r == z_r or -y_r == z_r:
        return "on_plane"
    if z_r > x_r and z_r > -x_r and z_r > y_r and z_r > -y_r:
        return "north_polar"
    elif z_r <= x_r and z_r <= -x_r and z_r <= y_r and z_r <= -y_r:
        return "south_polar"
    else:
        return "equatorial"

# projecting coastline with rosca_plonka projection
projected_lines = []
for linestring in coastline_geometry:
    # LineString to lines between each vertex (so that lines crossing phi_0 can be omitted from visualization)
    lines = pairwise(list(zip(list(linestring.xy[0]), list(linestring.xy[1]))))
    point_list = []
    for line in lines:
        point_1 = line[0]
        point_2 = line[1]
        point_1_region = region_from_point(point_1)
        point_2_region = region_from_point(point_2)
        if point_1[0] == 180.0:
            point_1 = (179.99999, point_1[1])
        if point_1[0] == -180.0:
            point_1 = (-179.99999, point_1[1])
        if point_2[0] == 180.0:
            point_2 = (179.99999, point_2[1])
        if point_2[0] == -180.0:
            point_2 = (-179.99999, point_2[1])
        if (point_1_region != point_2_region):
            point_1 = rosca_plonka_projection(point_1[0], point_1[1])
            point_list.append(point_1)
            if len(point_list) > 1:
                projected_lines.append(shapely.LineString(point_list))
            point_list = []
            continue
        # project points
        point_1 = rosca_plonka_projection(point_1[0], point_1[1])
        point_list.append(point_1)

    if len(point_list) > 1:
        point_2 = rosca_plonka_projection(point_2[0], point_2[1])
        point_list.append(point_2)
        projected_lines.append(shapely.LineString(point_list))

coastline_gdf = gpd.GeoDataFrame({"geometry": projected_lines})
coastline_gdf.plot(ax=ax, edgecolor="black", linewidth=2, zorder=2)

grid_gdf = gpd.read_file("grid_10deg.geojson")
grid_geometry = grid_gdf.geometry

# projecting grid with rosca_plonka projection
projected_lines = []
for linestring in grid_geometry:
    # LineString to lines between each vertex (so that lines crossing phi_0 can be omitted from visualization)
    lines = pairwise(list(zip(list(linestring.xy[0]), list(linestring.xy[1]))))
    point_list = []
    for line in lines:
        point_1 = line[0]
        point_2 = line[1]
        point_1_region = region_from_point(point_1)
        point_2_region = region_from_point(point_2)
        if (point_1_region != point_2_region):
            point_1 = rosca_plonka_projection(point_1[0], point_1[1])
            point_list.append(point_1)
            projected_lines.append(shapely.LineString(point_list))
            point_list = []
            continue
        # project points
        point_1 = rosca_plonka_projection(point_1[0], point_1[1])
        point_list.append(point_1)

    if len(point_list) > 1:
        projected_lines.append(shapely.LineString(point_list))

grid_gdf = gpd.GeoDataFrame({"geometry": projected_lines})
grid_gdf.plot(ax=ax, edgecolor=(0.75, 0.75, 0.75), linewidth=1.5, zorder=1)

border_geometry = [
    shapely.segmentize(
        shapely.LineString([(-180 + lon_0, phi_0 + 0.01), (180 + lon_0, phi_0 + 0.01)]),
        max_segment_length=0.05,
    ),
    shapely.segmentize(
        shapely.LineString([(-180 + lon_0, phi_0 - 0.01), (180 + lon_0, phi_0 - 0.01)]),
        max_segment_length=0.05,
    ),
    shapely.segmentize(
        shapely.LineString(
            [(-180 + lon_0, -phi_0 - 0.01), (180 + lon_0, -phi_0 - 0.01)]
        ),
        max_segment_length=0.05,
    ),
    shapely.segmentize(
        shapely.LineString(
            [(-180 + lon_0, -phi_0 + 0.01), (180 + lon_0, -phi_0 + 0.01)]
        ),
        max_segment_length=0.05,
    ),
    shapely.segmentize(
        shapely.LineString([(-180 + lon_0, phi_0), (-180 + lon_0, -phi_0)]),
        max_segment_length=0.05,
    ),
    shapely.segmentize(
        shapely.LineString(
            [(180 - 0.0001 + lon_0, phi_0), (180 - 0.0001 + lon_0, -phi_0)]
        ),
        max_segment_length=0.05,
    ),
]


rdggs = QPixDGGS(ellipsoid=ellipsoid, N_side=3)

border_points = [rdggs.ul_vertex["N"], (rdggs.ul_vertex["P"][0], rdggs.ul_vertex["N"][1]), rdggs.ul_vertex["P"], (rdggs.ul_vertex["R"][0] + (rdggs.ul_vertex["R"][0] - rdggs.ul_vertex["Q"][0]), rdggs.ul_vertex["P"][1]), (rdggs.ul_vertex["R"][0] + (rdggs.ul_vertex["R"][0] - rdggs.ul_vertex["Q"][0]), rdggs.ul_vertex["S"][1]), (rdggs.ul_vertex["P"][0], rdggs.ul_vertex["S"][1]), (rdggs.ul_vertex["P"][0], rdggs.ul_vertex["S"][1] - (rdggs.ul_vertex["O"][1] - rdggs.ul_vertex["S"][1])), (rdggs.ul_vertex["S"][0], rdggs.ul_vertex["S"][1] - (rdggs.ul_vertex["O"][1] - rdggs.ul_vertex["S"][1]))]
border_geom = shapely.LinearRing(border_points)

border_gdf = gpd.GeoDataFrame({"geometry": border_geom, "index":[1]})
border_gdf.plot(ax=ax, edgecolor="black", linewidth=2.5, zorder=3)


# Add font for figure
font_path = "fonts/Oswald-Regular.ttf"
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = prop.get_name()

# plotting DGGS cells

# resolution = 0
# cells = rdggs.grid(resolution=0)
# for cell in cells:
#     cell_boundary = cell.boundary(n=10, plane=True, interior=True)
#     cell_geometry = shapely.LinearRing(cell_boundary + [cell_boundary[0]])
#     cell_gdf = gpd.GeoDataFrame([{"geometry": cell_geometry}])
#     cell_gdf.plot(ax=ax, edgecolor="#0096a0", linewidth=2.5, zorder=4)
#     # backgroundcolor=(1,1,1,0.5)
#     plt.text(cell.nucleus(plane=True)[0], cell.nucleus(plane=True)[1], str(cell), family="sans-serif", ha="center",va="center", fontsize=40, color="#0096a0")

# resolution = 1
cells = rdggs.grid(resolution=1)
for cell in cells:
    cell_boundary = cell.boundary(n=10, plane=True, interior=True)
    cell_geometry = shapely.LinearRing(cell_boundary)
    cell_gdf = gpd.GeoDataFrame([{"geometry": cell_geometry}])
    cell_gdf.plot(ax=ax, edgecolor="#0096a0", linewidth=2.5, zorder=4)
    plt.text(
        cell.nucleus(plane=True)[0],
        cell.nucleus(plane=True)[1],
        str(cell),
        family="sans-serif",
        ha="center",
        va="center",
        fontsize=20,
        color="#0096a0",
        backgroundcolor=(1, 1, 1, 0.5),
    )

plt.savefig("qpix_world_level_1.svg")
plt.show()

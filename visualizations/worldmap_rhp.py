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
from rhealpixdggs.dggs import RHEALPixDGGS
from matplotlib import font_manager

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

rhealpix_projection = Projection(
    ellipsoid=ellipsoid,
    proj="rhealpix",
    north_square=north_square,
    south_square=south_square,
)
phi_0_auth = np.arcsin(2.0 / 3) * 180 / np.pi  # phi_0 on sphere
phi_0 = auth_lat(phi_0_auth, e=WGS84_ELLIPSOID.e, inverse=True)  # phi_0 on ellipsoid

# coastline from cartopy (Naturel Earth)
# coastline_geometry = list(cartopy.feature.NaturalEarthFeature(category="physical", name="coastline", scale="110m").geometries())

# coastline from cartopy exported and generalized/simplified in QGIS
coastline_geometry = gpd.read_file("coastline_simplified.geojson").geometry
coastline_gdf = gpd.GeoDataFrame({"geometry": coastline_geometry})

# projecting coastline with rhealpix projection
projected_lines = []
for linestring in coastline_geometry:
    # LineString to lines between each vertex (so that lines crossing phi_0 can be omitted from visualization)
    lines = pairwise(list(zip(list(linestring.xy[0]), list(linestring.xy[1]))))
    point_list = []
    for line in lines:
        point_1 = line[0]
        point_2 = line[1]
        if point_1[0] == 180.0:
            point_1 = (179.99999, point_1[1])
        if point_1[0] == -180.0:
            point_1 = (-179.99999, point_1[1])
        if point_2[0] == 180.0:
            point_2 = (179.99999, point_2[1])
        if point_2[0] == -180.0:
            point_2 = (-179.99999, point_2[1])
        if (
            (point_1[1] < phi_0 < point_2[1])
            or (point_1[1] > phi_0 > point_2[1])
            or (point_1[1] < -phi_0 < point_2[1])
            or (point_1[1] > -phi_0 > point_2[1])
        ):
            crossing_line = shapely.LineString([point_1, point_2])
            if point_1[1] < phi_0 < point_2[1]:
                phi_0_line = shapely.LineString(
                    [(-180, phi_0 - 0.00001), (180, phi_0 - 0.00001)]
                )
                add_to_lat_next_point = 0.00001 * 2
            if point_1[1] > phi_0 > point_2[1]:
                phi_0_line = shapely.LineString(
                    [(-180, phi_0 + 0.00001), (180, phi_0 + 0.00001)]
                )
                add_to_lat_next_point = -0.00001 * 2
            if point_1[1] < -phi_0 < point_2[1]:
                phi_0_line = shapely.LineString(
                    [(-180, -phi_0 - 0.00001), (180, -phi_0 - 0.00001)]
                )
                add_to_lat_next_point = 0.00001 * 2
            if point_1[1] > -phi_0 > point_2[1]:
                phi_0_line = shapely.LineString(
                    [(-180, -phi_0 + 0.00001), (180, -phi_0 + 0.00001)]
                )
                add_to_lat_next_point = -0.00001 * 2
            intersection_point = crossing_line.intersection(phi_0_line)
            point_1 = rhealpix_projection(point_1[0], point_1[1])
            point_list.append(point_1)
            intersection_point_projected = rhealpix_projection(
                intersection_point.x, intersection_point.y
            )
            point_list.append(intersection_point_projected)
            projected_lines.append(shapely.LineString(point_list))
            next_point_projected = rhealpix_projection(
                intersection_point.x, intersection_point.y + add_to_lat_next_point
            )
            point_list = [next_point_projected]
            continue
        # project points
        point_1 = rhealpix_projection(point_1[0], point_1[1])
        point_list.append(point_1)

    if len(point_list) > 1:
        point_2 = rhealpix_projection(point_2[0], point_2[1])
        point_list.append(point_2)
        projected_lines.append(shapely.LineString(point_list))

coastline_gdf = gpd.GeoDataFrame({"geometry": projected_lines})
coastline_gdf.plot(ax=ax, edgecolor="black", linewidth=2, zorder=2)

grid_gdf = gpd.read_file("grid_10deg.geojson")
grid_geometry = grid_gdf.geometry

# projecting grid with rhealpix projection
projected_lines = []
for linestring in grid_geometry:
    # LineString to lines between each vertex (so that lines crossing phi_0 can be omitted from visualization)
    lines = pairwise(list(zip(list(linestring.xy[0]), list(linestring.xy[1]))))
    point_list = []
    for line in lines:
        point_1 = line[0]
        point_2 = line[1]
        if (
            (point_1[1] < phi_0 < point_2[1])
            or (point_1[1] > phi_0 > point_2[1])
            or (point_1[1] < -phi_0 < point_2[1])
            or (point_1[1] > -phi_0 > point_2[1])
            or abs(point_1[1] == phi_0)
            or abs(point_2[1] == phi_0)
        ):
            point_1 = rhealpix_projection(point_1[0], point_1[1])
            point_list.append(point_1)
            projected_lines.append(shapely.LineString(point_list))
            point_list = []
            continue
        # project points
        point_1 = rhealpix_projection(point_1[0], point_1[1])
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

# projecting border with rhealpix projection
projected_lines = []
for linestring in border_geometry:
    # LineString to lines between each vertex
    lines = list(zip(list(linestring.xy[0]), list(linestring.xy[1])))
    point_list = []
    for point in lines:
        # project point
        point = rhealpix_projection(point[0], point[1])
        point_list.append(point)
    projected_lines.append(shapely.LineString(point_list))

border_gdf = gpd.GeoDataFrame({"geometry": projected_lines})
# border_gdf.plot(ax=ax, edgecolor="black", linewidth=2.5, zorder=3)


# Add font for figure
font_path = "fonts/Oswald-Regular.ttf"
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = prop.get_name()

# plotting DGGS cells
rdggs = RHEALPixDGGS(ellipsoid=ellipsoid, N_side=3)
# resolution = 0
# cells = rdggs.grid(resolution=0)
# for cell in cells:
#     cell_boundary = cell.boundary(n=10, plane=True, interior=True)
#     cell_geometry = shapely.LinearRing(cell_boundary + [cell_boundary[0]])
#     cell_gdf = gpd.GeoDataFrame([{"geometry": cell_geometry}])
#     cell_gdf.plot(ax=ax, edgecolor="#0096a0", linewidth=2.5, zorder=4)
#     # backgroundcolor=(1,1,1,0.5)
#     plt.text(cell.centroid(plane=True)[0], cell.centroid(plane=True)[1], str(cell), family="sans-serif", ha="center",va="center", fontsize=40, color="#0096a0")

# resolution = 1
cells = rdggs.grid(resolution=1)
for cell in cells:
    cell_boundary = cell.boundary(n=10, plane=True, interior=True)
    cell_geometry = shapely.LinearRing(cell_boundary)
    cell_gdf = gpd.GeoDataFrame([{"geometry": cell_geometry}])
    cell_gdf.plot(ax=ax, edgecolor="#0096a0", linewidth=2.5, zorder=4)
    plt.text(
        cell.centroid(plane=True)[0],
        cell.centroid(plane=True)[1],
        str(cell),
        family="sans-serif",
        ha="center",
        va="center",
        fontsize=20,
        color="#0096a0",
        backgroundcolor=(1, 1, 1, 0.5),
    )

plt.savefig("rhealpix_world_level_1.svg")
plt.show()

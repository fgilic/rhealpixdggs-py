import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from shapely import LineString, Polygon
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID, WGS84_ELLIPSOID_RADIANS
from mayavi import mlab
from rhealpixdggs.qpix_dggs import QPixDGGS
from rhealpixdggs.dggs import RHEALPixDGGS
import shapely
import geopandas as gpd
import mayavi

mayavi.mlab.figure(bgcolor=(1, 1, 1), size=(2000, 2000))
mayavi.mlab.move(up=1000)
# EPSG:7030
WGS84_A = WGS84_ELLIPSOID.a
WGS84_F = WGS84_ELLIPSOID.f
WGS84_B = WGS84_ELLIPSOID.b

north_square = 0
south_square = 0
lon_0 = 0
lat_0 = 0  # TODO
ellipsoid = WGS84_ELLIPSOID_RADIANS
ellipsoid.lon_0 = lon_0
ellipsoid.lat_0 = lat_0

# # plot sphere
# la = np.linspace(-np.pi, np.pi, 100)
# fi = np.linspace(-np.pi / 2, np.pi / 2, 100)
# x = np.outer(np.cos(fi), np.cos(la))
# y = np.outer(np.cos(fi), np.sin(la))
# z = np.outer(np.sin(fi), np.ones(np.size(la)))
# s = mlab.mesh(x, y, z, color = (0.5,0.5,0.5))

# plot WGS84 ellipsoid
a = WGS84_A
b = WGS84_A
e_2 = (a ** 2 - b ** 2) / a ** 2
la = np.linspace(-np.pi, np.pi, 100)
fi = np.linspace(-np.pi / 2, np.pi / 2, 100)
x = np.outer(a / np.sqrt(1 - e_2 * (np.sin(fi)) ** 2) * np.cos(fi), np.cos(la))
y = np.outer(a / np.sqrt(1 - e_2 * (np.sin(fi)) ** 2) * np.cos(fi), np.sin(la))
z = np.outer(
    a / np.sqrt(1 - e_2 * (np.sin(fi)) ** 2) * np.sin(fi), np.ones(np.size(la))
)
s = mlab.mesh(x, y, z, color=(1,1,1))  # color=(0.513, 0.596, 0.725)

rdggs = QPixDGGS(ellipsoid=ellipsoid, N_side=2)
cells = rdggs.grid(resolution=0)
for cell in cells:
    cell_boundary = cell.boundary(n=10, plane=False, interior=True)
    points_3d = []
    for point in cell_boundary:
        lam, phi = point
        n = a / math.sqrt(1 - e_2 * (math.sin(phi)) ** 2)
        x = n * math.cos(phi) * math.cos(lam)
        y = n * math.cos(phi) * math.sin(lam)
        z = n * math.sin(phi)
        points_3d.append((x,y,z))
    polygon_3d = Polygon(points_3d)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    mlab.plot3d(x, y, z, color=(0, 0.5882, 0.62745), tube_radius=130000)

# cells = rdggs.grid(resolution=1)
# for cell in cells:
#     cell_boundary = cell.boundary(n=10, plane=False, interior=True)
#     points_3d = []
#     for point in cell_boundary:
#         lam, phi = point
#         n = a / math.sqrt(1 - e_2 * (math.sin(phi)) ** 2)
#         x = n * math.cos(phi) * math.cos(lam)
#         y = n * math.cos(phi) * math.sin(lam)
#         z = n * math.sin(phi)
#         points_3d.append((x,y,z))
#     polygon_3d = Polygon(points_3d)
#     x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
#     mlab.plot3d(x, y, z, color=(0, 0.5882, 0.62745), tube_radius=90000)


# cells = rdggs.grid(resolution=2)
# for cell in cells:
#     cell_boundary = cell.boundary(n=10, plane=False, interior=True)
#     points_3d = []
#     for point in cell_boundary:
#         lam, phi = point
#         n = a / math.sqrt(1 - e_2 * (math.sin(phi)) ** 2)
#         x = n * math.cos(phi) * math.cos(lam)
#         y = n * math.cos(phi) * math.sin(lam)
#         z = n * math.sin(phi)
#         points_3d.append((x,y,z))
#     polygon_3d = Polygon(points_3d)
#     x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
#     mlab.plot3d(x, y, z, color=(0, 0.5882, 0.62745), tube_radius=65000)


# cells = rdggs.grid(resolution=3)
# for cell in cells:
#     cell_boundary = cell.boundary(n=3, plane=False, interior=True)
#     points_3d = []
#     for point in cell_boundary:
#         lam, phi = point
#         n = a / math.sqrt(1 - e_2 * (math.sin(phi)) ** 2)
#         x = n * math.cos(phi) * math.cos(lam)
#         y = n * math.cos(phi) * math.sin(lam)
#         z = n * math.sin(phi)
#         points_3d.append((x,y,z))
#     polygon_3d = Polygon(points_3d)
#     x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
#     mlab.plot3d(x, y, z, color=(0, 0.5882, 0.62745), tube_radius=40000)


# cells = rdggs.grid(resolution=4)
# for cell in cells:
#     cell_boundary = cell.boundary(n=3, plane=False, interior=True)
#     points_3d = []
#     for point in cell_boundary:
#         lam, phi = point
#         n = a / math.sqrt(1 - e_2 * (math.sin(phi)) ** 2)
#         x = n * math.cos(phi) * math.cos(lam)
#         y = n * math.cos(phi) * math.sin(lam)
#         z = n * math.sin(phi)
#         points_3d.append((x,y,z))
#     polygon_3d = Polygon(points_3d)
#     x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
#     mlab.plot3d(x, y, z, color=(0, 0.5882, 0.62745), tube_radius=23000)

coastline_geometry = gpd.read_file("coastline_simplified.geojson").geometry
coastline_gdf = gpd.GeoDataFrame({"geometry": coastline_geometry})

for linestring in coastline_geometry:
    points_3d = []
    for point in linestring.coords:
        lam, phi = point
        lam = lam * math.pi / 180
        phi = phi * math.pi / 180
        n = a / math.sqrt(1 - e_2 * (math.sin(phi)) ** 2)
        x = n * math.cos(phi) * math.cos(lam)
        y = n * math.cos(phi) * math.sin(lam)
        z = n * math.sin(phi)
        points_3d.append((x,y,z))
    if len(points_3d) < 4:
        polygon_3d = LineString(points_3d)
        x, y, z = list(zip(*list(polygon_3d.coords)))
    else:
        polygon_3d = Polygon(points_3d)
        x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    mlab.plot3d(x, y, z, color=(0, 0, 0), tube_radius=30000)

grid_geometry = gpd.read_file("grid_10deg.geojson").geometry
grid_gdf = gpd.GeoDataFrame({"geometry": coastline_geometry})

for linestring in grid_geometry:
    points_3d = []
    for point in linestring.coords:
        lam, phi = point
        lam = lam * math.pi / 180
        phi = phi * math.pi / 180
        n = a / math.sqrt(1 - e_2 * (math.sin(phi)) ** 2)
        x = n * math.cos(phi) * math.cos(lam)
        y = n * math.cos(phi) * math.sin(lam)
        z = n * math.sin(phi)
        points_3d.append((x,y,z))
    if len(points_3d) < 4:
        polygon_3d = LineString(points_3d)
        x, y, z = list(zip(*list(polygon_3d.coords)))
    else:
        polygon_3d = Polygon(points_3d)
        x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    mlab.plot3d(x, y, z, color=(0.6, 0.6, 0.6), tube_radius=15000)


mlab.gcf().scene.parallel_projection = True
mayavi.mlab.savefig("q2_0.png", magnification=2)

mlab.show()

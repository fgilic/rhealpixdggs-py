import pyproj
from rhealpixdggs import ellipsoids
from rhealpixdggs import dggs
from geographiclib.geodesic import Geodesic
from math import floor, log10

# define n_side
n_side = 3

# client = Client(processes=True, n_workers=8)

a_wgs84 = pyproj.list.get_ellps_map()["WGS84"]["a"]
f_wgs84 = 1 / pyproj.list.get_ellps_map()["WGS84"]["rf"]

wgs84_geod = Geodesic.WGS84

wgs84_crs = pyproj.CRS.from_epsg(4326)

wgs84_ellipsoid = ellipsoids.Ellipsoid(a=a_wgs84, f=f_wgs84)
rdggs = dggs.RHEALPixDGGS(
    ellipsoid=wgs84_ellipsoid, N_side=n_side, north_square=0, south_square=0, max_areal_resolution=0.5
)


def sig_figs(x: float, precision: int):
    x = float(x)
    precision = int(precision)

    return round(x, -int(floor(log10(abs(x)))) + (precision - 1))

def cell_perimeter(cell, num_points):
    # get cell boundary and construct shapely polygon
    cell_vertices = cell.boundary(n=num_points, plane=False)

    perimeter = 0
    for point in enumerate(cell_vertices):
        point_1 = point[1]
        try:
            point_2 = cell_vertices[point[0] + 1]
        except IndexError:
            point_2 = cell_vertices[0]
        length = wgs84_geod.Inverse(point_1[1], point_1[0], point_2[1], point_2[0], outmask=1025)
        perimeter += length["s12"]

    # print(perimeter)
    return perimeter

max_resolution = rdggs.max_resolution

# point = (0.000000000001, 40.7)
# point = (45.000000001, 0.0000000001)

point = (45.000000001, 0.0000000001)
for resolution in range(0, max_resolution+1):
    cell = rdggs.cell_from_point(resolution=resolution, p=point, plane=False)
    n_side = rdggs.N_side
    exp = 1
    perimeter_1 = cell_perimeter(cell, n_side ** exp)
    perimeter_2 = cell_perimeter(cell, n_side ** (exp + 1))
    while True:
        if sig_figs(perimeter_1, 10) == sig_figs(perimeter_2, 10):
            print(f"rHEALPix DGGS resolution {resolution}, n_side {n_side}")
            print(f"Cell: {cell.suid}")
            print(f"Point: {point}")
            print(f"num_of_perimeter_densification_points: {n_side}**{exp}")
            print("-----------")
            break
        exp += 1
        perimeter_1 = perimeter_2
        perimeter_2 = cell_perimeter(cell, n_side ** (exp + 1))
import pyproj
import shapely
from pyproj.crs.coordinate_operation import LambertAzimuthalEqualAreaConversion
from rhealpixdggs import ellipsoids
from rhealpixdggs import qpix_dggs
from shapely import Polygon
from shapely.ops import transform
from math import floor, log10

# define n_side
n_side = 3

a_wgs84 = pyproj.list.get_ellps_map()["WGS84"]["a"]
f_wgs84 = 1 / pyproj.list.get_ellps_map()["WGS84"]["rf"]

wgs84_crs = pyproj.CRS.from_epsg(4326)

wgs84_ellipsoid = ellipsoids.Ellipsoid(a=a_wgs84, f=f_wgs84)
rdggs = qpix_dggs.QPixDGGS(
    ellipsoid=wgs84_ellipsoid, N_side=n_side, north_square=0, south_square=0, max_areal_resolution=0.5
)
# function for rounding to a specified number of significant figures
# taken form https://mattgosden.medium.com/rounding-to-significant-figures-in-python-2415661b94c3
def sig_figs(x: float, precision: int):
    x = float(x)
    precision = int(precision)

    return round(x, -int(floor(log10(abs(x)))) + (precision - 1))

def cell_area(cell, num_points):
    cell_nucleus = cell.nucleus(plane=False)

    # get cell boundary in plane and construct shapely polygon
    cell_vertices = cell.boundary(n=num_points, plane=False)
    cell_polygon = Polygon(cell_vertices)

    laea_conversion = LambertAzimuthalEqualAreaConversion(
        latitude_natural_origin=cell_nucleus[1],
        longitude_natural_origin=cell_nucleus[0],
    )
    wgs84_laea = pyproj.crs.ProjectedCRS(
        conversion=laea_conversion, geodetic_crs=wgs84_crs
    )

    transformer = pyproj.Transformer.from_crs(
        crs_from=wgs84_crs, crs_to=wgs84_laea, always_xy=True, allow_ballpark=False
    )

    cell_polygon = shapely.ops.transform(transformer.transform, cell_polygon)
    # print(cell_polygon.area)
    return cell_polygon.area

max_resolution = rdggs.max_resolution

# point = (0.000000000001, 35.385453144)
# point = (45.000000001, 0.0000000001)

point = (0.000000000001, 35.385453144)
for resolution in range(0, max_resolution+1):
    cell = rdggs.cell_from_point(resolution=resolution, p=point, plane=False)
    n_side = rdggs.N_side
    exp = 1
    area_1 = cell_area(cell, n_side ** exp)
    area_2 = cell_area(cell, n_side ** (exp + 1))
    while True:
        if sig_figs(area_1, 10) == sig_figs(area_2, 10):
            print(f"QPIX DGGS resolution {resolution}, n_side {n_side}")
            print(f"Cell: {cell.suid}")
            print(f"Point: {point}")
            print(f"Theoretical area: {cell.area(plane=False)}")
            print(f"Calculated area: {area_1}")
            print(f"num_of_area_densification_points: {n_side}**{exp}")
            print("-----------")
            break
        exp += 1
        area_1 = area_2
        area_2 = cell_area(cell, n_side ** (exp + 1))
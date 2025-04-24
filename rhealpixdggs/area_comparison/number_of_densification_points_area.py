import pyproj
import shapely
from pyproj.crs.coordinate_operation import LambertAzimuthalEqualAreaConversion
from rhealpixdggs import ellipsoids
from rhealpixdggs import dggs
from shapely import Polygon, segmentize
import geopandas as gpd
import pandas as pd
from shapely.ops import transform
# from dask.distributed import Client
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

level = 1
print(f"Number of cells at level {level}: {rdggs.num_cells(res_1=level)}")
print(
    f"Cell area at level {level} (plane): {rdggs.cell_area(resolution=level, plane=True)} m2"
)
print(
    f"Cell area at level {level} (ellipsoid): {rdggs.cell_area(resolution=level, plane=False)} m2"
)

grid = rdggs.grid(level)

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
    print(cell_polygon.area)
    return cell_polygon.area

# def func(cell, n_side):
#     cell_suid = str(cell)
#     cell_region = cell.region()
#     cell_shape = cell.ellipsoidal_shape()
#     cell_area_theoretical = cell.area(plane=False)
#     cell_nucleus = cell.nucleus(plane=False)
#
#     cell_data = {
#         "Cell_suid": cell_suid,
#         "Cell_region": cell_region,
#         "Cell_shape": cell_shape,
#         "Theoretical_area_of_cell": cell_area_theoretical,
#         "Cell_nucleus": str(cell_nucleus),
#     }
#
#     # get cell boundary in plane and construct shapely polygon
#     cell_vertices = cell.boundary(n=num_points, plane=False)
#     cell_polygon = Polygon(cell_vertices)
#
#     laea_conversion = LambertAzimuthalEqualAreaConversion(
#         latitude_natural_origin=cell_nucleus[1],
#         longitude_natural_origin=cell_nucleus[0],
#     )
#     wgs84_laea = pyproj.crs.ProjectedCRS(
#         conversion=laea_conversion, geodetic_crs=wgs84_crs
#     )
#
#     transformer = pyproj.Transformer.from_crs(
#         crs_from=wgs84_crs, crs_to=wgs84_laea, always_xy=True, allow_ballpark=False
#     )
#
#     # project densified polygon from ellipsoid to plane with LAEA
#     cell_polygon = shapely.ops.transform(transformer.transform, cell_polygon)
#
#     cell_data["Calculated_area_of_cell"] = cell_polygon.area
#
#     perimeter = 0
#     for point in enumerate(cell_vertices):
#         point_1 = point[1]
#         try:
#             point_2 = cell_vertices[point[0] + 1]
#         except IndexError:
#             point_2 = cell_vertices[0]
#         length = wgs84_geod.Inverse(point_1[1], point_1[0], point_2[1], point_2[0], outmask=1025)
#         perimeter += length["s12"]
#
#     cell_data["Perimeter"] = perimeter
#     cell_data["Cell_vertices"] = cell.vertices(plane=False)
#
#     return cell_data


max_resolution = rdggs.max_resolution

point = (0, 40.9)
for resolution in range(0, max_resolution+1):
    cell = rdggs.cell_from_point(resolution=resolution, p=point, plane=False)
    n_side = rdggs.N_side
    exp = 1
    area_1 = cell_area(cell, n_side ** exp)
    area_2 = cell_area(cell, n_side ** (exp + 1))
    while True:
        if sig_figs(area_1, 10) == sig_figs(area_2, 10):
            print(f"Cell: {cell.suid}; shape: {cell.ellipsoidal_shape()}")
            print(f"Resolution {resolution}, n_side {n_side}: {n_side}**{exp}")
            print(f"Theoretical area: {cell.area(plane=False)}")
            break
        exp += 1
        area_1 = area_2
        area_2 = cell_area(cell, n_side ** (exp + 1))
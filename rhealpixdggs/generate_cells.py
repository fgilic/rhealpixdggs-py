from rhealpixdggs.dggs import WGS84_003
from shapely import Polygon
import geopandas as gpd
import pyproj
import shapely
from pyproj.crs.coordinate_operation import LambertAzimuthalEqualAreaConversion
import pandas as pd

wgs84_crs = pyproj.CRS.from_epsg(4326)

grid_0 = list(WGS84_003.grid(0))
grid_1 = list(WGS84_003.grid(1))
print(f"WGS84_003 level 1 cell area {WGS84_003.cell_area(resolution=1, plane=False)}")

index = [10000000]
columns = []
for cell in grid_1:
    columns.append(cell.suid)
data = []
for n in index:
    areas = []
    for cell in grid_1:
        points = cell.boundary(n=n, plane=False)
        cell = Polygon(points)
        centroid = (cell.centroid.x, cell.centroid.y)
        if centroid[1] > 90.0:
            centroid = (centroid[0],90.0)
            print(centroid[1])

        laea_conversion = LambertAzimuthalEqualAreaConversion(
            latitude_natural_origin=centroid[1],
            longitude_natural_origin=centroid[0],
        )

        wgs84_laea = pyproj.crs.ProjectedCRS(
            conversion=laea_conversion, geodetic_crs=wgs84_crs
        )
        transformer = pyproj.Transformer.from_crs(
            crs_from=wgs84_crs, crs_to=wgs84_laea, always_xy=True, allow_ballpark=False
        )
        cell_projected = shapely.ops.transform(transformer.transform, cell)
        area = cell_projected.area
        areas.append(area)
    data.append(areas)
    print(n)
pd.DataFrame(data=data, index=index, columns=columns).to_csv("test.csv")


# for cell in geometry:
#     centroid = (cell.centroid.x, cell.centroid.y)
#
#     laea_conversion = LambertAzimuthalEqualAreaConversion(
#         latitude_natural_origin=centroid[1],
#         longitude_natural_origin=centroid[0],
#     )
#
#     wgs84_laea = pyproj.crs.ProjectedCRS(
#         conversion=laea_conversion, geodetic_crs=wgs84_crs
#     )
#     transformer = pyproj.Transformer.from_crs(
#         crs_from=wgs84_crs, crs_to=wgs84_laea, always_xy=True, allow_ballpark=False
#     )
#
#     cell_projected = shapely.ops.transform(transformer.transform, cell)
#     area = cell_projected.area
#     perimeter = cell_projected.length
#     data.append({"area": area, "perimeter": perimeter})
#
# pd.DataFrame(data).to_csv("grid_3_rhp.csv")
#
#
#
# # gpd.GeoDataFrame({"geometry": geometry}).to_file("level_3_qpix_ellips.fgb")
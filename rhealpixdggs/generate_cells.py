from rhealpixdggs.dggs import WGS84_003
from shapely import Polygon
import geopandas as gpd
import pyproj
import shapely
from pyproj.crs.coordinate_operation import LambertAzimuthalEqualAreaConversion
import pandas as pd


grid_0 = WGS84_003.grid(0)
grid_1 = WGS84_003.grid(1)
grid_3 = WGS84_003.grid(3)

geometry = []
for cell in grid_3:
    points = cell.boundary(n=10, plane=False)
    boundary = Polygon(points)
    geometry.append(boundary)
wgs84_crs = pyproj.CRS.from_epsg(4326)

data = []

for cell in geometry:
    centroid = (cell.centroid.x, cell.centroid.y)

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
    perimeter = cell_projected.length
    data.append({"area": area, "perimeter": perimeter})

pd.DataFrame(data).to_csv("grid_3_rhp.csv")



# gpd.GeoDataFrame({"geometry": geometry}).to_file("level_3_qpix_ellips.fgb")
from rhealpixdggs.dggs import RHEALPixDGGS, WGS84_ELLIPSOID
from qpix_dggs import QPixDGGS
from shapely import Polygon
import geopandas as gpd
import pyproj
import shapely
from pyproj.crs.coordinate_operation import LambertAzimuthalEqualAreaConversion
import pandas as pd


dggs_type = "qpix" # qpix or rhp
n_side = 3
resolution = 11
region = "split"
num_of_boundary_dens_points = 3

if region == "japan":
    ul = (120.961, 52.766)
    dr = (157.652, 21.669)
elif region == "nz":
    ul = (163.427, -33.293)
    dr = (184.849, -51.729)
elif region == "s_america":
    ul = (-103.695, 22.897)
    dr = (-21.153, -67.164)
elif region == "north_polar":
    ul = (-180, 90)
    dr = (180, 50)
elif region == "gb":
    ul = (-20.522, 80.932)
    dr = (15.975, 47.551)
elif region == "cyprus":
    ul = (30.5008, 36.9223)
    dr = (36.0844, 28.7402)
elif region == "zagreb":
    ul = (18, 46.6)
    dr = (19, 45.9)
elif region == "split":
    ul = (18.58, 44.02)  # qpix
    dr = (18.62, 43.98)  # qpix
    # ul = (16.424380, 43.513748)  # rhp
    # dr = (16.449130, 43.497384)  # rhp


    pass
wgs84_crs = pyproj.CRS.from_epsg(4326)
if dggs_type == "rhp":
    dggs = RHEALPixDGGS(ellipsoid=WGS84_ELLIPSOID, N_side=n_side)
    cells = dggs.cells_from_region(ul=ul, dr=dr, resolution=resolution, plane=False)

    cells_data = []
    for row in cells:
        for cell in row:
            cell_boundary = cell.boundary(n=num_of_boundary_dens_points, plane=False)
            suid = str(cell)
            cell_geom = Polygon(cell_boundary)
            cells_data.append({"suid": suid, "geometry": cell_geom})
elif dggs_type == "qpix":
    dggs = QPixDGGS(ellipsoid=WGS84_ELLIPSOID, N_side=n_side)
    dggs_cells = RHEALPixDGGS(ellipsoid=WGS84_ELLIPSOID, N_side=n_side)
    cells = dggs_cells.cells_from_region(ul=ul, dr=dr, resolution=resolution, plane=False)

    cells_data = []
    for row in cells:
        for cell in row:
            cell_boundary = dggs.cell(cell.suid).boundary(n=num_of_boundary_dens_points, plane=False)
            suid = str(cell)
            cell_geom = Polygon(cell_boundary)
            cells_data.append({"suid": suid, "geometry": cell_geom})



# gpd.GeoSeries(data=cells_geom, crs="EPSG:4326").to_file(f"cells_{region}_{dggs_type}-{n_side}_res{resolution}.gpkg", driver="GPKG")

gpd.GeoDataFrame(data=cells_data, crs="EPSG:4326", geometry="geometry").to_file(f"cells_{region}_{dggs_type}-{n_side}_res{resolution}.gpkg", driver="GPKG")
#
# index = [10000000]
# columns = []
# for cell in grid_1:
#     columns.append(cell.suid)
# data = []
# for n in index:
#     areas = []
#     for cell in grid_1:
#         points = cell.boundary(n=n, plane=False)
#         cell = Polygon(points)
#         centroid = (cell.centroid.x, cell.centroid.y)
#         if centroid[1] > 90.0:
#             centroid = (centroid[0],90.0)
#             print(centroid[1])
#
#         laea_conversion = LambertAzimuthalEqualAreaConversion(
#             latitude_natural_origin=centroid[1],
#             longitude_natural_origin=centroid[0],
#         )
#
#         wgs84_laea = pyproj.crs.ProjectedCRS(
#             conversion=laea_conversion, geodetic_crs=wgs84_crs
#         )
#         transformer = pyproj.Transformer.from_crs(
#             crs_from=wgs84_crs, crs_to=wgs84_laea, always_xy=True, allow_ballpark=False
#         )
#         cell_projected = shapely.ops.transform(transformer.transform, cell)
#         area = cell_projected.area
#         areas.append(area)
#     data.append(areas)
#     print(n)
# pd.DataFrame(data=data, index=index, columns=columns).to_csv("test.csv")
#
#
# # for cell in geometry:
# #     centroid = (cell.centroid.x, cell.centroid.y)
# #
# #     laea_conversion = LambertAzimuthalEqualAreaConversion(
# #         latitude_natural_origin=centroid[1],
# #         longitude_natural_origin=centroid[0],
# #     )
# #
# #     wgs84_laea = pyproj.crs.ProjectedCRS(
# #         conversion=laea_conversion, geodetic_crs=wgs84_crs
# #     )
# #     transformer = pyproj.Transformer.from_crs(
# #         crs_from=wgs84_crs, crs_to=wgs84_laea, always_xy=True, allow_ballpark=False
# #     )
# #
# #     cell_projected = shapely.ops.transform(transformer.transform, cell)
# #     area = cell_projected.area
# #     perimeter = cell_projected.length
# #     data.append({"area": area, "perimeter": perimeter})
# #
# # pd.DataFrame(data).to_csv("grid_3_rhp.csv")
# #
# #
# #
# # # gpd.GeoDataFrame({"geometry": geometry}).to_file("level_3_qpix_ellips.fgb")
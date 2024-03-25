import rasterio
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from rhealpixdggs.dggs import RHEALPixDGGS
from pyproj import CRS, Transformer
from shapely import Polygon
import fiona
import geopandas as gpd
import pandas as pd
from shapely.geometry import mapping
import time

start_time = time.time()
raster_crs = CRS.from_epsg(32633)
dggs_crs = CRS.from_epsg(4326)
raster_dggs = Transformer.from_crs(raster_crs, dggs_crs)
dggs_raster = Transformer.from_crs(dggs_crs, raster_crs)


def read_RGB_raster(file_path):
    raster_file = rasterio.open(file_path)
    red_band = raster_file.read(1)
    green_band = raster_file.read(2)
    blue_band = raster_file.read(3)
    return (raster_file, red_band, green_band, blue_band)


def get_s2_cells(s2_tci):
    ul_utm = (s2_tci.bounds[0], s2_tci.bounds[3])
    lr_utm = (s2_tci.bounds[2], s2_tci.bounds[1])
    ul_wgs = raster_dggs.transform(*ul_utm)
    lr_wgs = raster_dggs.transform(*lr_utm)

    E = WGS84_ELLIPSOID

    rdggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=2, N_side=3)
    cells = rdggs.cells_from_region(
        12, (ul_wgs[1], ul_wgs[0]), (lr_wgs[1], lr_wgs[0]), plane=False
    )

    return cells


def s2_cell_nn(cell):
    vertices = cell.vertices(plane=False)
    cell_polygon = Polygon([vertices[0], vertices[1], vertices[2], vertices[3]])
    centroid_la, centroid_fi = cell.centroid(plane=False)
    centroid_utm = dggs_raster.transform(centroid_fi, centroid_la)
    row, column = s2_tci.index(*centroid_utm)

    return cell_polygon, row, column


def get_cell_rgb(row, column, red_band, green_band, blue_band):
    try:
        red = red_band[row, column]
        green = green_band[row, column]
        blue = blue_band[row, column]
    except IndexError:
        red = 0
        green = 0
        blue = 0

    return red, green, blue


def write_polygons(cells_polygons, fiona_schema):
    with fiona.open(
        "cells_S2_9.geojson",
        "w",
        driver="GeoJSON",
        crs=CRS.from_epsg(4326),
        schema=fiona_schema,
    ) as cells_geojson:
        for cell_polygon in cells_polygons:
            cell = {
                "geometry": mapping(cell_polygon[0]),
                "properties": dict([("color", cell_polygon[1])]),
            }
            cells_geojson.write(cell)


s2_tci, red_band, green_band, blue_band = read_RGB_raster("https://cloud-storage-fg.s3.eu-central-1.amazonaws.com/dggs/raster.tif")
cells = get_s2_cells(s2_tci)

cell_geometry = []
cell_rgba = []
for cell_group in cells:
    for cell in cell_group:
        cell_polygon, row, column = s2_cell_nn(cell)
        red, green, blue = get_cell_rgb(row, column, red_band, green_band, blue_band)
        cell_geometry.append(cell_polygon)
        cell_rgba.append(f"{red},{green},{blue},255")

#  write_polygons(cells_polygons, fiona_schema)

gdf = gpd.GeoDataFrame(
    data={"rgba": cell_rgba, "geometry": cell_geometry}, crs=CRS.from_epsg(4326)
)

gdf.to_file("cells_s2_12.fgb", engine="pyogrio")




print(f"Elapsed time: {(time.time() - start_time)} s")
"""
p = (0, 45)
c = rdggs.cell_from_point(13, p, plane=False)
print(c.area(plane=False))
"""
pass

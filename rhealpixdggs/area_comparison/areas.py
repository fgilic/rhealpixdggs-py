import pyproj
import shapely
from pyproj.crs.coordinate_operation import LambertAzimuthalEqualAreaConversion
from rhealpixdggs import ellipsoids
from rhealpixdggs import dggs
from shapely import Polygon, segmentize
import geopandas as gpd
import pandas as pd
from shapely.ops import transform
from dask.distributed import Client

if __name__ == "__main__":
    client = Client(processes=True, n_workers=16)

    a_wgs84 = pyproj.list.get_ellps_map()["WGS84"]["a"]
    f_wgs84 = 1 / pyproj.list.get_ellps_map()["WGS84"]["rf"]

    wgs84_crs = pyproj.CRS.from_epsg(4326)

    wgs84_ellipsoid = ellipsoids.Ellipsoid(a=a_wgs84, f=f_wgs84)
    rdggs = dggs.RHEALPixDGGS(
        ellipsoid=wgs84_ellipsoid, N_side=3, north_square=0, south_square=0
    )

    level = 1
    print(f"Number of cells at level {level}: {rdggs.num_cells(res_1=level)}")
    print(f"Cell area at level {level} (plane): {rdggs.cell_area(resolution=level, plane=True)} m2")
    print(
        f"Cell area at level {level} (ellipsoid): {rdggs.cell_area(resolution=level, plane=False)} m2"
    )

    grid = rdggs.grid(level)

    def func(cell):
        cell_suid = str(cell)
        cell_region = cell.region()
        cell_shape = cell.ellipsoidal_shape()
        cell_area_theoretical = cell.area(plane=False)
        cell_nucleus = cell.nucleus(plane=False)

        cell_data = {
            "Cell_suid": cell_suid,
            "Cell_region": cell_region,
            "Cell_shape": cell_shape,
            "Theoretical_area_of_cell": cell_area_theoretical,
            "Cell_nucleus": str(cell_nucleus),
        }

        # get cell vertices in plane and construct shapely polygon
        cell_vertices = cell.vertices(plane=True)
        cell_polygon = Polygon(cell_vertices)

        # densify cell boundary by factor 1200000
        cell_polygon = segmentize(
            cell_polygon, max_segment_length=cell_polygon.length / 1200000
        )

        # project densified cell vertices to ellipsoid and construct shapely polygon on ellipsoid
        cell_vertices_plane = list(zip(cell_polygon.boundary.xy[0], cell_polygon.boundary.xy[1]))
        cell_vertices_ellipsoid = [
            rdggs.rhealpix(*point, inverse=True, region=cell_region)
            for point in cell_vertices_plane
        ]
        cell_polygon = Polygon(cell_vertices_ellipsoid)

        cell_data["geometry"] = cell_polygon

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

        # project densified polygon from ellipsoid to LAEA plane
        cell_polygon = shapely.ops.transform(transformer.transform, cell_polygon)

        cell_data["Calculated_area_of_cell"] = cell_polygon.area

        return cell_data

    # # if input dggs grid is too large to be processed at once, divide it in ten parts and process
    # # each part individualy (creates 11 output files
    # grid = list(grid)
    # n = 0
    # for i in range(0, len(grid) + 1, int(len(grid) / 10)):
    #     futures = client.map(func, grid[i:i + int(len(grid) / 10)])
    #     results = client.gather(futures)
    #
    #     gdf = gpd.GeoDataFrame(data=results, geometry="geometry", crs=wgs84_crs)
    #     gdf.to_file(f"grid-4-44444_ellips-0-5-4_{n}.gpkg", driver="GPKG")
    #     n += 1
    #     client.restart()

    # # merge all output files in one (output in FlatGeobuf and in parquet (without geometry)
    # list_of_files = [f"grid-4-44444_ellips-0-5-4_{x}.gpkg" for x in range(0, 11)]
    #
    # numbe_of_rows = 0
    # dggs = None
    # for file_name in list_of_files:
    #     cells = gpd.read_file(file_name, driver="GPKG")
    #     numbe_of_rows += len(cells.index)
    #     dggs = pd.concat([dggs,cells])
    #     print("pass")
    # print(numbe_of_rows)
    # print(len(dggs.index))
    # pd.DataFrame(dggs[["Cell_suid", "Cell_region", "Cell_shape", "Theoretical_area_of_cell", "Cell_nucleus", "Calculated_area_of_cell"]]).to_parquet("grid-4-44444-0-5-4.parquet")
    # print("done")
    # dggs.to_file("grid-4-44444-ellips-0-5-4_merged.fgb", driver="FlatGeobuf")

    futures = client.map(func, list(grid))
    results = client.gather(futures)

    gdf = gpd.GeoDataFrame(data=results, geometry="geometry", crs=wgs84_crs)

    gdf.to_file("grid_1_0-5-4.fgb", driver="FlatGeobuf")

    pd.DataFrame(gdf[["Cell_suid", "Cell_region", "Cell_shape", "Theoretical_area_of_cell", "Cell_nucleus",
                           "Calculated_area_of_cell"]]).to_parquet("grid_1_0-5-4.parquet")
import pyproj
from rhealpixdggs import ellipsoids
from rhealpixdggs import dggs
from dask.distributed import Client
import shapely
from pyproj.crs.coordinate_operation import LambertAzimuthalEqualAreaConversion
from shapely import Polygon
from geographiclib.geodesic import Geodesic
import pandas as pd
import sys
from math import log

if __name__ == "__main__":
    client = Client(processes=True, n_workers=16)

    n_side = 3
    a_wgs84 = pyproj.list.get_ellps_map()["WGS84"]["a"]
    f_wgs84 = 1 / pyproj.list.get_ellps_map()["WGS84"]["rf"]
    wgs84_crs = pyproj.CRS.from_epsg(4326)
    wgs84_ellipsoid = ellipsoids.Ellipsoid(a=a_wgs84, f=f_wgs84)
    rdggs = dggs.RHEALPixDGGS(
        ellipsoid=wgs84_ellipsoid, north_square=0, south_square=0,
    )
    wgs84_geod = Geodesic.WGS84


    def func(arg):
        cell, num_points = arg
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

        # get cell boundary and construct shapely polygon
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

        cell_data["Calculated_area_of_cell"] = cell_polygon.area

        perimeter = 0
        for point in enumerate(cell_vertices):
            point_1 = point[1]
            try:
                point_2 = cell_vertices[point[0] + 1]
            except IndexError:
                point_2 = cell_vertices[0]
            length = wgs84_geod.Inverse(point_1[1], point_1[0], point_2[1], point_2[0], outmask=1025)
            perimeter += length["s12"]

        cell_data["Calculated_perimeter"] = cell_polygon.area
        cell_data["Number_of_densification_points"] = log(num_points, n_side)

        return cell_data



    for resolution in range(9):
        grid = rdggs.grid(resolution)
        if n_side == 2:
            num_points = max(2, 2 ** (17 - resolution))
        elif n_side == 3:
            num_points = max(3, 3 ** (11 - resolution))
        else:
            sys.exit("n_side not 2 or 3")

        grid_with_num_points = []
        for cell in grid:
            grid_with_num_points.append((cell, num_points))

        futures = client.map(func, grid_with_num_points)
        results = client.gather(futures)

        df = pd.DataFrame(data=results)

        df.to_file(f"area_perimeter_res-{resolution}_nside-{n_side}.fgb")
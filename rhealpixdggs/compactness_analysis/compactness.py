import pyproj
from rhealpixdggs import ellipsoids
from rhealpixdggs import dggs, qpix_dggs
from dask.distributed import Client
import shapely
import math
from pyproj.crs.coordinate_operation import LambertAzimuthalEqualAreaConversion
from shapely import Polygon
from geographiclib.geodesic import Geodesic
import pandas as pd
import sys
from math import log
from shapely.ops import transform

if __name__ == "__main__":
    client = Client(processes=True, n_workers=16)

    n_side = 3
    a_wgs84 = pyproj.list.get_ellps_map()["WGS84"]["a"]
    f_wgs84 = 1 / pyproj.list.get_ellps_map()["WGS84"]["rf"]
    r_wgs84_mean = a_wgs84 * (1 - f_wgs84 / 3)
    wgs84_crs = pyproj.CRS.from_epsg(4326)
    wgs84_ellipsoid = ellipsoids.Ellipsoid(a=a_wgs84, f=f_wgs84)
    rdggs = dggs.RHEALPixDGGS(
        ellipsoid=wgs84_ellipsoid, north_square=0, south_square=0, N_side=n_side
    )
    qdggs = qpix_dggs.QPixDGGS(
        ellipsoid=wgs84_ellipsoid, north_square=0, south_square=0, N_side=n_side
    )
    wgs84_geod = Geodesic.WGS84

    r_auth = wgs84_ellipsoid.R_A


    def func(arg):
        rcell, qcell, num_points = arg
        rcell_suid = str(rcell)
        rcell_region = rcell.region()
        rcell_shape = rcell.ellipsoidal_shape()
        rcell_area_theoretical = rcell.area(plane=False)
        rcell_nucleus = rcell.nucleus(plane=False)

        qcell_suid = str(qcell)
        qcell_area_theoretical = qcell.area(plane=False)
        qcell_nucleus = qcell.nucleus(plane=False)

        cell_data = {
            "rCell_suid": rcell_suid,
            "rCell_region": rcell_region,
            "rCell_shape": rcell_shape,
            "Theoretical_area_of_rcell": rcell_area_theoretical,
            "rCell_nucleus": str(rcell_nucleus),
            "qCell_suid": qcell_suid,
            "Theoretical_area_of_qcell": qcell_area_theoretical,
            "qCell_nucleus": str(qcell_nucleus),
        }

        # get cell boundary and construct shapely polygon
        rcell_vertices = rcell.boundary(n=num_points, plane=False)
        rcell_polygon = Polygon(rcell_vertices)
        qcell_vertices = qcell.boundary(n=num_points, plane=False)
        qcell_polygon = Polygon(qcell_vertices)

        rlaea_conversion = LambertAzimuthalEqualAreaConversion(
            latitude_natural_origin=rcell_nucleus[1],
            longitude_natural_origin=rcell_nucleus[0],
        )
        qlaea_conversion = LambertAzimuthalEqualAreaConversion(
            latitude_natural_origin=qcell_nucleus[1],
            longitude_natural_origin=qcell_nucleus[0],
        )

        rwgs84_laea = pyproj.crs.ProjectedCRS(
            conversion=rlaea_conversion, geodetic_crs=wgs84_crs
        )

        qwgs84_laea = pyproj.crs.ProjectedCRS(
            conversion=qlaea_conversion, geodetic_crs=wgs84_crs
        )

        rtransformer = pyproj.Transformer.from_crs(
            crs_from=wgs84_crs, crs_to=rwgs84_laea, always_xy=True, allow_ballpark=False
        )
        qtransformer = pyproj.Transformer.from_crs(
            crs_from=wgs84_crs, crs_to=qwgs84_laea, always_xy=True, allow_ballpark=False
        )

        rcell_polygon = transform(rtransformer.transform, rcell_polygon)
        qcell_polygon = transform(qtransformer.transform, qcell_polygon)

        rarea = rcell_polygon.area
        qarea = qcell_polygon.area

        cell_data["Calculated_area_of_rcell"] = rcell_polygon.area
        cell_data["Calculated_area_of_qcell"] = qcell_polygon.area

        rperimeter = 0
        for point in enumerate(rcell_vertices):
            point_1 = point[1]
            try:
                point_2 = rcell_vertices[point[0] + 1]
            except IndexError:
                point_2 = rcell_vertices[0]
            length = wgs84_geod.Inverse(point_1[1], point_1[0], point_2[1], point_2[0], outmask=1025)
            rperimeter += length["s12"]

        qperimeter = 0
        for point in enumerate(qcell_vertices):
            point_1 = point[1]
            try:
                point_2 = qcell_vertices[point[0] + 1]
            except IndexError:
                point_2 = qcell_vertices[0]
            length = wgs84_geod.Inverse(point_1[1], point_1[0], point_2[1], point_2[0], outmask=1025)
            qperimeter += length["s12"]

        cell_data["Calculated_perimeter_of_rcell"] = rperimeter
        cell_data["Calculated_perimeter_of_qcell"] = qperimeter
        cell_data["Number_of_densification_points"] = f"{n_side} ** {log(num_points, n_side)}"

        cell_data["rCell_ipq"] = 4 * math.pi * rarea / (rperimeter**2)
        cell_data["qCell_ipq"] = 4 * math.pi * qarea / (qperimeter**2)

        cell_data["rCell_c_sp"] = (4 * math.pi - rarea / r_auth**2) * rarea / rperimeter**2
        cell_data["qCell_c_sp"] = (4 * math.pi - qarea / r_auth**2) * qarea / qperimeter**2

        return cell_data



    for resolution in range(0,9):
        rgrid = list(rdggs.grid(resolution))
        qgrid = list(qdggs.grid(resolution))
        if n_side == 2:
            num_points = max(2, 2 ** (17 - resolution))
        elif n_side == 3:
            num_points = max(3, 3 ** (11 - resolution))
        else:
            sys.exit("n_side not 2 or 3")

        grid_with_num_points = []
        for ind, rcell in enumerate(rgrid):
            grid_with_num_points.append((rcell, qgrid[ind], num_points))

        futures = client.map(func, grid_with_num_points)
        results = client.gather(futures)

        df = pd.DataFrame(data=results)

        df.to_csv(f"compactness_res-{resolution}_nside-{n_side}.csv")
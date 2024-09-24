import pyproj
import shapely
from pyproj.crs.coordinate_operation import LambertAzimuthalEqualAreaConversion
from rhealpixdggs import ellipsoids
from rhealpixdggs import dggs
from shapely import Polygon, segmentize
import geopandas as gpd
import pyogrio
from shapely.ops import transform
from dask.distributed import Client

if __name__ == "__main__":
    client = Client(processes=True)

    a_wgs84 = pyproj.list.get_ellps_map()["WGS84"]["a"]
    f_wgs84 = 1 / pyproj.list.get_ellps_map()["WGS84"]["rf"]

    wgs84_crs = pyproj.CRS.from_epsg(4326)

    wgs84_ellipsoid = ellipsoids.Ellipsoid(a=a_wgs84, f=f_wgs84)
    rdggs = dggs.RHEALPixDGGS(
        ellipsoid=wgs84_ellipsoid, N_side=3, north_square=0, south_square=0
    )

    print(f"Number of cells at level 2: {rdggs.num_cells(res_1=2)}")
    print(f"Cell area at level 2 (plane): {rdggs.cell_area(resolution=2, plane=True)} m2")
    print(
        f"Cell area at level 2 (ellipsoid): {rdggs.cell_area(resolution=2, plane=False)} m2"
    )

    grid2 = rdggs.grid(2)

    cells_data = []


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

        # densify cell boundary by factor 10000
        cell_polygon = segmentize(
            cell_polygon, max_segment_length=cell_polygon.length / 40000
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


    futures = client.map(func, list(grid2))
    results = client.gather(futures)

    gdf = gpd.GeoDataFrame(data=results, geometry="geometry", crs=wgs84_crs)

    gdf.to_file("grid-2-10000_ellips-0-5-4.gpkg", driver="GPKG")
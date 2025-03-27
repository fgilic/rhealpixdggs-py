from rhealpixdggs.dggs import RHEALPixDGGS
from rhealpixdggs.ellipsoids import Ellipsoid
from shapely import Polygon
from pyproj.crs.coordinate_operation import (
    LambertAzimuthalEqualAreaConversion,
    AlbersEqualAreaConversion,
)
import pyproj
import shapely
import geopandas as gpd
import math

n_side = 3  # 2 or 3
projection = "laea"  # "aea" or "laea"
exponent_for_segmentation = 20

# points = [(20,1), (20,20), (20,41), (0, 42), (0, 60), (0, 85), (20,42), (20,60), (20,85)]

cell_suids_for_nside_3 = [
    "Q35303300000030333",
    "Q08036306033006603",
    "Q02006060306633306",
    "N22222242666424224",
    "N24466246246642262",
    "N44222442426464644",
    "N20222250876403105",
    "N23368146038840072",
    "N44202532308475755",
]
cell_suids_for_nside_2 = [
    "Q023332003130203110001132003",
    "Q003332203112023312203312001",
    "Q001110201312203330221132003",
    "N111111111121122221222112111",
    "N112122212111222121212222122",
    "N122211211211112212111221211",
    "N110001110020033221322002010",
    "N103033212011333031312233133",
    "N122201311210002312111320210",
]

wgs84_crs = pyproj.CRS.from_epsg(4326)

WGS84_A = pyproj.get_ellps_map()["WGS84"]["a"]  # 6378137.0
WGS84_F = 1 / pyproj.get_ellps_map()["WGS84"]["rf"]  # 298.257223563

ellipsoid = Ellipsoid(a=WGS84_A, f=WGS84_F)

rdggs_003 = RHEALPixDGGS(
    ellipsoid=ellipsoid,
    N_side=3,
    north_square=0,
    south_square=0,
    max_areal_resolution=0.01,
)
rdggs_002 = RHEALPixDGGS(
    ellipsoid=ellipsoid,
    N_side=2,
    north_square=0,
    south_square=0,
    max_areal_resolution=0.01,
)

max_resolution_3 = rdggs_003.max_resolution
max_resolution_2 = rdggs_002.max_resolution

cell_suid = []
cell_region = []
ellipsoidal_shape = []
nucleus = []
theoretical_area = []
calculated_area = []
error = []
geometry = []
resolution = []

if n_side == 3:
    rdggs = rdggs_003
    max_resolution = max_resolution_3
    cell_suids = cell_suids_for_nside_3
elif n_side == 2:
    rdggs = rdggs_002
    max_resolution = max_resolution_2
    cell_suids = cell_suids_for_nside_2

for suid in cell_suids:
    # cell = rdggs.cell_from_point(resolution=cell_resolution, p=point, plane=False)
    while len(suid) > 0:
        suid_initial = list(suid)
        suid = [suid_initial[0]]
        for item in suid_initial[1:]:
            suid.append(int(item))
        cell = rdggs.cell(tuple(suid))
        cell_suid.append(str(cell))
        cell_region.append(cell.region())

        # there is a bug when N_side is even
        if cell.suid in [("N",), ("S",)]:
            cell_ellipsoidal_shape = "cap"
        else:
            cell_ellipsoidal_shape = cell.ellipsoidal_shape()
        ellipsoidal_shape.append(cell_ellipsoidal_shape)

        cell_nucleus = cell.nucleus(plane=False)
        nucleus.append(cell_nucleus)
        cell_theoretical_area = cell.area(plane=False)
        theoretical_area.append(cell_theoretical_area)
        cell_resolution = cell.resolution
        resolution.append(cell_resolution)

        # it turned out commented part below is not needed
        # if cell_ellipsoidal_shape == "cap":
        #     # Errors in cell area calculation because of approximation of boundary by intermediate points is largest
        #     # for cap cells, therefore, area for cap cells is calculated by analytical formula from the bounding
        #     # phi of the cap cell.
        #     # Formula from Rapp 1991 (p. 42) https://kb.osu.edu/items/5cfd1585-a92a-5933-9469-366a5d6b6894
        #     # For quad cells in equatorial region, area could also be calculated in a similar way, but since
        #     # it is not as simple for dart and skew quad cells, we decided to calculate only cap cell area
        #     # from analytical formula.
        #     # Cap cell is bounded by one parallel and therefore accuracy of its calculation directly influences
        #     # accuracy of area calculation.
        #     phi = cell.nw_vertex(plane=False)[1] * math.pi / 180.0
        #     cell_calculated_area = math.pi * cell.ellipsoid.b**2 * (1 / (1 - cell.ellipsoid.e**2) + 1 / (2 * cell.ellipsoid.e) * math.log((1 + cell.ellipsoid.e) / (1 - cell.ellipsoid.e))) - math.pi * cell.ellipsoid.b**2 * (math.sin(phi) / (1 - cell.ellipsoid.e**2 * (math.sin(phi))**2) + 1 / (2*cell.ellipsoid.e) * math.log((1 + cell.ellipsoid.e * math.sin(phi))/(1 - cell.ellipsoid.e * math.sin(phi))))
        #     calculated_area.append(cell_calculated_area)
        # else:
        if cell_resolution == 0:
            n = 2**exponent_for_segmentation
        else:
            n = int(2**exponent_for_segmentation / (n_side * cell_resolution))
        boundary = cell.boundary(n=n, plane=False)

        if projection == "laea":
            laea_conversion = LambertAzimuthalEqualAreaConversion(
                latitude_natural_origin=cell_nucleus[1],
                longitude_natural_origin=cell_nucleus[0],
            )
            wgs84_laea = pyproj.crs.ProjectedCRS(
                conversion=laea_conversion, geodetic_crs=wgs84_crs
            )

            transformer = pyproj.Transformer.from_crs(
                crs_from=wgs84_crs,
                crs_to=wgs84_laea,
                always_xy=True,
                allow_ballpark=False,
            )

        elif projection == "aea":
            if -0.000001 < cell_nucleus[1] < 0.000001:
                if cell_nucleus[1] >= 0:
                    cell_nucleus = (cell_nucleus[0], 0.01)
                else:
                    cell_nucleus = (cell_nucleus[0], -0.01)

            aea_conversion = AlbersEqualAreaConversion(
                latitude_first_parallel=cell_nucleus[1],
                latitude_second_parallel=cell_nucleus[1],
                latitude_false_origin=cell_nucleus[1],
                longitude_false_origin=cell_nucleus[0],
            )

            wgs84_aea = pyproj.crs.ProjectedCRS(
                conversion=aea_conversion, geodetic_crs=wgs84_crs
            )

            transformer = pyproj.Transformer.from_crs(
                crs_from=wgs84_crs,
                crs_to=wgs84_aea,
                always_xy=True,
                allow_ballpark=False,
            )

        cell_calculated_area = shapely.ops.transform(
            transformer.transform, Polygon(boundary)
        ).area
        calculated_area.append(cell_calculated_area)

        area_error = cell_theoretical_area - cell_calculated_area
        error.append(area_error)

        print(str(cell), cell.ellipsoidal_shape(), cell_calculated_area, area_error)

        boundary_lq = Polygon(cell.boundary(n=10, plane=False))
        geometry.append(boundary_lq)
        suid = suid[:-1]

gdf = gpd.GeoDataFrame(
    {
        "geometry": geometry,
        "cell_suid": cell_suid,
        "resolution": resolution,
        "cell_region": cell_region,
        "ellipsoidal_shape": ellipsoidal_shape,
        "nucleus": nucleus,
        "theoretical_area": theoretical_area,
        "calculated_area": calculated_area,
        "error": error,
    }
)
gdf.to_csv(f"area_nside{n_side}_2-{exponent_for_segmentation}_{projection}.csv")
gdf.to_file(f"area_nside{n_side}_2-{exponent_for_segmentation}_{projection}.fgb")

from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from rhealpixdggs.dggs import RHEALPixDGGS
from shapely import Polygon, Point
import geopandas as gpd
import rasterio
from dask.distributed import Client
import rioxarray
from itertools import chain


# cells = rdggs.cells_from_region(resolution=5, ul=(16.45011,43.51516), dr=(16.46,43.51263), plane=False)
# cells_data = []
# for row in cells:
#     for cell in row:
#         geometry = Polygon(cell.vertices(plane=False))
#         cell_data = {"id":str(cell), "geometry":geometry}
#         cells_data.append(cell_data)
#
# gdf = gpd.GeoDataFrame(cells_data, crs="EPSG:4326")
# gdf.to_parquet("test.parquet")

if __name__ == "__main__":
    client = Client(processes=True, n_workers=8)

    def func(strips):
        ul, lr = strips
        cells = rdggs.cells_from_region(5, ul=ul, dr=lr, plane=False)
        centroids = []

        for row in cells:
            for cell in row:
                centroid = {}
                centroid["geometry"] = Point(cell.nucleus(plane=False))
                value_2010 = rds.sel(
                    x=centroid["geometry"].x,
                    y=centroid["geometry"].y,
                    band=11,
                    method="nearest",
                ).item()
                value_2020 = rds.sel(
                    x=centroid["geometry"].x,
                    y=centroid["geometry"].y,
                    band=21,
                    method="nearest",
                ).item()

                if value_2010 == 190 or value_2020 == 190:
                    centroid["suid"] = str(cell)
                    centroid["built_up_2010"] = 0
                    centroid["built_up_2020"] = 0
                    if value_2010 == 190:
                        centroid["built_up_2010"] = 1
                    if value_2020 == 190:
                        centroid["built_up_2020"] = 1
                    centroids.append(centroid)
        return centroids

    ellipsoid = WGS84_ELLIPSOID
    rdggs = RHEALPixDGGS(ellipsoid=ellipsoid, N_side=13, north_square=0, south_square=0)



    file = "https://cloud-storage-fg.s3.eu-central-1.amazonaws.com/medgu-24/GLC_FCS30D_clip.tif"
    rds = rioxarray.open_rasterio(file)

    ul = (rds.coords["x"].min().item(), rds.coords["y"].max().item())
    lr = (rds.coords["x"].max().item(), rds.coords["y"].min().item())

    number_of_strips = 600
    strip_height = (ul[1] - lr[1]) / number_of_strips

    strips = [((ul[0], ul[1]), (lr[0], ul[1] - strip_height))]
    for x in range(1, number_of_strips):
        ul = [ul[0], ul[1] - strip_height]
        lr = [lr[0], ul[1] - strip_height]
        strips.append([ul, lr])

    strip_groups = []
    for x in range(0, 600, 60):
        strip_groups.append(strips[x : x + 60])

    n = 0
    for strip_group in strip_groups:
        futures = client.map(func, strip_group)
        results = client.gather(futures)
        results = list(chain.from_iterable(results))
        if len(results) == 0:
            n += 1
            continue
        gdf = gpd.GeoDataFrame(data=results, geometry="geometry", crs="EPSG:4326")
        gdf.to_parquet(f"adriatic_strips/nuclei_{n}.parquet")
        n += 1

    # TODO merge all files >= cell_nuclei_built-up_raw.parquet

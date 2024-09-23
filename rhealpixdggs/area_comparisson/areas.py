import pyproj
from rhealpixdggs import ellipsoids
from rhealpixdggs import dggs

a_wgs84 = pyproj.list.get_ellps_map()["WGS84"]["a"]
f_wgs84 = 1 / pyproj.list.get_ellps_map()["WGS84"]["rf"]

wgs84_ellipsoid = ellipsoids.Ellipsoid(a=a_wgs84, f=f_wgs84)
rdggs = dggs.RHEALPixDGGS(ellipsoid=wgs84_ellipsoid, N_side=3, north_square=0, south_square=0)

print(f"Number of cells at level 2: {rdggs.num_cells(res_1=2)}")
print(f"Cell area at level 2 (plane): {rdggs.cell_area(resolution=2, plane=True)} m2")
print(f"Cell area at level 2 (ellipsoid): {rdggs.cell_area(resolution=2, plane=False)} m2")

grid2 = rdggs.grid(2)

for cell in grid2:
    print(cell)

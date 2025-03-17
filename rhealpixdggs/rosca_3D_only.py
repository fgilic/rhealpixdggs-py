import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from shapely import LineString, Polygon
from shapely.plotting import plot_line, plot_points
import itertools
from mayavi import mlab

# 'curved' is curved square
# 'square' is rectangular square with straight sides

# EPSG:7030
WGS84_A = 6378137
WGS84_F = 1 / 298.257223563
WGS84_B = WGS84_A - WGS84_A * WGS84_F


def map_to_curved(point_xy):
    x, y = point_xy
    if x == 0 and y == 0:
        return (0, 0)

    if abs(y) <= abs(x):
        x_c = (
            pow(2, 1 / 4)
            * x
            / math.sqrt(math.pi / 6)
            * (math.sqrt(2) * math.cos(y * math.pi / (12 * x)) - 1)
            / math.sqrt(math.sqrt(2) - math.cos(y * math.pi / (12 * x)))
        )
        y_c = (
            pow(2, 1 / 4)
            * x
            / math.sqrt(math.pi / 6)
            * math.sqrt(2)
            * math.sin(y * math.pi / (12 * x))
            / math.sqrt(math.sqrt(2) - math.cos(y * math.pi / (12 * x)))
        )
    else:
        x_c = (
            pow(2, 1 / 4)
            * y
            / math.sqrt(math.pi / 6)
            * math.sqrt(2)
            * math.sin(x * math.pi / (12 * y))
            / math.sqrt(math.sqrt(2) - math.cos(x * math.pi / (12 * y)))
        )
        y_c = (
            pow(2, 1 / 4)
            * y
            / math.sqrt(math.pi / 6)
            * (math.sqrt(2) * math.cos(x * math.pi / (12 * y)) - 1)
            / math.sqrt(math.sqrt(2) - math.cos(x * math.pi / (12 * y)))
        )
    return (x_c, y_c)


def map_to_square(point_xy):
    x, y = point_xy
    if x == 0 and y == 0:
        return (0, 0)

    if abs(y) <= abs(x):
        x_s = (
            math.sqrt(math.pi / 6)
            / math.sqrt(2)
            * math.copysign(1, x)
            * math.pow(2 * x**2 + y**2, 1 / 4)
            * math.sqrt(abs(x) + math.sqrt(2 * x**2 + y**2))
        )
        y_s = (
            math.sqrt(2)
            / math.sqrt(math.pi / 6)
            * math.pow(2 * x**2 + y**2, 1 / 4)
            * math.sqrt(abs(x) + math.sqrt(2 * x**2 + y**2))
            * (
                math.copysign(1, x) * math.atan(y / x)
                - math.atan(y / math.sqrt(2 * x**2 + y**2))
            )
        )
    else:
        x_s = (
            math.sqrt(2)
            / math.sqrt(math.pi / 6)
            * math.pow(2 * y**2 + x**2, 1 / 4)
            * math.sqrt(abs(y) + math.sqrt(2 * y**2 + x**2))
            * (
                math.copysign(1, y) * math.atan(x / y)
                - math.atan(x / math.sqrt(2 * y**2 + x**2))
            )
        )
        y_s = (
            math.sqrt(math.pi / 6)
            / math.sqrt(2)
            * math.copysign(1, y)
            * math.pow(2 * y**2 + x**2, 1 / 4)
            * math.sqrt(abs(y) + math.sqrt(2 * y**2 + x**2))
        )
    return (x_s, y_s)


def densify_square(square, parts_along_sides=10):
    max_segment_length = square.length / (4 * parts_along_sides)
    densified_square = square.segmentize(max_segment_length)
    return densified_square


def square_to_curved(square):
    densified_square = densify_square(square)
    densified_points_on_square = list(zip(*densified_square.boundary.xy))

    points_on_curved = []
    for point in densified_points_on_square:
        points_on_curved.append(map_to_curved(point))

    curved_square = Polygon(points_on_curved)
    return curved_square


def construct_square(ul_coords, side):
    ul_x, ul_y = ul_coords

    ur_coords = (ul_x + side, ul_y)
    lr_coords = (ul_x + side, ul_y - side)
    ll_coords = (ul_x, ul_y - side)

    return Polygon([ul_coords, ur_coords, lr_coords, ll_coords])


def decompose_square_4(square):
    minx, miny, maxx, maxy = square.bounds
    new_square_size = (maxx - minx) / 2
    ul_square = construct_square((minx, maxy), new_square_size)
    ur_square = construct_square((minx + new_square_size, maxy), new_square_size)
    lr_square = construct_square(
        (minx + new_square_size, miny + new_square_size), new_square_size
    )
    ll_square = construct_square((minx, miny + new_square_size), new_square_size)

    return [ul_square, ur_square, lr_square, ll_square]


def build_square_grid(max_level, aperture=4):
    # TODO implement different apertures
    a = math.sqrt(2 * math.pi / 3)
    square_0 = Polygon(
        [(-a / 2, a / 2), (a / 2, a / 2), (a / 2, -a / 2), (-a / 2, -a / 2)]
    )
    grid = {"level_0": [square_0]}
    for level in range(1, max_level + 1):
        squares_previous_level = grid[f"level_{level - 1}"]
        squares_current_level = []
        for square in squares_previous_level:
            squares_current_level.extend(decompose_square_4(square))
        grid[f"level_{level}"] = squares_current_level
    return grid


def lambert_inverse(point_xy, tangent_plane, geographic=False, radians=True):
    x, y = point_xy
    # TODO catch error if wrong tangent plane
    if tangent_plane == (0, 0, 1):
        x_sphere = math.sqrt(1 - (x**2 + y**2) / 4) * x
        y_sphere = math.sqrt(1 - (x**2 + y**2) / 4) * y
        z_sphere = 1 - (x**2 + y**2) / 2

    elif tangent_plane == (1, 0, 0):
        x_sphere = 1 - (x**2 + y**2) / 2
        y_sphere = math.sqrt(1 - (x**2 + y**2) / 4) * x
        z_sphere = math.sqrt(1 - (x**2 + y**2) / 4) * y

    elif tangent_plane == (0, 1, 0):
        x_sphere = math.sqrt(1 - (x**2 + y**2) / 4) * x
        y_sphere = 1 - (x**2 + y**2) / 2
        z_sphere = math.sqrt(1 - (x**2 + y**2) / 4) * y

    elif tangent_plane == (-1, 0, 0):
        x_sphere = -(1 - (x**2 + y**2) / 2)
        y_sphere = math.sqrt(1 - (x**2 + y**2) / 4) * x
        z_sphere = math.sqrt(1 - (x**2 + y**2) / 4) * y

    elif tangent_plane == (0, -1, 0):
        x_sphere = math.sqrt(1 - (x**2 + y**2) / 4) * x
        y_sphere = -(1 - (x**2 + y**2) / 2)
        z_sphere = math.sqrt(1 - (x**2 + y**2) / 4) * y

    elif tangent_plane == (0, 0, -1):
        x_sphere = math.sqrt(1 - (x**2 + y**2) / 4) * x
        y_sphere = math.sqrt(1 - (x**2 + y**2) / 4) * y
        z_sphere = -(1 - (x**2 + y**2) / 2)

    # rotate by 45deg around z axis
    rotation_matrix = np.array(
        [[math.cos(math.pi / 4), -math.sin(math.pi / 4), 0], [math.sin(math.pi / 4), math.cos(math.pi / 4), 0],
         [0, 0, 1]])
    vector = np.array([x_sphere, y_sphere, z_sphere])
    x_sphere, y_sphere, z_sphere = np.matmul(rotation_matrix, vector)

    if not geographic:
        return (x_sphere, y_sphere, z_sphere)
    else:
        long = math.atan2(y_sphere, x_sphere)
        lat = math.asin(z_sphere)

        if not radians:
            lat = lat * 180 / math.pi
            long = long * 180 / math.pi

        return lat, long


def map_sphere_to_ellips(point_lat_long, return_geographic=True, radians=True):
    # authalic latitude and (common) longitude
    auth_lat, long = point_lat_long

    if not radians:
        auth_lat = math.radians(auth_lat)
        long = math.radians(long)

    a = 6378137
    b = 6356752.3141
    n = (a - b) / (a + b)

    c_auth_to_phi = [
        [4 / 3, 4 / 45, -16 / 35, -2582 / 14175, 60136 / 467775, 28112932 / 212837625],
        [0, 46 / 45, 152 / 945, -11966 / 14175, -21016 / 51975, 251310128 / 638512875],
        [0, 0, 3044 / 2835, 3802 / 14175, -94388 / 66825, -8797648 / 10945935],
        [0, 0, 0, 6059 / 4725, 41072 / 93555, -1472637812 / 638512875],
        [0, 0, 0, 0, 768272 / 467775, 455935736 / 638512875],
        [0, 0, 0, 0, 0, 4210684958 / 1915538625],
    ]

    mat_s = [
        math.sin(2 * auth_lat),
        math.sin(4 * auth_lat),
        math.sin(6 * auth_lat),
        math.sin(8 * auth_lat),
        math.sin(10 * auth_lat),
        math.sin(12 * auth_lat),
    ]

    mat_p = np.transpose([n, n**2, n**3, n**4, n**5, n**6])

    common_lat = auth_lat + np.matmul(np.matmul(mat_s, c_auth_to_phi), mat_p)

    if return_geographic == True:
        return (common_lat, long)
    else:
        e_2 = (a**2 - b**2) / a**2
        n = a / math.sqrt(1 - e_2 * (math.sin(common_lat)) ** 2)
        x = n * math.cos(common_lat) * math.cos(long)
        y = n * math.cos(common_lat) * math.sin(long)
        z = n * math.sin(common_lat)
        return x, y, z


square_grid_0 = build_square_grid(1, aperture=9)["level_0"]
square_grid_1 = build_square_grid(1, aperture=9)["level_1"]
square_grid_2 = build_square_grid(2, aperture=9)["level_2"]
square_grid_3 = build_square_grid(3, aperture=9)["level_3"]
square_grid_4 = build_square_grid(4, aperture=9)["level_4"]

# # plot sphere
# la = np.linspace(-np.pi, np.pi, 100)
# fi = np.linspace(-np.pi / 2, np.pi / 2, 100)
# x = np.outer(np.cos(fi), np.cos(la))
# y = np.outer(np.cos(fi), np.sin(la))
# z = np.outer(np.sin(fi), np.ones(np.size(la)))
# s = mlab.mesh(x, y, z, color = (0.5,0.5,0.5))

# plot WGS84 ellipsoid
a = 6378137
b = 6356752.3141
e_2 = (a**2 - b**2) / a**2
la = np.linspace(-np.pi, np.pi, 100)
fi = np.linspace(-np.pi / 2, np.pi / 2, 100)
x = np.outer(a / np.sqrt(1 - e_2 * (np.sin(fi)) ** 2) * np.cos(fi), np.cos(la))
y = np.outer(a / np.sqrt(1 - e_2 * (np.sin(fi)) ** 2) * np.cos(fi), np.sin(la))
z = np.outer(
    a / np.sqrt(1 - e_2 * (np.sin(fi)) ** 2) * np.sin(fi), np.ones(np.size(la))
)
s = mlab.mesh(x, y, z, color=(0.5, 0.5, 0.5))


for square in square_grid_2:
    x_coords_s, y_coords_s = densify_square(square).boundary.xy
    points_s = list(zip(x_coords_s, y_coords_s))

    curved_square = square_to_curved(square)
    x_coords_c, y_coords_c = curved_square.boundary.xy
    points_c = list(zip(x_coords_c, y_coords_c))

    # # on sphere
    # points_3d_s = []
    # for point in points_s:
    #     points_3d_s.append(lambert_inverse(point, tangent_plane=(0, 0, 1)))
    # polygon_3d = Polygon(points_3d_s)
    # x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    # mlab.plot3d(x, y, z, color=(0.1, 0.1, 0.8), tube_radius=0.015)

    points_3d_c = []
    for point in points_c:
        # points_3d_c.append(lambert_inverse(point, tangent_plane=(0, 0, 1)))
        points_3d_c.append(
            map_sphere_to_ellips(
                lambert_inverse(point, tangent_plane=(0, 0, 1), geographic=True),
                return_geographic=False,
                radians=True,
            )
        )
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    # mlab.plot3d(x, y, z, color=(0.8, 0.1, 0.1), tube_radius=0.015)
    mlab.plot3d(x, y, z, color=(0.8, 0.1, 0.1), tube_radius=120000)

    points_3d_c = []
    for point in points_c:
        # points_3d_c.append(lambert_inverse(point, tangent_plane=(0, 0, -1)))
        points_3d_c.append(
            map_sphere_to_ellips(
                lambert_inverse(point, tangent_plane=(0, 0, -1), geographic=True),
                return_geographic=False,
                radians=True,
            )
        )
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    mlab.plot3d(x, y, z, color=(0.8, 0.1, 0.1), tube_radius=120000)

    points_3d_c = []
    for point in points_c:
        # points_3d_c.append(lambert_inverse(point, tangent_plane=(1, 0, 0)))
        points_3d_c.append(
            map_sphere_to_ellips(
                lambert_inverse(point, tangent_plane=(1, 0, 0), geographic=True),
                return_geographic=False,
                radians=True,
            )
        )
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    mlab.plot3d(x, y, z, color=(0.8, 0.1, 0.1), tube_radius=120000)

    points_3d_c = []
    for point in points_c:
        # points_3d_c.append(lambert_inverse(point, tangent_plane=(-1, 0, 0)))
        points_3d_c.append(
            map_sphere_to_ellips(
                lambert_inverse(point, tangent_plane=(-1, 0, 0), geographic=True),
                return_geographic=False,
                radians=True,
            )
        )
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    mlab.plot3d(x, y, z, color=(0.8, 0.1, 0.1), tube_radius=120000)

    points_3d_c = []
    for point in points_c:
        # points_3d_c.append(lambert_inverse(point, tangent_plane=(0, -1, 0)))
        points_3d_c.append(
            map_sphere_to_ellips(
                lambert_inverse(point, tangent_plane=(0, -1, 0), geographic=True),
                return_geographic=False,
                radians=True,
            )
        )
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    mlab.plot3d(x, y, z, color=(0.8, 0.1, 0.1), tube_radius=120000)

    points_3d_c = []
    for point in points_c:
        # points_3d_c.append(lambert_inverse(point, tangent_plane=(0, 1, 0)))
        points_3d_c.append(
            map_sphere_to_ellips(
                lambert_inverse(point, tangent_plane=(0, 1, 0), geographic=True),
                return_geographic=False,
                radians=True,
            )
        )
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    mlab.plot3d(x, y, z, color=(0.8, 0.1, 0.1), tube_radius=120000)


for square in square_grid_3:
    x_coords_s, y_coords_s = densify_square(square).boundary.xy
    points_s = list(zip(x_coords_s, y_coords_s))

    curved_square = square_to_curved(square)
    x_coords_c, y_coords_c = curved_square.boundary.xy
    points_c = list(zip(x_coords_c, y_coords_c))

    # # on sphere
    # points_3d_s = []
    # for point in points_s:
    #     points_3d_s.append(lambert_inverse(point, tangent_plane=(0, 0, 1)))
    # polygon_3d = Polygon(points_3d_s)
    # x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    # mlab.plot3d(x, y, z, color=(0.1, 0.1, 0.8), tube_radius=0.015)

    points_3d_c = []
    for point in points_c:
        # points_3d_c.append(lambert_inverse(point, tangent_plane=(0, 0, 1)))
        points_3d_c.append(
            map_sphere_to_ellips(
                lambert_inverse(point, tangent_plane=(0, 0, 1), geographic=True),
                return_geographic=False,
                radians=True,
            )
        )
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    # mlab.plot3d(x, y, z, color=(0.8, 0.1, 0.1), tube_radius=0.015)
    mlab.plot3d(x, y, z, color=(0.1, 0.8, 0.1), tube_radius=100000)

    points_3d_c = []
    for point in points_c:
        # points_3d_c.append(lambert_inverse(point, tangent_plane=(0, 0, -1)))
        points_3d_c.append(
            map_sphere_to_ellips(
                lambert_inverse(point, tangent_plane=(0, 0, -1), geographic=True),
                return_geographic=False,
                radians=True,
            )
        )
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    mlab.plot3d(x, y, z, color=(0.1, 0.8, 0.1), tube_radius=100000)

    points_3d_c = []
    for point in points_c:
        # points_3d_c.append(lambert_inverse(point, tangent_plane=(1, 0, 0)))
        points_3d_c.append(
            map_sphere_to_ellips(
                lambert_inverse(point, tangent_plane=(1, 0, 0), geographic=True),
                return_geographic=False,
                radians=True,
            )
        )
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    mlab.plot3d(x, y, z, color=(0.1, 0.8, 0.1), tube_radius=100000)

    points_3d_c = []
    for point in points_c:
        # points_3d_c.append(lambert_inverse(point, tangent_plane=(-1, 0, 0)))
        points_3d_c.append(
            map_sphere_to_ellips(
                lambert_inverse(point, tangent_plane=(-1, 0, 0), geographic=True),
                return_geographic=False,
                radians=True,
            )
        )
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    mlab.plot3d(x, y, z, color=(0.1, 0.8, 0.1), tube_radius=100000)

    points_3d_c = []
    for point in points_c:
        # points_3d_c.append(lambert_inverse(point, tangent_plane=(0, -1, 0)))
        points_3d_c.append(
            map_sphere_to_ellips(
                lambert_inverse(point, tangent_plane=(0, -1, 0), geographic=True),
                return_geographic=False,
                radians=True,
            )
        )
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    mlab.plot3d(x, y, z, color=(0.1, 0.8, 0.1), tube_radius=100000)

    points_3d_c = []
    for point in points_c:
        # points_3d_c.append(lambert_inverse(point, tangent_plane=(0, 1, 0)))
        points_3d_c.append(
            map_sphere_to_ellips(
                lambert_inverse(point, tangent_plane=(0, 1, 0), geographic=True),
                return_geographic=False,
                radians=True,
            )
        )
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    mlab.plot3d(x, y, z, color=(0.1, 0.8, 0.1), tube_radius=100000)








mlab.show()

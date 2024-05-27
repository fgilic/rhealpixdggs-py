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


# t = np.linspace(-1/math.sqrt(3), 1/math.sqrt(3), 100)
# x = np.sqrt(2 - 2 * t ** 2) / np.sqrt(2 + np.sqrt(2 - 2 * t ** 2))
# y = 2 * t / np.sqrt(2 + np.sqrt(2 - 2 * t ** 2))
# z = np.ones(100)
#
# fig = plt.figure()
# ax_2d = fig.add_subplot(121)
# ax_2d.plot(x, y, color = 'blue')
# ax_2d.set_aspect('equal')
# plt.show()
#
# x_curved = np.append(np.append(np.append(x, np.flip(y)), np.flip(-x)), y)
# y_curved = np.append(np.append(np.append(np.flip(y), -x), np.flip(-y)), x)
# plt.style.use('bmh')

def map_to_curved(point_xy):
    x, y = point_xy
    if x == 0 and y == 0:
        return (0, 0)

    if abs(y) <= abs(x):
        x_c = pow(2, 1/4) * x / math.sqrt(math.pi / 6) * (math.sqrt(2) * math.cos(y * math.pi / (12 * x)) - 1) / math.sqrt(math.sqrt(2) - math.cos(y * math.pi / (12 * x)))
        y_c = pow(2, 1/4) * x / math.sqrt(math.pi / 6) * math.sqrt(2) * math.sin(y * math.pi / (12 * x)) / math.sqrt(math.sqrt(2) - math.cos(y * math.pi / (12 * x)))
    else:
        x_c = pow(2, 1/4) * y / math.sqrt(math.pi / 6) * math.sqrt(2) * math.sin(x * math.pi / (12 * y)) / math.sqrt(math.sqrt(2) - math.cos(x * math.pi / (12 * y)))
        y_c = pow(2, 1/4) * y / math.sqrt(math.pi / 6) * (math.sqrt(2) * math.cos(x * math.pi / (12 * y)) - 1) / math.sqrt(math.sqrt(2) - math.cos(x * math.pi / (12 * y)))
    return (x_c, y_c)

def map_to_square(point_xy):
    x, y = point_xy
    if x == 0 and y == 0:
        return (0, 0)

    if abs(y) <= abs(x):
        x_s = math.sqrt(math.pi / 6) / math.sqrt(2) * math.copysign(1, x) * math.pow(2 * x ** 2 + y ** 2, 1/4) * math.sqrt(abs(x) + math.sqrt(2 * x ** 2 + y ** 2))
        y_s = math.sqrt(2) / math.sqrt(math.pi / 6) * math.pow(2 * x ** 2 + y ** 2, 1/4) * math.sqrt(abs(x) + math.sqrt(2 * x ** 2 + y ** 2)) * (math.copysign(1, x) * math.atan(y / x) - math.atan(y / math.sqrt(2 * x ** 2 + y ** 2)))
    else:
        x_s = math.sqrt(2) / math.sqrt(math.pi / 6) * math.pow(2 * y ** 2 + x ** 2, 1/4) * math.sqrt(abs(y) + math.sqrt(2 * y ** 2 + x ** 2)) * (math.copysign(1, y) * math.atan(x / y) - math.atan(x / math.sqrt(2 * y ** 2 + x ** 2)))
        y_s = math.sqrt(math.pi / 6) / math.sqrt(2) * math.copysign(1, y) * math.pow(2 * y ** 2 + x ** 2, 1/4) * math.sqrt(abs(y) + math.sqrt(2 * y ** 2 + x ** 2))
    return (x_s, y_s)


def densify_square(square, parts_along_sides = 10):
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
    lr_square = construct_square((minx + new_square_size, miny + new_square_size), new_square_size)
    ll_square = construct_square((minx, miny + new_square_size), new_square_size)

    return [ul_square, ur_square, lr_square, ll_square]

def build_square_grid(max_level, aperture = 4):
    a = math.sqrt(2 * math.pi / 3)
    square_0 = Polygon([(- a / 2, a / 2), (a / 2, a / 2), (a / 2, - a / 2), (- a / 2, - a / 2)])
    grid = {'level_0': [square_0]}
    for level in range(1, max_level + 1):
        squares_previous_level = grid[f'level_{level - 1}']
        squares_current_level = []
        for square in squares_previous_level:
            squares_current_level.extend(decompose_square_4(square))
        grid[f'level_{level}'] = squares_current_level
    return grid

def plot_polygon(polygon, color="blue"):
    ax_2d.plot(*polygon.boundary.xy, color = color)

def lambert_inverse(point_xy, tangent_plane, geographic=False):
    x, y = point_xy

    if tangent_plane == (0, 0, 1):
        x_sphere = math.sqrt(1 - (x ** 2 + y ** 2) / 4) * x
        y_sphere = math.sqrt(1 - (x ** 2 + y ** 2) / 4) * y
        z_sphere = 1 - (x ** 2 + y ** 2) / 2

    if tangent_plane == (1, 0, 0):
        x_sphere = 1 - (x ** 2 + y ** 2) / 2
        y_sphere = math.sqrt(1 - (x ** 2 + y ** 2) / 4) * x
        z_sphere = math.sqrt(1 - (x ** 2 + y ** 2) / 4) * y

    if tangent_plane == (0, 1, 0):
        x_sphere = math.sqrt(1 - (x ** 2 + y ** 2) / 4) * x
        y_sphere = 1 - (x ** 2 + y ** 2) / 2
        z_sphere = math.sqrt(1 - (x ** 2 + y ** 2) / 4) * y

    if tangent_plane == (-1, 0, 0):
        x_sphere = -(1 - (x ** 2 + y ** 2) / 2)
        y_sphere = math.sqrt(1 - (x ** 2 + y ** 2) / 4) * x
        z_sphere = math.sqrt(1 - (x ** 2 + y ** 2) / 4) * y

    if tangent_plane == (0, -1, 0):
        x_sphere = math.sqrt(1 - (x ** 2 + y ** 2) / 4) * x
        y_sphere = -(1 - (x ** 2 + y ** 2) / 2)
        z_sphere = math.sqrt(1 - (x ** 2 + y ** 2) / 4) * y

    if tangent_plane == (0, 0, -1):
        x_sphere = math.sqrt(1 - (x ** 2 + y ** 2) / 4) * x
        y_sphere = math.sqrt(1 - (x ** 2 + y ** 2) / 4) * y
        z_sphere = -(1 - (x ** 2 + y ** 2) / 2)

    if not geographic:
        return (x_sphere, y_sphere, z_sphere)
    else:
        fi = math.atan2(y_sphere/x_sphere) * 180 / math.pi
        la = math.asin(z_sphere) * 180 / math.pi
        return (fi, la)


# a = math.sqrt(2 * math.pi / 3)
# square_0 = Polygon([(- a / 2, a / 2), (a / 2, a / 2), (a / 2, - a / 2), (- a / 2, - a / 2)])

fig = plt.figure()
ax_2d = fig.add_subplot(121)
ax_3d = fig.add_subplot(122, projection='3d')
square_grid_4 = build_square_grid(4, aperture=4)["level_0"]

# plot sphere
la = np.linspace(-np.pi, np.pi, 100)
fi = np.linspace(-np.pi / 2, np.pi / 2, 100)
x = np.outer(np.cos(fi), np.cos(la))
y = np.outer(np.cos(fi), np.sin(la))
z = np.outer(np.sin(fi), np.ones(np.size(la)))
s = mlab.mesh(x, y, z, color = (0.5,0.5,0.5))



for square in square_grid_4:
    ax_2d.plot(*(square).boundary.xy, color="blue")
    curved_square = square_to_curved(square)
    ax_2d.plot(*curved_square.boundary.xy, color="red")

    x_coords_s, y_coords_s = densify_square(square).boundary.xy
    points_s = list(zip(x_coords_s, y_coords_s))

    x_coords_c, y_coords_c = curved_square.boundary.xy
    points_c = list(zip(x_coords_c, y_coords_c))

    points_3d_s = []
    for point in points_s:
        points_3d_s.append(lambert_inverse(point, tangent_plane=(0, 0, 1)))
    polygon_3d = Polygon(points_3d_s)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    ax_3d.plot3D(x, y, z, 'blue')


    points_3d_c = []
    for point in points_c:
        points_3d_c.append(lambert_inverse(point, tangent_plane=(0, 0, 1)))
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    ax_3d.plot3D(x, y, z, 'red')
    mlab.plot3d(x, y, z, color=(0,0.8,0.1))


    # curve_c1_s = LineString([(math.sqrt(2 * math.pi / 3) / 2, math.sqrt(2 * math.pi / 3) / 2), (math.sqrt(2 * math.pi / 3) / 2,- math.sqrt(2 * math.pi / 3) / 2)])
    curve_c1_s = LineString([(math.sqrt(2 * math.pi / 3) / 2, math.sqrt(2 * math.pi / 3) / 2), (math.sqrt(2 * math.pi / 3) / 2, 0)])
    curve_c1_s = curve_c1_s.segmentize(curve_c1_s.length / 10)
    points_curve_c1_s = list(zip(*curve_c1_s.coords.xy))
    points_curve_c1_c = []
    for point in points_curve_c1_s:
        points_curve_c1_c.append(map_to_curved(point))
    curve_c1_c = LineString(points_curve_c1_c)
    ax_2d.plot(*curve_c1_s.coords.xy, color="green")
    ax_2d.plot(*curve_c1_c.coords.xy, color="black")
    x_point = ()
    y_point = ()
    z_point = ()
    for point in points_curve_c1_s:
        x, y, z = lambert_inverse(point, tangent_plane=(-1, 0, 0))
        x_point = x_point + (x,)
        y_point = y_point + (y,)
        z_point = z_point + (z,)
    ax_3d.plot3D(x_point, y_point, z_point, 'green')
    x_point = ()
    y_point = ()
    z_point = ()
    for point in points_curve_c1_c:
        x, y, z = lambert_inverse(point, tangent_plane=(-1, 0, 0))
        x_point = x_point + (x,)
        y_point = y_point + (y,)
        z_point = z_point + (z,)
    ax_3d.plot3D(x_point, y_point, z_point, 'black', linewidth=3)
    ax_2d.set_xlabel("x")
    ax_2d.set_ylabel("y")
    ax_3d.set_xlabel("x")
    ax_3d.set_ylabel("y")
    ax_3d.set_zlabel("z")

    points_3d_c = []
    for point in points_c:
        points_3d_c.append(lambert_inverse(point, tangent_plane=(1, 0, 0)))
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    ax_3d.plot3D(x, y, z, 'red')
    mlab.plot3d(x, y, z, color=(0, 0.8, 0.1))

    points_3d_c = []
    for point in points_c:
        points_3d_c.append(lambert_inverse(point, tangent_plane=(0, 1, 0)))
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    ax_3d.plot3D(x, y, z, 'red')
    mlab.plot3d(x, y, z, color=(0, 0.8, 0.1))

    points_3d_c = []
    for point in points_c:
        points_3d_c.append(lambert_inverse(point, tangent_plane=(-1, 0, 0)))
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    ax_3d.plot3D(x, y, z, 'red')
    mlab.plot3d(x, y, z, color=(0, 0.8, 0.1))

    points_3d_c = []
    for point in points_c:
        points_3d_c.append(lambert_inverse(point, tangent_plane=(0, -1, 0)))
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    ax_3d.plot3D(x, y, z, 'red')
    mlab.plot3d(x, y, z, color=(0, 0.8, 0.1))

    points_3d_c = []
    for point in points_c:
        points_3d_c.append(lambert_inverse(point, tangent_plane=(0, 0, -1)))
    polygon_3d = Polygon(points_3d_c)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    ax_3d.plot3D(x, y, z, 'red')
    mlab.plot3d(x, y, z, color=(0, 0.8, 0.1))

    # points_3d = []
    # for point in points:
    #     points_3d.append(lambert_inverse(point, tangent_plane=(0, -1, 0)))
    # polygon_3d = Polygon(points_3d)
    # x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    # ax_3d.plot3D(x, y, z, 'blue')
    #
    # points_3d = []
    # for point in points:
    #     points_3d.append(lambert_inverse(point, tangent_plane=(0, 0, 1)))
    # polygon_3d = Polygon(points_3d)
    # x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    # ax_3d.plot3D(x, y, z, 'blue')
    #
    # points_3d = []
    # for point in points:
    #     points_3d.append(lambert_inverse(point, tangent_plane=(0, 0, -1)))
    # polygon_3d = Polygon(points_3d)
    # x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    # ax_3d.plot3D(x, y, z, 'blue')

# plot_polygon(square_0, color='black')
# plot_polygon(square_to_curved(square_0), color='blue')

# squares_1 = decompose_square_4(square_0)
# for square in squares_1:
#     plot_polygon(square, color='black')
#     plot_polygon(square_to_curved(square), color='blue')
#
# all_points_square = []
# all_points_curved = []
# for square in squares_1:
#     squares_2 = decompose_square_4(square)
#     for square in squares_2:
#         squares_3 = decompose_square_4(square)
#         for square in squares_3:
#             # plot_polygon(square, color='black')
#             plot_polygon(square_to_curved(square), color='red')
#             all_points_square.append(list(zip(*square.boundary.xy)))
#             all_points_curved.append(list(zip(*square_to_curved(square).boundary.xy)))
#
#
# all_points_square = list(itertools.chain(*all_points_square))
# all_points_curved = list(itertools.chain(*all_points_curved))
#
# for square in squares_1:
#     squares_2 = decompose_square_4(square)
#     for square in squares_2:
#         plot_polygon(square, color='black')
#         plot_polygon(square_to_curved(square), color='green')
#
#
# for square in squares_1:
#     plot_polygon(square, color='black')
#     plot_polygon(square_to_curved(square), color='blue')
#
#
# # line = LineString([[-math.sqrt(math.pi / 6), -math.sqrt(math.pi / 6)], [math.sqrt(math.pi / 6), -math.sqrt(math.pi / 6)]])
# # plot_line(line, ax=ax, add_points=False, color='green', lw=1)
# # line_densified = line.segmentize(0.01)
#
# ax.set_aspect('equal')
# ax = fig.add_subplot(122, projection='3d')
#
# xs = []
# ys = []
# zs = []
# for point in all_points_square:
#     x_s, y_s, z_s = lambert_inverse(point)
#     xs.append(x_s)
#     ys.append(y_s)
#     zs.append(z_s)
#
# xc = []
# yc = []
# zc = []
# for point in all_points_curved:
#     x_c, y_c, z_c = lambert_inverse(point)
#     xc.append(x_c)
#     yc.append(y_c)
#     zc.append(z_c)

# ax.plot3D(np.array(xs), np.array(ys), np.array(zs), 'black')
# ax.plot3D(np.array(xc), np.array(yc), np.array(zc), 'blue')

mlab.show()
ax_2d.set_aspect('equal')
ax_3d.set_aspect('equal')
plt.show()
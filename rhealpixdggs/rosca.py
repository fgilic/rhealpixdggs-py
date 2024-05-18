import math
import numpy as np
import matplotlib.pyplot as plt
from shapely import LineString, Polygon
from shapely.plotting import plot_line, plot_points
import itertools

# t = np.linspace(-1/math.sqrt(3), 1/math.sqrt(3), 100)
# x = np.sqrt(2 - 2 * t ** 2) / np.sqrt(2 + np.sqrt(2 - 2 * t ** 2))
# y = 2 * t / np.sqrt(2 + np.sqrt(2 - 2 * t ** 2))
# z = np.ones(100)
#
# x_curved = np.append(np.append(np.append(x, np.flip(y)), np.flip(-x)), y)
# y_curved = np.append(np.append(np.append(np.flip(y), -x), np.flip(-y)), x)
#
# plt.style.use('bmh')

def mapping_to_curved(point_xy, inverse=False):
    x, y = point_xy
    try:
        # TODO check boundary cases (x, y = 0)
        if not inverse:
            if abs(y) <= abs(x):
                u = pow(2, 1/4) * x / math.sqrt(math.pi / 6) * (math.sqrt(2) * math.cos(y * math.pi / (12 * x)) - 1) / math.sqrt(math.sqrt(2) - math.cos(y * math.pi / (12 * x)))
                v = pow(2, 1/4) * x / math.sqrt(math.pi / 6) * math.sqrt(2) * math.sin(y * math.pi / (12 * x)) / math.sqrt(math.sqrt(2) - math.cos(y * math.pi / (12 * x)))
            else:
                u = pow(2, 1/4) * y / math.sqrt(math.pi / 6) * math.sqrt(2) * math.sin(x * math.pi / (12 * y)) / math.sqrt(math.sqrt(2) - math.cos(x * math.pi / (12 * y)))
                v = pow(2, 1/4) * y / math.sqrt(math.pi / 6) * (math.sqrt(2) * math.cos(x * math.pi / (12 * y)) - 1) / math.sqrt(math.sqrt(2) - math.cos(x * math.pi / (12 * y)))
        else:
            if abs(y) <= abs(x):
                u = math.sqrt(math.pi / 6) / math.sqrt(2) * math.copysign(1, x) * math.pow(2 * x ** 2 + y ** 2, 1/4) * math.sqrt(abs(x) + math.sqrt(2 * x ** 2 + y ** 2))
                v = math.sqrt(2) / math.sqrt(math.pi / 6) * math.pow(2 * x ** 2 + y ** 2, 1/4) * math.sqrt(abs(x) + math.sqrt(2 * x ** 2 + y ** 2)) * (math.copysign(1, x) * math.atan(y / x) - math.atan(y / math.sqrt(2 * x ** 2 + y ** 2)))
            else:
                u = math.sqrt(2) / math.sqrt(math.pi / 6) * math.pow(2 * y ** 2 + x ** 2, 1/4) * math.sqrt(abs(y) + math.sqrt(2 * y ** 2 + x ** 2)) * (math.copysign(1, y) * math.atan(x / y) - math.atan(x / math.sqrt(2 * y ** 2 + x ** 2)))
                v = math.sqrt(math.pi / 6) / math.sqrt(2) * math.copysign(1, y) * math.pow(2 * y ** 2 + x ** 2, 1/4) * math.sqrt(abs(y) + math.sqrt(2 * y ** 2 + x ** 2))
        return (u, v)
    except ZeroDivisionError:
        return (x, y)

def densify_square(square, parts_along_sides = 10):
    max_segment_length = square.length / (4 * parts_along_sides)
    densified_square = square.segmentize(max_segment_length)
    return densified_square

def square_to_curved(square):
    densified_square = densify_square(square)
    densified_points_on_square = list(zip(*densified_square.boundary.xy))

    points_on_curved = []
    for point in densified_points_on_square:
        points_on_curved.append(mapping_to_curved(point))

    curved_square = Polygon(points_on_curved)
    return curved_square

def construct_square(ul_coords, side):
    ur_coords = (ul_coords[0] + side, ul_coords[1])
    lr_coords = (ur_coords[0], ur_coords[1] - side)
    ll_coords = (lr_coords[0] - side, lr_coords[1])

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
    dggs = {'level_0': [square_0]}
    for level in range(1, max_level + 1):
        squares_previous_level = dggs[f'level_{level - 1}']
        squares_current_level = []
        for square in squares_previous_level:
            squares_current_level.extend(decompose_square_4(square))
        dggs[f'level_{level}'] = squares_current_level
    return dggs

def plot_polygon(polygon, color="blue"):
    ax_2d.plot(*polygon.boundary.xy, color = color)

def lambert_inverse(point_xy, tangent_plane = (0, 0, 1), geographic=False):
    x, y = point_xy
    x_l = math.sqrt(1 - (x ** 2 + y ** 2) / 4) * x
    y_l = math.sqrt(1 - (x ** 2 + y ** 2) / 4) * y
    z_l = 1 - (x ** 2 + y ** 2) / 2
    if tangent_plane[2] == 1:
        pass
    elif tangent_plane[2] == -1:
        z_l = - z_l
    elif tangent_plane[1] == 1:
        y_l = z_l
        z_l = x_l
        x_l = y_l
    elif tangent_plane[1] == -1:
        y_l = - z_l
        z_l = x_l
        x_l = y_l
    if not geographic:
        return (x_l, y_l, z_l)
    else:
        fi = math.atan2(y_l/x_l) * 180 / math.pi
        la = math.asin(z_l) * 180 / math.pi
        return (fi, la)


# a = math.sqrt(2 * math.pi / 3)
# square_0 = Polygon([(- a / 2, a / 2), (a / 2, a / 2), (a / 2, - a / 2), (- a / 2, - a / 2)])

fig = plt.figure()
ax_2d = fig.add_subplot(121)
ax_3d = fig.add_subplot(122, projection='3d')
square_grid_4 = build_square_grid(4, aperture = 4)["level_4"]

for square in square_grid_4:
    plot_polygon(square)
    curved_square = square_to_curved(square)
    plot_polygon(curved_square, color = "red")
    x_coords, y_coords = curved_square.boundary.xy
    points = list(zip(x_coords, y_coords))
    points_3d = []
    for point in points:
        points_3d.append(lambert_inverse(point, tangent_plane=(0,1,0)))
    polygon_3d = Polygon(points_3d)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    ax_3d.plot3D(x, y, z, 'blue')

    points_3d = []
    for point in points:
        points_3d.append(lambert_inverse(point, tangent_plane=(0, -1, 0)))
    polygon_3d = Polygon(points_3d)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    ax_3d.plot3D(x, y, z, 'blue')

    points_3d = []
    for point in points:
        points_3d.append(lambert_inverse(point, tangent_plane=(0, 0, 1)))
    polygon_3d = Polygon(points_3d)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    ax_3d.plot3D(x, y, z, 'blue')

    points_3d = []
    for point in points:
        points_3d.append(lambert_inverse(point, tangent_plane=(0, 0, -1)))
    polygon_3d = Polygon(points_3d)
    x, y, z = list(zip(*list(polygon_3d.boundary.coords)))
    ax_3d.plot3D(x, y, z, 'blue')

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
ax_2d.set_aspect('equal')
ax_3d.set_aspect('equal')
plt.show()
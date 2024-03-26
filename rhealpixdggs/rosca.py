import math
import numpy as np
import matplotlib.pyplot as plt
from shapely import LineString
from shapely.plotting import plot_line, plot_points

t = np.linspace(-1/math.sqrt(3), 1/math.sqrt(3), 100)
x = np.sqrt(2 - 2 * t ** 2) / np.sqrt(2 + np.sqrt(2 - 2 * t ** 2))
y = 2 * t / np.sqrt(2 + np.sqrt(2 - 2 * t ** 2))
z = np.ones(100)

x_curved = np.append(np.append(np.append(x, np.flip(y)), np.flip(-x)), y)
y_curved = np.append(np.append(np.append(np.flip(y), -x), np.flip(-y)), x)

a = math.sqrt(2 * math.pi / 3)
x_square = [a/2, -a/2, -a/2, a/2, a/2]
y_square = [a/2, a/2, -a/2, -a/2, a/2]

plt.style.use('bmh')

fig, ax = plt.subplots()
ax.plot(x_square, y_square, color='black', lw=1)
ax.plot(x_curved, y_curved, color='blue', lw=1)



def mapping_sq_to_curved(x, y, inverse=False):
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


line = LineString([[-math.sqrt(math.pi / 6), -math.sqrt(math.pi / 6)], [math.sqrt(math.pi / 6), -math.sqrt(math.pi / 6)]])
# plot_line(line, ax=ax, add_points=False, color='green', lw=1)
line_densified = line.segmentize(0.01)

line_densified_curved_points = []
for point in line_densified.coords:
    line_densified_curved_points.append(mapping_sq_to_curved(point[0], point[1]))

line_densified_curved = LineString(line_densified_curved_points)
ax.plot(line.xy[0], line.xy[1], color='green', lw=1)
ax.plot(line_densified_curved.xy[0], line_densified.xy[1], color='red', lw=1)

line_densified_straight_points = []
for point in line_densified_curved.coords:
    line_densified_straight_points.append(mapping_sq_to_curved(point[0], point[1], inverse=True))

line_densified_straight= LineString(line_densified_straight_points)
ax.plot(line_densified_straight.xy[0], line_densified_straight.xy[1], color='yellow', lw=1)

ax.set_aspect('equal')
plt.show()
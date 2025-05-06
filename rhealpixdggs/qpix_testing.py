import math
import numpy as np
import geopandas as gpd
import shapely
from pyproj import Proj
from itertools import pairwise

def map_to_curved(point_xy):
    # takes tuple, returns tuple
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
    # takes tuple, returns tuple
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

def lambert_inverse(point_xy, tangent_plane, geographic=False, radians=True):
    # takes (x, y) tuple, tangent plane tuple
    # returns (lat, long) or (x, y, z) tuple
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

    # rotate for 45deg around z axis
    angle = math.pi / 4
    rotation_matrix = np.array([[math.cos(angle), -math.sin(angle), 0], [math.sin(angle), math.cos(angle), 0], [0, 0, 1]])
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

        return (long, lat)

def lambert(point, tangent_plane, geographic=False, radians=True):
    if geographic:
        long, lat = point
        if not radians:
            long = long * math.pi / 180
            lat = lat * math.pi / 180

        long = long - math.pi/4

        x_sphere = math.cos(lat) * math.cos(long)
        y_sphere = math.cos(lat) * math.sin(long)
        z_sphere = math.sin(lat)
    else:
        x_sphere, y_sphere, z_sphere = point

    # # rotate for -45deg around z axis
    # angle = -math.pi / 4
    # rotation_matrix = np.array([[math.cos(angle), -math.sin(angle), 0], [math.sin(angle), math.cos(angle), 0], [0, 0, 1]])
    # vector = np.array([x_sphere, y_sphere, z_sphere])
    # x_sphere, y_sphere, z_sphere = np.matmul(rotation_matrix, vector)

    if tangent_plane == (0, 0, 1):
        x = math.sqrt(2 / (1 + z_sphere)) * x_sphere
        y = math.sqrt(2 / (1 + z_sphere)) * y_sphere
    if tangent_plane == (0, 0, -1):
        x = math.sqrt(2 / (1 + -z_sphere)) * x_sphere
        y = math.sqrt(2 / (1 + -z_sphere)) * y_sphere
    if tangent_plane == (1, 0, 0):
        x = math.sqrt(2 / (1 + x_sphere)) * y_sphere
        y = math.sqrt(2 / (1 + x_sphere)) * z_sphere
    if tangent_plane == (-1, 0, 0):
        x = math.sqrt(2 / (1 + -x_sphere)) * y_sphere
        y = math.sqrt(2 / (1 + -x_sphere)) * z_sphere
    if tangent_plane == (0, 1, 0):
        x = math.sqrt(2 / (1 + y_sphere)) * x_sphere
        y = math.sqrt(2 / (1 + y_sphere)) * z_sphere
    if tangent_plane == (0, 1, 0):
        x = math.sqrt(2 / (1 + -y_sphere)) * x_sphere
        y = math.sqrt(2 / (1 + -y_sphere)) * z_sphere

    return (x,y)

def rosca_plonka_sphere(lam, phi, north_square=0, south_square=0):

    # rotate by -pi/4 around z-axes for calculating totated x, y, z for
    # determining which tangent plane to use
    lam_r = lam - math.pi / 4
    x_r = math.cos(phi) * math.cos(lam_r)
    y_r = math.cos(phi) * math.sin(lam_r)
    z_r = math.sin(phi)

    if (0.0 <= lam <= math.pi / 2.0) and z_r <= x_r and z_r >= -x_r:
        laea = Proj("+proj=laea +lon_0=45.0 +R=1.0")
        x, y = laea(longitude=lam, latitude=phi, inverse=False, radians=True)
        x_translation = math.sqrt(2 * math.pi / 3) / 2
        x, y = (map_to_square((x, y)))
        level_0 = "Q"
        return ((x + x_translation, y), level_0)
    elif (math.pi / 2 < lam <= math.pi) and z_r <= y_r and z_r >= -y_r:
        laea = Proj("+proj=laea +lon_0=135.0 +R=1.0")
        x, y = laea(longitude=lam, latitude=phi, inverse=False, radians=True)
        x_translation = 3 * math.sqrt(2 * math.pi / 3) / 2
        x, y = (map_to_square((x, y)))
        level_0 = "R"
        return ((x + x_translation, y), level_0)
    elif (-math.pi / 2 <= lam < 0) and z_r >= y_r and z_r <= -y_r:
        laea = Proj("+proj=laea +lon_0=-45.0 +R=1.0")
        x, y = laea(longitude=lam, latitude=phi, inverse=False, radians=True)
        x_translation = -math.sqrt(2 * math.pi / 3) / 2
        x, y = (map_to_square((x, y)))
        level_0 = "P"
        return ((x + x_translation, y), level_0)
    elif (-math.pi <= lam < -math.pi / 2) and z_r >= x_r and z_r <= -x_r:
        laea = Proj("+proj=laea +lon_0=-135.0 +R=1.0")
        x, y = laea(longitude=lam, latitude=phi, inverse=False, radians=True)
        x_translation = -3 * math.sqrt(2 * math.pi / 3) / 2
        x, y = (map_to_square((x, y)))
        level_0 = "O"
        return ((x + x_translation, y), level_0)
    elif z_r > x_r and z_r > -x_r and z_r > y_r and z_r > -y_r:
        if north_square == 0:
            x_translation = -3 * math.sqrt(2 * math.pi / 3) / 2
            lon_0 = -135.0
        elif north_square == 1:
            x_translation = -math.sqrt(2 * math.pi / 3) / 2
            lon_0 = -45.0
        elif north_square == 2:
            x_translation = math.sqrt(2 * math.pi / 3) / 2
            lon_0 = 45.0
        elif north_square == 3:
            x_translation = 3 * math.sqrt(2 * math.pi / 3) / 2
            lon_0 = 135.0
        laea = Proj(f"+proj=laea +lon_0={lon_0} +lat_0=90.0 +R=1.0")
        x, y = laea(longitude=lam, latitude=phi, inverse=False, radians=True)
        x_translation = x_translation
        y_translation = math.sqrt(2 * math.pi / 3)
        x, y = (map_to_square((x, y)))
        level_0 = "N"
        return ((x + x_translation, y + y_translation), level_0)
    else:
        if south_square == 0:
            x_translation = -3 * math.sqrt(2 * math.pi / 3) / 2
            lon_0 = -135.0
        elif south_square == 1:
            x_translation = -math.sqrt(2 * math.pi / 3) / 2
            lon_0 = -45.0
        elif south_square == 2:
            x_translation = math.sqrt(2 * math.pi / 3) / 2
            lon_0 = 45.0
        elif south_square == 3:
            x_translation = 3 * math.sqrt(2 * math.pi / 3) / 2
            lon_0 = 135.0
        laea = Proj(f"+proj=laea +lon_0={lon_0} +lat_0=-90.0 +R=1.0")
        x, y = laea(longitude=lam, latitude=phi, inverse=False, radians=True)
        x_translation = x_translation
        y_translation = -math.sqrt(2 * math.pi / 3)
        x, y = (map_to_square((x, y)))
        level_0 = "S"
        return ((x + x_translation, y + y_translation), level_0)


def rosca_plonka_sphere_inverse(x, y, north_square=0, south_square=0):
    if (0.0 <= x <= math.sqrt(2 * math.pi / 3)) and -math.sqrt(2 * math.pi / 3) / 2 <= y <= math.sqrt(2 * math.pi / 3) / 2:
        # level 0 cell Q
        x_translation = -math.sqrt(2 * math.pi / 3) / 2
        x = x + x_translation
        x, y = map_to_curved((x, y))
        laea = Proj("+proj=laea +lon_0=45.0 +R=1.0")
        lam, phi = laea(x, y, inverse=True, radians=False)
        return lam, phi
    elif math.sqrt(2 * math.pi / 3) < x <= 2 * (math.sqrt(2 * math.pi / 3)) and -math.sqrt(2 * math.pi / 3) / 2 <= y <= math.sqrt(2 * math.pi / 3) / 2:
        # level 0 cell R
        x_translation = -3 * math.sqrt(2 * math.pi / 3) / 2
        x = x + x_translation
        x, y = map_to_curved((x, y))
        laea = Proj("+proj=laea +lon_0=135.0 +R=1.0")
        lam, phi = laea(x, y, inverse=True, radians=False)
        return lam, phi
    elif -math.sqrt(2 * math.pi / 3) <= x < 0.0 and -math.sqrt(2 * math.pi / 3) / 2 <= y <= math.sqrt(2 * math.pi / 3) / 2:
        # level 0 cell P
        x_translation = math.sqrt(2 * math.pi / 3) / 2
        x = x + x_translation
        x, y = map_to_curved((x, y))
        laea = Proj("+proj=laea +lon_0=-45.0 +R=1.0")
        lam, phi = laea(x, y, inverse=True, radians=False)
        return lam, phi
    elif -2 * math.sqrt(2 * math.pi / 3) <= x < -math.sqrt(2 * math.pi / 3) and -math.sqrt(2 * math.pi / 3) / 2 <= y <= math.sqrt(2 * math.pi / 3) / 2:
        # level 0 cell O
        x_translation = 3 * math.sqrt(2 * math.pi / 3) / 2
        x = x + x_translation
        x, y = map_to_curved((x, y))
        laea = Proj("+proj=laea +lon_0=-135.0 +R=1.0")
        lam, phi = laea(x, y, inverse=True, radians=False)
        return lam, phi
    elif y > math.sqrt(2 * math.pi / 3) / 2:
        # level 0 cell N
        # add more checks in condition
        x_translation = 3 * math.sqrt(2 * math.pi / 3) / 2 - north_square % 4 * math.sqrt(2 * math.pi / 3)
        y_translation = -math.sqrt(2 * math.pi / 3)
        lon_0 = -135.0 + north_square % 4 * 90.0
        x = x + x_translation
        y = y + y_translation
        x, y = map_to_curved((x, y))
        laea = Proj(f"+proj=laea +lon_0={lon_0} +lat_0=90.0 +R=1.0")
        lam, phi = laea(x, y, inverse=True, radians=False)
        return lam, phi
    elif y < -math.sqrt(2 * math.pi / 3) / 2:
        # level 0 cell S
        # add more checks in condition
        x_translation = 3 * math.sqrt(2 * math.pi / 3) / 2 - south_square % 4 * math.sqrt(2 * math.pi / 3)
        y_translation = math.sqrt(2 * math.pi / 3)
        lon_0 = -135.0 + south_square % 4 * 90.0
        x = x + x_translation
        y = y + y_translation
        x, y = map_to_curved((x, y))
        laea = Proj(f"+proj=laea +lon_0={lon_0} +lat_0=-90.0 +R=1.0")
        lam, phi = laea(x, y, inverse=True, radians=False)
        return lam, phi

gdf = gpd.read_file("coastline_proj.fgb")
geometry = gdf.geometry
projected_lines = []
for linestring in geometry:
    projected_points = []
    lines = pairwise(list(zip(list(linestring.xy[0]), list(linestring.xy[1]))))
    for line in lines:
        point_1 = line[0]
        point_2 = line[1]
        projected_point_1 = rosca_plonka_sphere_inverse(point_1[0], point_1[1], north_square=2, south_square=1)
        projected_point_2 = rosca_plonka_sphere_inverse(point_2[0], point_2[1], north_square=2, south_square=1)

        # if (projected_point_1[1] in ["O", "P", "Q", "R"] and projected_point_2[1] in ["N", "S"]) or (projected_point_2[1] in ["O", "P", "Q", "R"] and projected_point_1[1] in ["N", "S"]):
        #     projected_points.append(projected_point_1[0])
        #     if len(projected_points) > 1:
        #         projected_lines.append(shapely.LineString(projected_points))
        #     projected_points = []
        #     continue
        # projected_points.append(projected_point_1[0])
        if projected_point_1 is None:
            continue
        projected_points.append(projected_point_1)



    if len(projected_points) > 1:
        try:
            projected_lines.append(shapely.LineString(projected_points))
        except:
            pass
    # geometry_projected.append(shapely.LineString(points_projected))
x_coastline_projected = gpd.GeoDataFrame({"geometry": projected_lines})
x_coastline_projected.to_file("coastline_proj_inverse.fgb")
x_coastline_projected.plot()
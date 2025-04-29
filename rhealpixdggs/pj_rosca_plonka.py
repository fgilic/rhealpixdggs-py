# I copied the code below from pj_healpix.py and modified it
# DOI: 10.1016/j.cam.2011.07.009
"""
This Python 3.11 module implements the rosca_plonka map projection.

- Alexander Raichev (AR), 2013-01-26: Refactored code from release 0.3.

NOTE:

All lengths are measured in meters and all angles are measured in radians
unless indicated otherwise.
By 'ellipsoid' below, I mean an oblate ellipsoid of revolution.
"""

# *****************************************************************************
#       Copyright (C) 2013 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Import third-party modules.
from pyproj import Proj

# Import standard modules.
import math

# Import my modules.
from rhealpixdggs.utils import auth_rad, auth_lat


def map_to_curved(point_xy):
    # takes (x, y) tuple, returns (x, y) tuple
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
    return x_c, y_c


def map_to_square(point_xy):
    # takes (x, y) tuple, returns (x, y) tuple
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
    return x_s, y_s


def rosca_plonka_sphere(lam, phi, north_square=0, south_square=0):
    # rotate by -pi/4 around z-axes for calculating rotated x, y, z for
    # determining which tangent plane to use
    # (rotating so that prime meridian projection lays between P and Q level 0 cells - as in rhealpix)
    lam_r = lam - math.pi / 4
    x_r = math.cos(phi) * math.cos(lam_r)
    y_r = math.cos(phi) * math.sin(lam_r)
    z_r = math.sin(phi)

    if 0.0 <= lam < math.pi / 2.0 and x_r >= z_r >= -x_r:
        # level 0 cell Q
        laea = Proj("+proj=laea +lon_0=45.0 +R=1.0")
        x, y = laea(longitude=lam, latitude=phi, inverse=False, radians=True)
        x_translation = math.sqrt(2 * math.pi / 3) / 2
        x, y = map_to_square((x, y))
        return x + x_translation, y
    elif math.pi / 2 <= lam <= math.pi and y_r >= z_r >= -y_r:
        # level 0 cell R
        laea = Proj("+proj=laea +lon_0=135.0 +R=1.0")
        x, y = laea(longitude=lam, latitude=phi, inverse=False, radians=True)
        x_translation = 3 * math.sqrt(2 * math.pi / 3) / 2
        x, y = map_to_square((x, y))
        return x + x_translation, y
    elif -math.pi / 2 <= lam < 0.0 and y_r <= z_r <= -y_r:
        # level 0 cell P
        laea = Proj("+proj=laea +lon_0=-45.0 +R=1.0")
        x, y = laea(longitude=lam, latitude=phi, inverse=False, radians=True)
        x_translation = -math.sqrt(2 * math.pi / 3) / 2
        x, y = map_to_square((x, y))
        return x + x_translation, y
    elif -math.pi <= lam < -math.pi / 2 and x_r <= z_r <= -x_r:
        # level 0 cell O
        laea = Proj("+proj=laea +lon_0=-135.0 +R=1.0")
        x, y = laea(longitude=lam, latitude=phi, inverse=False, radians=True)
        x_translation = -3 * math.sqrt(2 * math.pi / 3) / 2
        x, y = map_to_square((x, y))
        return x + x_translation, y
    elif z_r > x_r and z_r > -x_r and z_r > y_r and z_r > -y_r:
        # level 0 cell N
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
        x, y = map_to_square((x, y))
        return x + x_translation, y + y_translation
    else:
        # level 0 cell S
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
        x, y = map_to_square((x, y))
        return x + x_translation, y + y_translation


def rosca_plonka_sphere_inverse(x, y, north_square=0, south_square=0):
    if y > math.sqrt(2 * math.pi / 3) / 2:
        # level 0 cell N
        # TODO add more checks in condition
        if -2 * math.sqrt(2 * math.pi / 3) + (north_square * math.sqrt(2 * math.pi / 3)) <= x < -math.sqrt(2 * math.pi / 3) + (north_square * math.sqrt(2 * math.pi / 3)):
            x_translation = 3 * math.sqrt(
                2 * math.pi / 3
            ) / 2 - north_square * math.sqrt(2 * math.pi / 3)
            y_translation = -math.sqrt(2 * math.pi / 3)
            lon_0 = -135.0 + north_square * 90.0
            x = x + x_translation
            y = y + y_translation
            x, y = map_to_curved((x, y))
            laea = Proj(f"+proj=laea +lon_0={lon_0} +lat_0=90.0 +R=1.0")
            lam, phi = laea(x, y, inverse=True, radians=True)
            return lam, phi
        else:
            y = math.sqrt(2 * math.pi / 3) / 2
        # take care of rounding errors for y
    elif y < -math.sqrt(2 * math.pi / 3) / 2:
        # level 0 cell S
        # TODO add more checks in condition
        if -2 * math.sqrt(2 * math.pi / 3) + (north_square * math.sqrt(2 * math.pi / 3)) <= x < -math.sqrt(2 * math.pi / 3) + (north_square * math.sqrt(2 * math.pi / 3)):
            x_translation = 3 * math.sqrt(
                2 * math.pi / 3
            ) / 2 - south_square * math.sqrt(2 * math.pi / 3)
            y_translation = math.sqrt(2 * math.pi / 3)
            lon_0 = -135.0 + south_square * 90.0
            x = x + x_translation
            y = y + y_translation
            x, y = map_to_curved((x, y))
            laea = Proj(f"+proj=laea +lon_0={lon_0} +lat_0=-90.0 +R=1.0")
            lam, phi = laea(x, y, inverse=True, radians=True)
            return lam, phi
        else:
            # take care of rounding errors for y
            y = -math.sqrt(2 * math.pi / 3) / 2
    if (0.0 <= x < math.sqrt(2 * math.pi / 3)) and -math.sqrt(
        2 * math.pi / 3
    ) / 2 <= y <= math.sqrt(2 * math.pi / 3) / 2:
        # level 0 cell Q
        x_translation = -math.sqrt(2 * math.pi / 3) / 2
        x = x + x_translation
        x, y = map_to_curved((x, y))
        laea = Proj("+proj=laea +lon_0=45.0 +R=1.0")
        lam, phi = laea(x, y, inverse=True, radians=True)
        return lam, phi
    elif (
        math.sqrt(2 * math.pi / 3) <= x <= 2 * (math.sqrt(2 * math.pi / 3))
        and -math.sqrt(2 * math.pi / 3) / 2 <= y <= math.sqrt(2 * math.pi / 3) / 2
    ):
        # level 0 cell R
        x_translation = -3 * math.sqrt(2 * math.pi / 3) / 2
        x = x + x_translation
        x, y = map_to_curved((x, y))
        laea = Proj("+proj=laea +lon_0=135.0 +R=1.0")
        lam, phi = laea(x, y, inverse=True, radians=True)
        return lam, phi
    elif (
        -math.sqrt(2 * math.pi / 3) <= x < 0.0
        and -math.sqrt(2 * math.pi / 3) / 2 <= y <= math.sqrt(2 * math.pi / 3) / 2
    ):
        # level 0 cell P
        x_translation = math.sqrt(2 * math.pi / 3) / 2
        x = x + x_translation
        x, y = map_to_curved((x, y))
        laea = Proj("+proj=laea +lon_0=-45.0 +R=1.0")
        lam, phi = laea(x, y, inverse=True, radians=True)
        return lam, phi
    elif (
        -2 * math.sqrt(2 * math.pi / 3) <= x < -math.sqrt(2 * math.pi / 3)
        and -math.sqrt(2 * math.pi / 3) / 2 <= y <= math.sqrt(2 * math.pi / 3) / 2
    ):
        # level 0 cell O
        x_translation = 3 * math.sqrt(2 * math.pi / 3) / 2
        x = x + x_translation
        x, y = map_to_curved((x, y))
        laea = Proj("+proj=laea +lon_0=-135.0 +R=1.0")
        lam, phi = laea(x, y, inverse=True, radians=True)
        return lam, phi



def rosca_plonka_ellipsoid(lam, phi, e, north_square, south_square):
    beta = auth_lat(phi, e, radians=True, inverse=False)
    return rosca_plonka_sphere(lam, beta, north_square, south_square)


def rosca_plonka_ellipsoid_inverse(x, y, e, north_square, south_square):
    lam, beta = rosca_plonka_sphere_inverse(x, y, north_square, south_square)
    phi = auth_lat(beta, e, radians=True, inverse=True)
    return lam, phi


def rosca_plonka(a, e, north_square, south_square):
    r_a = auth_rad(a, e)

    def f(u, v, radians=False, inverse=False):
        if not inverse:
            lam, phi = u, v
            if not radians:
                # Convert to radians.
                lam = lam * math.pi / 180.0
                phi = phi * math.pi / 180.0
            x, y = rosca_plonka_ellipsoid(
                lam, phi, e=e, north_square=north_square, south_square=south_square
            )
            return r_a * x, r_a * y
        else:
            # Scale down to r_a = 1.
            x, y = u / r_a, v / r_a
            lam, phi = rosca_plonka_ellipsoid_inverse(
                x, y, e=e, north_square=north_square, south_square=south_square
            )
            if not radians:
                # Convert to degrees.
                lam = lam * 180.0 / math.pi
                phi = phi * 180.0 / math.pi
            return lam, phi

    return f

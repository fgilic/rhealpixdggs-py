import timeit

setup_wrap_long_master = """
from numpy import pi, floor, deg2rad, rad2deg
from typing import Any
def wrap_longitude(lam: float, radians: bool = False) -> float:
    if not radians:
        # Convert to radians.
        lam = deg2rad(lam)
    if lam < -pi or lam >= pi:
        result = lam - 2 * pi * floor(lam / (2 * pi))  # x mod 2*pi
        if result >= pi:
            result = result - 2 * pi
    else:
        result = lam
    if not radians:
        # Convert to degrees.
        result = rad2deg(result)
    return result
    """

exec_time = (
    timeit.timeit(
        stmt="wrap_longitude(pi, radians=False)", setup=setup_wrap_long_master, number=1000000
    )
    * 1000
)
print(f"Exec. time (wrap_longitude() master branch, degrees) = {exec_time}")

exec_time = (
    timeit.timeit(
        stmt="wrap_longitude(pi, radians=True)", setup=setup_wrap_long_master, number=1000000
    )
    * 1000
)
print(f"Exec. time (wrap_longitude() master branch, radians) = {exec_time}")

setup_wrap_long_pr = """
from math import pi
from typing import Any
def wrap_longitude(lam: float, radians: bool = False) -> float:
    if radians:
        if lam < -pi or lam >= pi:
            result = lam % (2 * pi)
            if result >= pi:
                result = result - 2 * pi
        else:
            result = lam
        return result

    else:
        if lam < -180 or lam >= 180:
            result = lam % (360)
            if result >= 180:
                result = result - 360
        else:
            result = lam
        return result"""


exec_time = (
    timeit.timeit(
        stmt="wrap_longitude(pi, radians=False)", setup=setup_wrap_long_pr, number=1000000
    )
    * 1000
)
print(f"Exec. time (wrap_longitude() PR branch, degrees) = {exec_time}")

exec_time = (
    timeit.timeit(
        stmt="wrap_longitude(pi, radians=True)", setup=setup_wrap_long_pr, number=1000000
    )
    * 1000
)
print(f"Exec. time (wrap_longitude() PR branch, radians) = {exec_time}")

setup_wrap_long_pr_new = """
from math import pi
from typing import Any
def wrap_longitude(lam: float, radians: bool = False) -> float:
    if radians:
        half_range = pi
    else:
        half_range = 180

    if lam < -half_range or lam >= half_range:
        result = lam % (2 * half_range)
        if result >= half_range:
            result = result - 2 * half_range
    else:
        result = lam
    return result"""


exec_time = (
    timeit.timeit(
        stmt="wrap_longitude(pi, radians=False)", setup=setup_wrap_long_pr_new, number=1000000
    )
    * 1000
)
print(f"Exec. time (wrap_longitude() PR branch new, degrees) = {exec_time}")

exec_time = (
    timeit.timeit(
        stmt="wrap_longitude(pi, radians=True)", setup=setup_wrap_long_pr_new, number=1000000
    )
    * 1000
)
print(f"Exec. time (wrap_longitude() PR branch new, radians) = {exec_time}")
print()


# -------------------------------------------------------

setup_wrap_lat_master = """
from numpy import pi, floor, deg2rad, rad2deg, sign
from typing import Any
def wrap_longitude(lam: float, radians: bool = False) -> float:
    if not radians:
        # Convert to radians.
        lam = deg2rad(lam)
    if lam < -pi or lam >= pi:
        result = lam - 2 * pi * floor(lam / (2 * pi))  # x mod 2*pi
        if result >= pi:
            result = result - 2 * pi
    else:
        result = lam
    if not radians:
        # Convert to degrees.
        result = rad2deg(result)
    return result
def wrap_latitude(phi: float, radians: bool = False) -> float:
    if not radians:
        # Convert to radians.
        phi = deg2rad(phi)
    # Put phi in range -pi <= phi < pi.
    phi = wrap_longitude(phi, radians=True)
    if abs(phi) <= pi / 2:
        result = phi
    else:
        result = phi - sign(phi) * pi
    if not radians:
        # Convert to degrees.
        result = rad2deg(result)
    return result
    """

exec_time = (
    timeit.timeit(
        stmt="wrap_latitude(pi, radians=False)", setup=setup_wrap_lat_master, number=1000000
    )
    * 1000
)
print(f"Exec. time (wrap_latitude() master branch, degrees) = {exec_time}")

exec_time = (
    timeit.timeit(
        stmt="wrap_latitude(pi, radians=False)", setup=setup_wrap_lat_master, number=1000000
    )
    * 1000
)
print(f"Exec. time (wrap_latitude() master branch, radians) = {exec_time}")

setup_wrap_lat_pr = """
from math import pi, copysign
from typing import Any
def wrap_longitude(lam: float, radians: bool = False) -> float:
    if radians:
        if lam < -pi or lam >= pi:
            result = lam % (2 * pi)
            if result >= pi:
                result = result - 2 * pi
        else:
            result = lam
        return result

    else:
        if lam < -180 or lam >= 180:
            result = lam % (360)
            if result >= 180:
                result = result - 360
        else:
            result = lam
        return result

def wrap_latitude(phi: float, radians: bool = False) -> float:
# Put phi in range -pi <= phi < pi.
    phi = wrap_longitude(phi, radians=radians)

    if radians:
        if abs(phi) <= pi / 2:
            result = phi
        else:
            result = phi - copysign(pi, phi)
        return result
    else:
        if abs(phi) <= 180 / 2:
            result = phi
        else:
            result = phi - copysign(180, phi)
        return result        
"""


exec_time = (
    timeit.timeit(
        stmt="wrap_latitude(pi, radians=False)", setup=setup_wrap_lat_pr, number=1000000
    )
    * 1000
)
print(f"Exec. time (wrap_latitude() PR branch, degrees) = {exec_time}")

exec_time = (
    timeit.timeit(
        stmt="wrap_latitude(pi, radians=True)", setup=setup_wrap_lat_pr, number=1000000
    )
    * 1000
)
print(f"Exec. time (wrap_latitude() PR branch, radians) = {exec_time}")

setup_wrap_lat_pr_new = """
from math import pi, copysign
from typing import Any
def wrap_longitude(lam: float, radians: bool = False) -> float:
    if radians:
        half_range = pi
    else:
        half_range = 180

    if lam < -half_range or lam >= half_range:
        result = lam % (2 * half_range)
        if result >= half_range:
            result = result - 2 * half_range
    else:
        result = lam
    return result
def wrap_latitude(phi: float, radians: bool = False) -> float:
# Put phi in range -pi <= phi < pi.
    phi = wrap_longitude(phi, radians=radians)

    if radians:
        half_range = pi
    else:
        half_range = 180

    if abs(phi) <= half_range / 2:
        result = phi
    else:
        result = phi - copysign(half_range, phi)
    return result
    """


exec_time = (
    timeit.timeit(
        stmt="wrap_latitude(pi, radians=False)", setup=setup_wrap_lat_pr_new, number=1000000
    )
    * 1000
)
print(f"Exec. time (wrap_latitude() PR branch new, degrees) = {exec_time}")

exec_time = (
    timeit.timeit(
        stmt="wrap_latitude(pi, radians=True)", setup=setup_wrap_lat_pr_new, number=1000000
    )
    * 1000
)
print(f"Exec. time (wrap_latitude() PR branch new, radians) = {exec_time}")

print()

# -------------------------------------------------------

setup_auth_lat_master = """
from numpy import sign, rad2deg, deg2rad, pi, sin, arcsin, log
def auth_lat(
    phi: float, e: float, inverse: bool = False, radians: bool = False
) -> float:
    if e == 0:
        return phi
    if not radians:
        # Convert to radians to do calculations below.
        phi = deg2rad(phi)
    if not inverse:
        # Compute authalic latitude from latitude phi.
        q = ((1 - e**2) * sin(phi)) / (1 - (e * sin(phi)) ** 2) - (1 - e**2) / (
            2.0 * e
        ) * log((1 - e * sin(phi)) / (1 + e * sin(phi)))
        qp = 1 - (1 - e**2) / (2.0 * e) * log((1.0 - e) / (1.0 + e))
        ratio = q / qp
        # Avoid rounding errors.
        if abs(ratio) > 1:
            # Make abs(ratio) = 1
            ratio = sign(ratio)
        result = arcsin(ratio)
    else:
        # Compute an approximation of latitude from authalic latitude phi.
        result = (
            phi
            + (e**2 / 3.0 + 31 * e**4 / 180.0 + 517 * e**6 / 5040.0)
            * sin(2 * phi)
            + (23 * e**4 / 360.0 + 251 * e**6 / 3780.0) * sin(4 * phi)
            + (761 * e**6 / 45360.0) * sin(6 * phi)
        )
    if not radians:
        # Convert back to degrees.
        result = rad2deg(result)
    return result"""

exec_time = (
        timeit.timeit(
            stmt="auth_lat(phi=pi/4, e=0.08181919104281579, radians=False)", setup=setup_auth_lat_master, number=100000
        )
        * 1000
)
print(f"Exec. time (auth_lat() master, degrees) = {exec_time}")

exec_time = (
        timeit.timeit(
            stmt="auth_lat(phi=pi/4, e=0.08181919104281579, radians=True)", setup=setup_auth_lat_master, number=100000
        )
        * 1000
)
print(f"Exec. time (auth_lat() master, radians) = {exec_time}")



setup_auth_lat_pr = """
from math import sqrt, pi, sin, log, asin
def auth_lat(
    phi: float, e: float, inverse: bool = False, radians: bool = False
) -> float:
    if e == 0:
        return phi
    # Compute flattening f and third flattening n from eccentricity e.
    f = 1 - sqrt(1 - e**2)
    n = (1 - sqrt(1 - e**2)) / (1 + sqrt(1 - e**2))

    if not inverse:
        # Compute authalic latitude from common latitude phi.
        # For large flattenings (f > 1/150) use direct formula,
        # for small flattenings (f <= 1/150) use power series.
        if abs(f) > 1 / 150:
            # Use direct formula for large flattenings.
            if not radians:
                # Convert to radians to do calculations below.
                phi = phi * pi / 180
            # Compute authalic latitude from latitude phi.
            q = ((1 - e**2) * sin(phi)) / (1 - (e * sin(phi)) ** 2) - (1 - e**2) / (
                2.0 * e
            ) * log((1 - e * sin(phi)) / (1 + e * sin(phi)))
            qp = 1 - (1 - e**2) / (2.0 * e) * log((1.0 - e) / (1.0 + e))
            ratio = q / qp
            # Avoid rounding errors.
            if abs(ratio) > 1:
                # Make abs(ratio) = 1
                ratio = copysign(1, ratio)
            result = asin(ratio)
            if not radians:
                result = result * 180 / pi
            return result
        else:
            # Use power series approximation for small flattenings (f <= 1/150).
            # Power series expansion taken from https://doi.org/10.48550/arXiv.2212.05818 (Equation A19)
            if radians:
                authalic_lat = phi + (
                    (
                        -4 / 3 * n
                        - 4 / 45 * n**2
                        + 88 / 315 * n**3
                        + 538 / 4725 * n**4
                        + 20824 / 467775 * n**5
                        - 44732 / 2837835 * n**6
                    )
                    * sin(2 * phi)
                    + (
                        34 / 45 * n**2
                        + 8 / 105 * n**3
                        - 2482 / 14175 * n**4
                        - 37192 / 467775 * n**5
                        - 12467764 / 212837625 * n**6
                    )
                    * sin(4 * phi)
                    + (
                        -1532 / 2835 * n**3
                        - 898 / 14175 * n**4
                        + 54968 / 467775 * n**5
                        + 100320856 / 1915538625 * n**6
                    )
                    * sin(6 * phi)
                    + (
                        6007 / 14175 * n**4
                        + 24496 / 467775 * n**5
                        - 5884124 / 70945875 * n**6
                    )
                    * sin(8 * phi)
                    + (-23356 / 66825 * n**5 - 839792 / 19348875 * n**6) * sin(10 * phi)
                    + (570284222 / 1915538625 * n**6) * sin(12 * phi)
                )
                return authalic_lat
            else:
                authalic_lat = phi * pi / 180 + (
                    (
                        -4 / 3 * n
                        - 4 / 45 * n**2
                        + 88 / 315 * n**3
                        + 538 / 4725 * n**4
                        + 20824 / 467775 * n**5
                        - 44732 / 2837835 * n**6
                    )
                    * sin(2 * phi * pi / 180)
                    + (
                        34 / 45 * n**2
                        + 8 / 105 * n**3
                        - 2482 / 14175 * n**4
                        - 37192 / 467775 * n**5
                        - 12467764 / 212837625 * n**6
                    )
                    * sin(4 * phi * pi / 180)
                    + (
                        -1532 / 2835 * n**3
                        - 898 / 14175 * n**4
                        + 54968 / 467775 * n**5
                        + 100320856 / 1915538625 * n**6
                    )
                    * sin(6 * phi * pi / 180)
                    + (
                        6007 / 14175 * n**4
                        + 24496 / 467775 * n**5
                        - 5884124 / 70945875 * n**6
                    )
                    * sin(8 * phi * pi / 180)
                    + (-23356 / 66825 * n**5 - 839792 / 19348875 * n**6)
                    * sin(10 * phi * pi / 180)
                    + (570284222 / 1915538625 * n**6) * sin(12 * phi * pi / 180)
                )
                return authalic_lat * 180 / pi
    else:
        # Compute common latitude from authalic latitude phi.
        # Power series expansion taken from https://doi.org/10.48550/arXiv.2212.05818 (Equation A20)
        if radians:
            common_lat = phi + (
                (
                    4 / 3 * n
                    + 4 / 45 * n**2
                    - 16 / 35 * n**3
                    - 2582 / 14175 * n**4
                    + 60136 / 467775 * n**5
                    + 28112932 / 212837625 * n**6
                )
                * sin(2 * phi)
                + (
                    46 / 45 * n**2
                    + 152 / 945 * n**3
                    - 11966 / 14175 * n**4
                    - 21016 / 51975 * n**5
                    + 251310128 / 638512875 * n**6
                )
                * sin(4 * phi)
                + (
                    3044 / 2835 * n**3
                    + 3802 / 14175 * n**4
                    - 94388 / 66825 * n**5
                    - 8797648 / 10945935 * n**6
                )
                * sin(6 * phi)
                + (
                    6059 / 4725 * n**4
                    + 41072 / 93555 * n**5
                    - 1472637812 / 638512875 * n**6
                )
                * sin(8 * phi)
                + (768272 / 467775 * n**5 + 455935736 / 638512875 * n**6)
                * sin(10 * phi)
                + (4210684958 / 1915538625 * n**6) * sin(12 * phi)
            )
            return common_lat
        else:
            common_lat = phi * pi / 180 + (
                (
                    4 / 3 * n
                    + 4 / 45 * n**2
                    - 16 / 35 * n**3
                    - 2582 / 14175 * n**4
                    + 60136 / 467775 * n**5
                    + 28112932 / 212837625 * n**6
                )
                * sin(2 * phi * pi / 180)
                + (
                    46 / 45 * n**2
                    + 152 / 945 * n**3
                    - 11966 / 14175 * n**4
                    - 21016 / 51975 * n**5
                    + 251310128 / 638512875 * n**6
                )
                * sin(4 * phi * pi / 180)
                + (
                    3044 / 2835 * n**3
                    + 3802 / 14175 * n**4
                    - 94388 / 66825 * n**5
                    - 8797648 / 10945935 * n**6
                )
                * sin(6 * phi * pi / 180)
                + (
                    6059 / 4725 * n**4
                    + 41072 / 93555 * n**5
                    - 1472637812 / 638512875 * n**6
                )
                * sin(8 * phi * pi / 180)
                + (768272 / 467775 * n**5 + 455935736 / 638512875 * n**6)
                * sin(10 * phi * pi / 180)
                + (4210684958 / 1915538625 * n**6) * sin(12 * phi * pi / 180)
            )
            return common_lat * 180 / pi
    """

exec_time = (
    timeit.timeit(
        stmt="auth_lat(phi=pi/4, e=0.08181919104281579, radians=False)", setup=setup_auth_lat_pr, number=100000
    )
    * 1000
)
print(f"Exec. time (auth_lat() PR branch, degrees) = {exec_time}")

exec_time = (
    timeit.timeit(
        stmt="auth_lat(phi=pi/4, e=0.08181919104281579, radians=True)", setup=setup_auth_lat_pr, number=100000
    )
    * 1000
)
print(f"Exec. time (auth_lat() PR branch, radians) = {exec_time}")


setup_auth_lat_pr_new = """
from math import sqrt, pi, sin, log, asin
def auth_lat(
    phi: float, e: float, inverse: bool = False, radians: bool = False
) -> float:
    if e == 0:
        return phi
    # Compute flattening f and third flattening n from eccentricity e.
    f = 1 - sqrt(1 - e**2)
    n = (1 - sqrt(1 - e**2)) / (1 + sqrt(1 - e**2))

    if not inverse:
        # Compute authalic latitude from common latitude phi.
        # For large flattenings (f > 1/150) use direct formula,
        # for small flattenings (f <= 1/150) use power series.
        if abs(f) > 1 / 150:
            # Use direct formula for large flattenings.
            if not radians:
                # Convert to radians to do calculations below.
                phi = phi * pi / 180
            # Compute authalic latitude from latitude phi.
            q = ((1 - e**2) * sin(phi)) / (1 - (e * sin(phi)) ** 2) - (1 - e**2) / (
                2.0 * e
            ) * log((1 - e * sin(phi)) / (1 + e * sin(phi)))
            qp = 1 - (1 - e**2) / (2.0 * e) * log((1.0 - e) / (1.0 + e))
            ratio = q / qp
            # Avoid rounding errors.
            if abs(ratio) > 1:
                # Make abs(ratio) = 1
                ratio = copysign(1, ratio)
            result = asin(ratio)
            if not radians:
                result = result * 180 / pi
            return result
        else:
            # Use power series approximation for small flattenings (f <= 1/150).
            # Power series expansion taken from https://doi.org/10.48550/arXiv.2212.05818 (Equation A19)
            if not radians:
                phi = phi * pi / 180

            authalic_lat = phi + (
                n
                * (
                    -4 / 3
                    + n
                    * (
                        -4 / 45
                        + n
                        * (
                            88 / 315
                            + n
                            * (
                                538 / 4725
                                + n * (20824 / 467775 + n * (-44732 / 2837835))
                            )
                        )
                    )
                )
                * sin(2 * phi)
                + n
                * (
                    n
                    * (
                        34 / 45
                        + n
                        * (
                            8 / 105
                            + n
                            * (
                                -2482 / 14175
                                + n * (-37192 / 467775 + n * (-12467764 / 212837625))
                            )
                        )
                    )
                )
                * sin(4 * phi)
                + n
                * (
                    n
                    * (
                        n
                        * (
                            -1532 / 2835
                            + n
                            * (
                                -898 / 14175
                                + n * (54968 / 467775 + n * 100320856 / 1915538625)
                            )
                        )
                    )
                )
                * sin(6 * phi)
                + n
                * (
                    n
                    * (
                        n
                        * (
                            n
                            * (
                                6007 / 14175
                                + n * (24496 / 467775 + n * (-5884124 / 70945875))
                            )
                        )
                    )
                )
                * sin(8 * phi)
                + n
                * (n * (n * (n * (n * (-23356 / 66825 + n * (-839792 / 19348875))))))
                * sin(10 * phi)
                + n
                * (n * (n * (n * (n * (n * 570284222 / 1915538625)))))
                * sin(12 * phi)
            )

            if not radians:
                authalic_lat = authalic_lat * 180 / pi

            return authalic_lat
    else:
        # Compute common latitude from authalic latitude phi.
        # Power series expansion taken from https://doi.org/10.48550/arXiv.2212.05818 (Equation A20)
        if not radians:
            phi = phi * pi / 180

        common_lat = phi + (
            n
            * (
                4 / 3
                + n
                * (
                    4 / 45
                    + n
                    * (
                        -16 / 35
                        + n
                        * (
                            -2582 / 14175
                            + n * (60136 / 467775 + n * 28112932 / 212837625)
                        )
                    )
                )
            )
            * sin(2 * phi)
            + n
            * (
                n
                * (
                    46 / 45
                    + n
                    * (
                        152 / 945
                        + n
                        * (
                            -11966 / 14175
                            + n * (-21016 / 51975 + n * 251310128 / 638512875)
                        )
                    )
                )
            )
            * sin(4 * phi)
            + n
            * (
                n
                * (
                    n
                    * (
                        3044 / 2835
                        + n
                        * (
                            3802 / 14175
                            + n * (-94388 / 66825 + n * (-8797648 / 10945935))
                        )
                    )
                )
            )
            * sin(6 * phi)
            + n
            * (
                n
                * (
                    n
                    * (
                        n
                        * (
                            6059 / 4725
                            + n * (41072 / 93555 + n * (-1472637812 / 638512875))
                        )
                    )
                )
            )
            * sin(8 * phi)
            + n
            * (n * (n * (n * (n * (768272 / 467775 + n * 455935736 / 638512875)))))
            * sin(10 * phi)
            + n * (n * (n * (n * (n * (n * 4210684958 / 1915538625))))) * sin(12 * phi)
        )

        if not radians:
            common_lat = common_lat * 180 / pi

        return common_lat
    """

exec_time = (
    timeit.timeit(
        stmt="auth_lat(phi=pi/4, e=0.08181919104281579, radians=False)", setup=setup_auth_lat_pr_new, number=100000
    )
    * 1000
)
print(f"Exec. time (auth_lat() PR branch new, degrees) = {exec_time}")

exec_time = (
    timeit.timeit(
        stmt="auth_lat(phi=pi/4, e=0.08181919104281579, radians=True)", setup=setup_auth_lat_pr_new, number=100000
    )
    * 1000
)
print(f"Exec. time (auth_lat() PR branch new, radians) = {exec_time}")



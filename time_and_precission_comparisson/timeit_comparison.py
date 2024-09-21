import math, timeit
import numpy as np


def auth_lat_old(
    phi: float, e: float, inverse: bool = False, radians: bool = False
) -> float:
    if e == 0:
        return phi
    if not radians:
        # Convert to radians to do calculations below.
        phi = np.deg2rad(phi)
    if not inverse:
        # Compute authalic latitude from latitude phi.
        q = ((1 - e**2) * np.sin(phi)) / (1 - (e * np.sin(phi)) ** 2) - (1 - e**2) / (
            2.0 * e
        ) * np.log((1 - e * np.sin(phi)) / (1 + e * np.sin(phi)))
        qp = 1 - (1 - e**2) / (2.0 * e) * np.log((1.0 - e) / (1.0 + e))
        ratio = q / qp
        # Avoid rounding errors.
        if abs(ratio) > 1:
            # Make abs(ratio) = 1
            ratio = np.sign(ratio)
        result = np.arcsin(ratio)
    else:
        # Compute an approximation of latitude from authalic latitude phi.
        result = (
            phi
            + (e**2 / 3.0 + 31 * e**4 / 180.0 + 517 * e**6 / 5040.0)
            * np.sin(2 * phi)
            + (23 * e**4 / 360.0 + 251 * e**6 / 3780.0) * np.sin(4 * phi)
            + (761 * e**6 / 45360.0) * np.sin(6 * phi)
        )
    if not radians:
        # Convert back to degrees.
        result = np.rad2deg(result)
    return result


latitudes = [x for x in range(-90,91)]

for x in latitudes:
    exec_time = (
            timeit.timeit(
                f"auth_lat_old(phi={x}, e=0.08181919104281579, radians=False)", globals=globals(), number=100000)
    )
    print(f"{x} {exec_time}")

latitudes = [x*np.pi/180 for x in range(-90,91)]

for x in latitudes:
    exec_time = (
            timeit.timeit(
                f"auth_lat_old(phi={x}, e=0.08181919104281579, radians=True)", globals=globals(), number=100000)
    )
    print(f"{x} {exec_time}")


def auth_lat_new(
    phi: float, e: float, inverse: bool = False, radians: bool = False
) -> float:
    if e == 0:
        return phi
    # Compute flattening f and third flattening n from eccentricity e.
    f = 1 - math.sqrt(1 - e**2)
    n = (1 - math.sqrt(1 - e**2)) / (1 + math.sqrt(1 - e**2))

    if not inverse:
        # Compute authalic latitude from common latitude phi.
        # For large flattenings (f > 1/150) use direct formula,
        # for small flattenings (f <= 1/150) use power series.
        if abs(f) > 1 / 150:
            # Use direct formula for large flattenings.
            if not radians:
                # Convert to radians to do calculations below.
                phi = phi * math.pi / 180
            # Compute authalic latitude from latitude phi.
            q = ((1 - e**2) * math.sin(phi)) / (1 - (e * math.sin(phi)) ** 2) - (1 - e**2) / (
                2.0 * e
            ) * math.log((1 - e * math.sin(phi)) / (1 + e * math.sin(phi)))
            qp = 1 - (1 - e**2) / (2.0 * e) * math.log((1.0 - e) / (1.0 + e))
            ratio = q / qp
            # Avoid rounding errors.
            if abs(ratio) > 1:
                # Make abs(ratio) = 1
                ratio = math.copysign(1, ratio)
            result = math.asin(ratio)
            if not radians:
                result = result * 180 / math.pi
            return result
        else:
            # Use power series approximation for small flattenings (f <= 1/150).
            # Power series expansion taken from https://doi.org/10.48550/arXiv.2212.05818 (Equation A19)
            if not radians:
                phi = phi * math.pi / 180

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
                * math.sin(2 * phi)
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
                * math.sin(4 * phi)
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
                * math.sin(6 * phi)
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
                * math.sin(8 * phi)
                + n
                * (n * (n * (n * (n * (-23356 / 66825 + n * (-839792 / 19348875))))))
                * math.sin(10 * phi)
                + n
                * (n * (n * (n * (n * (n * 570284222 / 1915538625)))))
                * math.sin(12 * phi)
            )

            if not radians:
                authalic_lat = authalic_lat * 180 / math.pi

            return authalic_lat
    else:
        # Compute common latitude from authalic latitude phi.
        # Power series expansion taken from https://doi.org/10.48550/arXiv.2212.05818 (Equation A20)
        if not radians:
            phi = phi * math.pi / 180

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
            * math.sin(2 * phi)
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
            * math.sin(4 * phi)
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
            * math.sin(6 * phi)
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
            * math.sin(8 * phi)
            + n
            * (n * (n * (n * (n * (768272 / 467775 + n * 455935736 / 638512875)))))
            * math.sin(10 * phi)
            + n * (n * (n * (n * (n * (n * 4210684958 / 1915538625))))) * math.sin(12 * phi)
        )

        if not radians:
            common_lat = common_lat * 180 / math.pi

        return common_lat


def auth_lat_new_clenshaw(
    phi: float, e: float, inverse: bool = False, radians: bool = False
) -> float:

    if e == 0:
        return phi
    # Compute flattening f and third flattening n from eccentricity e.
    f = 1 - math.sqrt(1 - e**2)
    n = (1 - math.sqrt(1 - e**2)) / (1 + math.sqrt(1 - e**2))

    if not inverse:
        # Compute authalic latitude from common latitude phi.
        # For large flattenings (f > 1/150) use direct formula,
        # for small flattenings (f <= 1/150) use power series.
        if abs(f) > 1 / 150:
            # Use direct formula for large flattenings.
            if not radians:
                # Convert to radians to do calculations below.
                phi = phi * math.pi / 180
            # Compute authalic latitude from latitude phi.
            q = ((1 - e**2) * math.sin(phi)) / (1 - (e * math.sin(phi)) ** 2) - (1 - e**2) / (
                2.0 * e
            ) * math.log((1 - e * math.sin(phi)) / (1 + e * math.sin(phi)))
            qp = 1 - (1 - e**2) / (2.0 * e) * math.log((1.0 - e) / (1.0 + e))
            ratio = q / qp
            # Avoid rounding errors.
            if abs(ratio) > 1:
                # Make abs(ratio) = 1
                ratio = math.copysign(1, ratio)
            result = math.asin(ratio)
            if not radians:
                result = result * 180 / math.pi
            return result
        else:
            # Use power series approximation for small flattenings (f <= 1/150).
            # Power series expansion taken from https://doi.org/10.48550/arXiv.2212.05818 (Equation A19)
            if not radians:
                phi = phi * math.pi / 180

            c0 = n * (-4 / 3 + n * (-4 / 45 + n * (88 / 315 + n * (538 / 4725 + n * (20824 / 467775 + n * (-44732 / 2837835))))))
            c1 = n * (n * (34 / 45 + n * (8 / 105 + n * (-2482 / 14175 + n * (-37192 / 467775 + n * (-12467764 / 212837625))))))
            c2 = n * (n * (n * (-1532 / 2835 + n * (-898 / 14175 + n * (54968 / 467775 + n * 100320856 / 1915538625)))))
            c3 = n * (n * (n * (n * (6007 / 14175 + n * (24496 / 467775 + n * (-5884124 / 70945875))))))
            c4 = n * (n * (n * (n * (n * (-23356 / 66825 + n * (-839792 / 19348875))))))
            c5 = n * (n * (n * (n * (n * (n * 570284222 / 1915538625)))))

            c = [c0, c1, c2, c3, c4, c5]

            authalic_lat = clenshaw(phi, c)

            if not radians:
                authalic_lat = authalic_lat * 180 / math.pi

            return authalic_lat
    else:
        # Compute common latitude from authalic latitude phi.
        # Power series expansion taken from https://doi.org/10.48550/arXiv.2212.05818 (Equation A20)
        if not radians:
            phi = phi * math.pi / 180

        c0 = n * (4 / 3 + n * (4 / 45 + n * (-16 / 35 + n * (-2582 / 14175 + n * (60136 / 467775 + n * 28112932 / 212837625)))))
        c1 = n * (n * (46 / 45 + n * (152 / 945 + n * (-11966 / 14175 + n * (-21016 / 51975 + n * 251310128 / 638512875)))))
        c2 = n * (n * (n * (3044 / 2835 + n * (3802 / 14175 + n * (-94388 / 66825 + n * (-8797648 / 10945935))))))
        c3 = n * (n * (n * (n * (6059 / 4725 + n * (41072 / 93555 + n * (-1472637812 / 638512875))))))
        c4 = n * (n * (n * (n * (n * (768272 / 467775 + n * 455935736 / 638512875)))))
        c5 = n * (n * (n * (n * (n * (n * 4210684958 / 1915538625)))))

        c = [c0, c1, c2, c3, c4, c5]

        common_lat = clenshaw(phi, c)

        if not radians:
            common_lat = common_lat * 180 / math.pi

        return common_lat


def clenshaw(phi: float, c: list) -> float:
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)
    x = (cos_phi + sin_phi) * (cos_phi - sin_phi)

    u5 = c[5]
    u4 = 2 * x * u5 + c[4]
    u3 = 2 * x * u4 - u5 + c[3]
    u2 = 2 * x * u3 - u4 + c[2]
    u1 = 2 * x * u2 - u3 + c[1]
    u0 = 2 * x * u1 - u2 + c[0]

    return phi + 2 * u0 * sin_phi * cos_phi


latitudes = [x for x in range(-90,91)]

for x in latitudes:
    exec_time = (
            timeit.timeit(
                f"auth_lat_new(phi={x}, e=0.08181919104281579, radians=False)", globals=globals(), number=100000)
    )
    print(f"{x} {exec_time}")

for x in latitudes:
    exec_time = (
            timeit.timeit(
                f"auth_lat_new(phi={x}, e=0.08181919104281579, radians=True)", globals=globals(), number=100000)
    )
    print(f"{x} {exec_time}")

for x in latitudes:
    exec_time = (
            timeit.timeit(
                f"auth_lat_new_clenshaw(phi={x}, e=0.08181919104281579, radians=False)", globals=globals(), number=100000)
    )
    print(f"{x} {exec_time}")

for x in latitudes:
    exec_time = (
            timeit.timeit(
                f"auth_lat_new_clenshaw(phi={x}, e=0.08181919104281579, radians=True)", globals=globals(), number=100000)
    )
    print(f"{x} {exec_time}")
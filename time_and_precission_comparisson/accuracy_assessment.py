import csv
import math

import mpmath as mp
import numpy as np

# set precision to 70 decimal places
mp.mp.dps = 70


def auth_lat_direct_mpmath(fi, e2):
    fi = mp.mpf(str(fi))
    e2 = mp.mpf(str(e2))
    return mp.degrees(
        mp.asin(
            (
                (1 - e2)
                * (
                    mp.sin(mp.radians(fi))
                    / (1 - e2 * mp.power(mp.sin(mp.radians(fi)), 2))
                    - 1
                    / (2 * mp.sqrt(e2))
                    * mp.log(
                        (1 - mp.sqrt(e2) * mp.sin(mp.radians(fi)))
                        / (1 + mp.sqrt(e2) * mp.sin(mp.radians(fi)))
                    )
                )
            )
            / (
                1
                - (1 - e2)
                / (2 * mp.sqrt(e2))
                * mp.log((1 - mp.sqrt(e2)) / (1 + mp.sqrt(e2)))
            )
        )
    )


def auth_lat_old(phi, e, inverse=False, radians=False):
    if e == 0:
        return phi
    if not radians:
        # Convert to radians to do calculations below.
        phi = np.deg2rad(phi)
    if not inverse:
        # Compute authalic latitude from latitude phi.
        q = ((1 - e**2) * np.sin(phi)) / (1 - (e * np.sin(phi)) ** 2) - (
            1 - e**2
        ) / (2.0 * e) * np.log((1 - e * np.sin(phi)) / (1 + e * np.sin(phi)))
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
            + (e**2 / 3.0 + 31 * e**4 / 180.0 + 517 * e**6 / 5040.0) * np.sin(2 * phi)
            + (23 * e**4 / 360.0 + 251 * e**6 / 3780.0) * np.sin(4 * phi)
            + (761 * e**6 / 45360.0) * np.sin(6 * phi)
        )
    if not radians:
        # Convert back to degrees.
        result = np.rad2deg(result)
    return result


def auth_lat_old_mpmath(phi, e, inverse=False, radians=False):
    phi = mp.mpf(str(phi))
    e = mp.mpf(str(e))
    if e == 0:
        return phi
    if not radians:
        # Convert to radians to do calculations below.
        phi = mp.radians(phi)
    if not inverse:
        # Compute authalic latitude from latitude phi.
        q = ((1 - mp.power(e, 2)) * mp.sin(phi)) / (
            1 - mp.power(e * mp.sin(phi), 2)
        ) - (1 - mp.power(e, 2)) / (2 * e) * mp.log(
            (1 - e * mp.sin(phi)) / (1 + e * mp.sin(phi))
        )
        qp = 1 - (1 - mp.power(e, 2)) / (2 * e) * mp.log((1 - e) / (1 + e))
        ratio = q / qp
        # Avoid rounding errors.
        if abs(ratio) > 1:
            # Make abs(ratio) = 1
            ratio = np.sign(ratio)
        result = mp.asin(ratio)
    else:
        # Compute an approximation of latitude from authalic latitude phi.
        result = (
            phi
            + (
                mp.power(e, 2) / 3
                + 31 * mp.power(e, 4) / 180
                + 517 * mp.power(e, 6) / 5040
            )
            * mp.sin(2 * phi)
            + (23 * mp.power(e, 4) / 360 + 251 * mp.power(e, 6) / 3780)
            * mp.sin(4 * phi)
            + (761 * mp.power(e, 6) / 45360) * mp.sin(6 * phi)
        )
    if not radians:
        # Convert back to degrees.
        result = mp.degrees(result)
    return result


def auth_lat_new(phi, e, inverse=False, radians=False):
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
            q = ((1 - e**2) * math.sin(phi)) / (1 - (e * math.sin(phi)) ** 2) - (
                1 - e**2
            ) / (2.0 * e) * math.log((1 - e * math.sin(phi)) / (1 + e * math.sin(phi)))
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
            + n
            * (n * (n * (n * (n * (n * 4210684958 / 1915538625)))))
            * math.sin(12 * phi)
        )

        if not radians:
            common_lat = common_lat * 180 / math.pi

        return common_lat


def auth_lat_new_mpmath(phi, e, inverse=False, radians=False):
    phi = mp.mpf(str(phi))
    e = mp.mpf(str(e))

    if e == mp.mpf("0"):
        return phi
    # Compute flattening f and third flattening n from eccentricity e.
    f = 1 - mp.sqrt(1 - mp.power(e, 2))
    n = (1 - mp.sqrt(1 - mp.power(e, 2))) / (1 + mp.sqrt(1 - mp.power(e, 2)))

    if not inverse:
        # Compute authalic latitude from common latitude phi.
        # For large flattenings (f > 1/150) use direct formula,
        # for small flattenings (f <= 1/150) use power series.
        if abs(f) > 1 / 150:
            # Use direct formula for large flattenings.
            if not radians:
                # Convert to radians to do calculations below.
                phi = mp.radians(phi)
            # Compute authalic latitude from latitude phi.
            q = ((1 - mp.power(e, 2)) * mp.sin(phi)) / (
                1 - mp.power(e * mp.sin(phi), 2)
            ) - (1 - mp.power(e, 2)) / (2 * e) * mp.log(
                (1 - e * mp.sin(phi)) / (1 + e * mp.sin(phi))
            )
            qp = 1 - (1 - mp.power(e)) / (2 * e) * mp.log((1 - e) / (1 + e))
            ratio = q / qp
            # Avoid rounding errors.
            if abs(ratio) > 1:
                # Make abs(ratio) = 1
                ratio = math.copysign(1, ratio)
            result = mp.asin(ratio)
            if not radians:
                result = mp.degrees(result)
            return result
        else:
            # Use power series approximation for small flattenings (f <= 1/150).
            # Power series expansion taken from https://doi.org/10.48550/arXiv.2212.05818 (Equation A19)
            if not radians:
                phi = mp.radians(phi)

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
                * mp.sin(2 * phi)
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
                * mp.sin(4 * phi)
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
                * mp.sin(6 * phi)
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
                * mp.sin(8 * phi)
                + n
                * (n * (n * (n * (n * (-23356 / 66825 + n * (-839792 / 19348875))))))
                * mp.sin(10 * phi)
                + n
                * (n * (n * (n * (n * (n * 570284222 / 1915538625)))))
                * mp.sin(12 * phi)
            )

            if not radians:
                authalic_lat = mp.degrees(authalic_lat)

            return authalic_lat
    else:
        # Compute common latitude from authalic latitude phi.
        # Power series expansion taken from https://doi.org/10.48550/arXiv.2212.05818 (Equation A20)
        if not radians:
            phi = mp.radians(phi)

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
            * mp.sin(2 * phi)
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
            * mp.sin(4 * phi)
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
            * mp.sin(6 * phi)
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
            * mp.sin(8 * phi)
            + n
            * (n * (n * (n * (n * (768272 / 467775 + n * 455935736 / 638512875)))))
            * mp.sin(10 * phi)
            + n
            * (n * (n * (n * (n * (n * 4210684958 / 1915538625)))))
            * mp.sin(12 * phi)
        )

        if not radians:
            common_lat = mp.degrees(common_lat)

        return common_lat


# WGS 84 (EPSG:7030)
# https://epsg.org/ellipsoid_7030/WGS-84.html
eccentricity = 0.08181919084262149
eccentricity_2 = 0.006694379990141316


common_latitudes = [x / 10 for x in range(-900, 901)]

with open("accuracy_results.csv", "w", newline="") as csvfile:
    csv_writer = csv.writer(csvfile, delimiter=",")
    csv_writer.writerow(
        [
            "Common latitude",
            "Authalic (direct mpmath)",
            "Authalic (OLD)",
            "Authalic (NEW)",
            "Diff (auth OLD - auth NEW) E-12",
            "Diff (auth OLD - auth direct mpmath) E-12",
            "Diff (auth NEW - auth direct mpmath) E-12",
            "Authalic inverse (OLD)",
            "Authalic inverse (NEW)",
            "Diff (auth inv OLD - auth inv NEW) E-12",
            "Diff (auth inv OLD - common) E-12",
            "Diff (auth inv NEW - common) E-12",
        ]
    )
    for phi in common_latitudes:
        auth_latitude_mpmath = auth_lat_direct_mpmath(phi, eccentricity_2)
        auth_latitude_old = auth_lat_old(phi, eccentricity)
        auth_latitude_new = auth_lat_new(phi, eccentricity)

        diff_old_new = str(
            (mp.mpf(str(auth_latitude_old)) - mp.mpf(str(auth_latitude_new))) * 10**12
        )
        diff_old_mpmath = str(
            (mp.mpf(str(auth_latitude_old)) - mp.mpf(str(auth_latitude_mpmath)))
            * 10**12
        )
        diff_new_mpmath = str(
            (mp.mpf(str(auth_latitude_new)) - mp.mpf(str(auth_latitude_mpmath)))
            * 10**12
        )

        auth_latitude_inverse_old = auth_lat_old(
            float(auth_latitude_mpmath), eccentricity, inverse=True
        )
        auth_latitude_inverse_new = auth_lat_new(
            float(auth_latitude_mpmath), eccentricity, inverse=True
        )

        diff_inverse_old_new = str(
            (
                mp.mpf(str(auth_latitude_inverse_old))
                - mp.mpf(str(auth_latitude_inverse_new))
            )
            * 10**12
        )
        diff_inverse_old_common = str(
            (mp.mpf(str(auth_latitude_inverse_old)) - mp.mpf(str(phi))) * 10**12
        )
        diff_inverse_new_common = str(
            (mp.mpf(str(auth_latitude_inverse_new)) - mp.mpf(str(phi))) * 10**12
        )

        write_list = [
            phi,
            auth_latitude_mpmath,
            auth_latitude_old,
            auth_latitude_new,
            diff_old_new,
            diff_old_mpmath,
            diff_new_mpmath,
            auth_latitude_inverse_old,
            auth_latitude_inverse_new,
            diff_inverse_old_new,
            diff_inverse_old_common,
            diff_inverse_new_common,
        ]

        csv_writer.writerow(write_list)

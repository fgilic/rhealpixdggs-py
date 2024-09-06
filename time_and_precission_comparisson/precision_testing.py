import csv
import math

import mpmath as mp
import numpy as np

# precision to 70 decimal places
mp.mp.dps = 70


def auth_lat_mpmath(fi, e2):
    fi = mp.mpf(str(fi))
    return mp.degrees(
        mp.asin(
            (
                (1 - e2)
                * (
                    mp.sin(mp.radians(fi)) / (1 - e2 * (mp.sin(mp.radians(fi))) ** 2)
                    - 1
                    / (2 * mp.sqrt(e2))
                    * mp.log(
                        (1 - mp.sqrt(e2) * mp.sin(mp.radians(fi)))
                        / (1 + mp.sqrt(e2) * mp.sin(mp.radians(fi)))
                    )
                )
            )
            / (1 - (1 - e2) / (2 * mp.sqrt(e2)) * mp.log((1 - mp.sqrt(e2)) / (1 + mp.sqrt(e2))))
        )
    )


def auth_lat_old(phi, e, inverse = False, radians = False):
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


def auth_lat_new(phi, e, inverse = False, radians = False):
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


# GRS 80 (EPSG:7019)
a = mp.mpf("6378137")
b = mp.mpf("6356752.31414")
e2 = (a**2 - b**2) / a**2
# e2 = mp.mpf("0.08181919104281579") ** 2


# GRS 80 (EPSG:7019)
a = 6378137
b = 6356752.31414
e = math.sqrt((a**2 - b**2) / a**2)


common_latitudes = [x/10 for x in range(-900,902)]
auth_latitudes_mpmath = [auth_lat_mpmath(x, e2) for x in common_latitudes]
auth_latitudes_old = [auth_lat_old(x, e) for x in common_latitudes]
auth_latitudes_new = [auth_lat_new(x, e) for x in common_latitudes]

auth_latitudes_inverse_old = [auth_lat_old(float(x), e, inverse=True) for x in auth_latitudes_mpmath]
auth_latitudes_inverse_new = [auth_lat_new(float(x), e, inverse=True) for x in auth_latitudes_mpmath]

with open("results_2.csv", "w", newline="") as csvfile:
    csv_writer = csv.writer(csvfile, delimiter=",")
    csv_writer.writerow(["Common latitude", "Authalic (mpmath)", "Authalic (old code)", "Authalic (new code)",
                         "Diff. (old - new) E-10", "Diff. (new - mpmath) E-10", "Diff. (old - mpmath)E-10",
                         "Authalic inverse (old code)", "Authalic inverse (new code)", "Diff. (auth_inv_new - common)E-10",
                         "Diff. (auth_inv_old - common)E-10"])
    for i in range(0, len(auth_latitudes_old) - 1):
        write_list = [common_latitudes[i], str(auth_latitudes_mpmath[i]), auth_latitudes_old[i], auth_latitudes_new[i]]

        diff_old_new = str((mp.mpf(str(auth_latitudes_old[i])) - mp.mpf(str(auth_latitudes_new[i])))*10000000000)
        diff_new_mpmath = str((mp.mpf(str(auth_latitudes_new[i])) - auth_latitudes_mpmath[i])*10000000000)
        diff_old_mpmath = str((mp.mpf(str(auth_latitudes_old[i])) - auth_latitudes_mpmath[i])*10000000000)

        diff_invnew_common = str((mp.mpf(str(auth_latitudes_inverse_new[i])) - mp.mpf(str(common_latitudes[i]))) * 10000000000)
        diff_invold_common = str((mp.mpf(str(auth_latitudes_inverse_old[i])) - mp.mpf(str(common_latitudes[i])))*10000000000)

        write_list.extend([diff_old_new, diff_new_mpmath, diff_old_mpmath, auth_latitudes_inverse_old[i], auth_latitudes_inverse_new[i], diff_invnew_common, diff_invold_common])

        csv_writer.writerow(write_list)

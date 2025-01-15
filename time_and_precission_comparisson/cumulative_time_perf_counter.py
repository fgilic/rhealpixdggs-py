import math, time, timeit
import numpy as np
import gc

# disable garbage collector
gc.disable()

def auth_lat_current_numpy(
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

def auth_lat_current_math(
    phi: float, e: float, inverse: bool = False, radians: bool = False
) -> float:
    if e == 0:
        return phi
    if not radians:
        # Convert to radians to do calculations below.
        phi = phi * math.pi / 180
    if not inverse:
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
    else:
        # Compute an approximation of latitude from authalic latitude phi.
        result = (
            phi
            + (e**2 / 3.0 + 31 * e**4 / 180.0 + 517 * e**6 / 5040.0)
            * math.sin(2 * phi)
            + (23 * e**4 / 360.0 + 251 * e**6 / 3780.0) * math.sin(4 * phi)
            + (761 * e**6 / 45360.0) * math.sin(6 * phi)
        )
    if not radians:
        # Convert back to degrees.
        result = result * 180 / math.pi
    return result


def auth_lat_modified(phi: float, e: float, inverse: bool = False, radians: bool = False) -> float:
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


# for points in num_points:
#     latitudes = np.linspace(-90, 90, points)
#     cumulative_time = 0.0
#     for x in latitudes:
#         exec_time = (
#                 timeit.timeit(
#                     f"auth_lat_new(phi={x}, e=0.08181919104281579, radians=False)", globals=globals(), number=1)
#         )
#         cumulative_time += exec_time
#     print(f"{points} {cumulative_time}")




# num_points = [180+1, 1800+1,18000+1, 180000+1, 18000000+1]
num_points = [180+1, 1080+1,10800+1, 64800+1, 648000+1, 6480000+1]

# for points in num_points:
#     latitudes = np.linspace(-90, 90, points)
#     cumulative_time = 0.0
#     for x in latitudes:
#         exec_time = (
#                 timeit.timeit(
#                     f"auth_lat_current_numpy(phi={x}, e=0.08181919084262149, radians=False)", globals=globals(), number=1)
#         )
#         cumulative_time += exec_time
#     print(f"{points} {cumulative_time}")




#########################################START#########################################
# current numpy
print("Current NumPy")
for points in num_points:
    latitudes = np.linspace(-90, 90, points)

    if len(latitudes) < 70000:
        # if less than 70000 input latitudes repeat calculation 1000x and take minimum time
        cumulative_times_list = []
        for i in range(1000):
            time_start = time.perf_counter()
            for x in latitudes:
                auth_lat_current_numpy(phi=x, e=0.08181919084262149, inverse=False, radians=False)
            time_end = time.perf_counter()
            time_duration = time_end - time_start
            cumulative_times_list.append(time_duration)
        time_duration = min(cumulative_times_list)
    else:
        time_start = time.perf_counter()
        for x in latitudes:
            auth_lat_current_numpy(phi=x, e=0.08181919084262149, inverse=False, radians=False)
        time_end = time.perf_counter()
        time_duration = time_end - time_start
    print(f"{points} {time_duration}")
#
# current math
print("\nCurrent math")
for points in num_points:
    latitudes = np.linspace(-90, 90, points)

    if len(latitudes) < 70000:
        # if less than 70000 input latitudes repeat calculation 1000x and take minimum time
        cumulative_times_list = []
        for i in range(1000):
            time_start = time.perf_counter()
            for x in latitudes:
                auth_lat_current_math(phi=x, e=0.08181919084262149, inverse=False, radians=False)
            time_end = time.perf_counter()
            time_duration = time_end - time_start
            cumulative_times_list.append(time_duration)
        time_duration = min(cumulative_times_list)
    else:
        time_start = time.perf_counter()
        for x in latitudes:
            auth_lat_current_math(phi=x, e=0.08181919084262149, inverse=False, radians=False)
        time_end = time.perf_counter()
        time_duration = time_end - time_start
    print(f"{points} {time_duration}")
#
# modified
print("\nModified")
for points in num_points:
    latitudes = np.linspace(-90, 90, points)

    if len(latitudes) < 70000:
        # if less than 70000 input latitudes repeat calculation 1000x and take minimum time
        cumulative_times_list = []
        for i in range(1000):
            time_start = time.perf_counter()
            for x in latitudes:
                auth_lat_modified(phi=x, e=0.08181919084262149, inverse=False, radians=False)
            time_end = time.perf_counter()
            time_duration = time_end - time_start
            cumulative_times_list.append(time_duration)
        time_duration = min(cumulative_times_list)
    else:
        time_start = time.perf_counter()
        for x in latitudes:
            auth_lat_modified(phi=x, e=0.08181919084262149, inverse=False, radians=False)
        time_end = time.perf_counter()
        time_duration = time_end - time_start
    print(f"{points} {time_duration}")
#
# current inverse numpy
print("\nCurrent inverse NumPy")
for points in num_points:
    latitudes = np.linspace(-90, 90, points)

    if len(latitudes) < 70000:
        # if less than 70000 input latitudes repeat calculation 1000x and take minimum time
        cumulative_times_list = []
        for i in range(1000):
            time_start = time.perf_counter()
            for x in latitudes:
                auth_lat_current_numpy(phi=x, e=0.08181919084262149, inverse=True, radians=False)
            time_end = time.perf_counter()
            time_duration = time_end - time_start
            cumulative_times_list.append(time_duration)
        time_duration = min(cumulative_times_list)
    else:
        time_start = time.perf_counter()
        for x in latitudes:
            auth_lat_current_numpy(phi=x, e=0.08181919084262149, inverse=True, radians=False)
        time_end = time.perf_counter()
        time_duration = time_end - time_start
    print(f"{points} {time_duration}")

# current inverse math
print("\nCurrent inverse math")
for points in num_points:
    latitudes = np.linspace(-90, 90, points)

    if len(latitudes) < 70000:
        # if less than 70000 input latitudes repeat calculation 1000x and take minimum time
        cumulative_times_list = []
        for i in range(1000):
            time_start = time.perf_counter()
            for x in latitudes:
                auth_lat_current_math(phi=x, e=0.08181919084262149, inverse=True, radians=False)
            time_end = time.perf_counter()
            time_duration = time_end - time_start
            cumulative_times_list.append(time_duration)
        time_duration = min(cumulative_times_list)
    else:
        time_start = time.perf_counter()
        for x in latitudes:
            auth_lat_current_math(phi=x, e=0.08181919084262149, inverse=True, radians=False)
        time_end = time.perf_counter()
        time_duration = time_end - time_start
    print(f"{points} {time_duration}")


# modified inverse
print("\nModified inverse")
for points in num_points:
    latitudes = np.linspace(-90, 90, points)

    if len(latitudes) < 70000:
        # if less than 70000 input latitudes repeat calculation 1000x and take minimum time
        cumulative_times_list = []
        for i in range(1000):
            time_start = time.perf_counter()
            for x in latitudes:
                auth_lat_modified(phi=x, e=0.08181919084262149, inverse=True, radians=False)
            time_end = time.perf_counter()
            time_duration = time_end - time_start
            cumulative_times_list.append(time_duration)
        time_duration = min(cumulative_times_list)
    else:
        time_start = time.perf_counter()
        for x in latitudes:
            auth_lat_modified(phi=x, e=0.08181919084262149, inverse=True, radians=False)
        time_end = time.perf_counter()
        time_duration = time_end - time_start
    print(f"{points} {time_duration}")
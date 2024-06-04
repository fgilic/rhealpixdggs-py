import math
from rhealpixdggs.utils import auth_rad, auth_lat
import numpy as np

a = 6378137
b = 6356752.314140356
fi = 1

e_2 = (a ** 2 - b ** 2) / a ** 2
e = math.sqrt((a ** 2 - b ** 2) / a ** 2)
e_crtano_2 = (a ** 2 - b ** 2) / b ** 2
c = a ** 2 / b
n = (a - b) / (a + b)

# authalic radius
R_auth_moritz = c * (
    1
    - 2 / 3 * e_crtano_2
    + 26 / 45 * e_crtano_2**2
    - 100 / 189 * e_crtano_2**3
    + 7034 / 14175 * e_crtano_2**4
)
R_auth_snyder = a * math.sqrt(
    (1 - (1 - e_2) / (2 * e) * math.log((1 - e) / (1 + e))) / 2
)
# Karnely from (13)
R_auth_karney = math.sqrt(a ** 2 / 2 + b ** 2 / 2 * math.atanh(e) / e)
# Francula from (3.77)
R_auth_francula = a * (1 - 2 / 3 * n + 26 / 45 * n ** 2 - 374 / 945 * n ** 3 + 722 / 2025 * n ** 4)
print("Authalic radius:")
print(f"Moritz = {R_auth_moritz}")
print(f"Snyder = {R_auth_snyder}")
print(f"rHEALPix (Snyder) = {auth_rad(a, e)}")
print(f"Karney (13) = {R_auth_karney}")
print(f"Francula = {R_auth_francula}")


# authalic latitude
# Snyder
q = (1 - e_2) * (
    math.sin(fi) / (1 - e_2 * (math.sin(fi)) ** 2)
    - 1 / (2 * e) * math.log((1 - e * math.sin(fi)) / (1 + e * math.sin(fi)))
)
q_p = 1 - (1 - e_2) / (2 * e) * math.log((1 - e) / (1 + e))
auth_lat_snyder = math.asin(q / q_p)

#China
auth_lat_china = (
    fi
    - 1 / 3 * (e_2 + 31 / 60 * e_2) * math.sin(2 * fi)
    + 17 / 360 * e_2**2 * math.sin(4 * fi)
)

#rHEALPix
auth_lat_rhealpix = auth_lat(fi, e, radians=True)


# Karney
# c_auth_to_phi = [
#     [4 / 3, 4 / 45, -16 / 35, -2582 / 14175, 60136 / 467775, 28112932 / 212837625],
#     [0, 46 / 45, 152 / 945, -11966 / 14175, -21016 / 51975, 251310128 / 638512875],
#     [0, 0, 3044 / 2835, 3802 / 14175, -94388 / 66825, -8797648 / 10945935],
#     [0, 0, 0, 6059 / 4725, 41072 / 93555, -1472637812 / 638512875],
#     [0, 0, 0, 0, 768272 / 467775, 455935736 / 638512875],
#     [0, 0, 0, 0, 0, 4210684958 / 1915538625],
# ]


# c_phi_to_auth = [
#     [-4 / 3, -4 / 45, 88 / 315, 538 / 4725, 20824 / 467775, -44732 / 2837835],
#     [0, 34 / 45, 8 / 105, -2482 / 14175, -37192 / 467775, -12467764 / 212837625],
#     [0, 0, -1532 / 2835, -898 / 14175, 54968 / 467775, 100320856 / 1915538625],
#     [0, 0, 0, 6007 / 14175, 24496 / 467775, -5884124 / 70945875],
#     [0, 0, 0, 0, -23356 / 66825, -839792 / 19348875],
#     [0, 0, 0, 0, 0, 570284222 / 1915538625],
# ]
#
# mat_s = [math.sin(2 * fi), math.sin(4 * fi), math.sin(6 * fi), math.sin(8 * fi), math.sin(10 * fi), math.sin(12 * fi)]
#
# mat_p = np.transpose([n, n**2, n**3, n**4, n**5, n**6])

auth_lat_karney = fi + (
    (- 4 / 3 * n - 4 / 45 * n ** 2 + 88 / 315 * n ** 3 + 538 / 4725 * n ** 4 + 20824 / 467775 * n ** 5 - 44732 / 2837835 * n ** 6) * math.sin(2 * fi) +
    (34 / 45 * n ** 2 + 8 / 105 * n ** 3 - 2482 / 14175 * n ** 4 - 37192 / 467775 * n ** 5 - 12467764 / 212837625 * n ** 6) * math.sin(4 * fi) +
    (-1532 / 2835 * n ** 3 - 898 / 14175 * n ** 4 + 54968 / 467775 * n ** 5 + 100320856 / 1915538625 * n ** 6) * math.sin(6 * fi) +
    (6007 / 14175 * n ** 4 + 24496 / 467775 * n ** 5 - 5884124 / 70945875 * n ** 6) * math.sin(8 * fi) +
    (-23356 / 66825 * n ** 5 - 839792 / 19348875 * n ** 6) * math.sin(10 * fi) +
    (570284222 / 1915538625 * n ** 6) * math.sin(12 * fi)
    )

print("\nAuthalic latitudes:")
print(f"fi common = {math.degrees(fi)}")
print("---")
print(f"fi auth (Snyder) = {math.degrees(auth_lat_snyder)}")
print(f"fi auth (rHEALPix) = {math.degrees(auth_lat_rhealpix)}")
print(f"fi auth (China) = {math.degrees(auth_lat_china)}")
print(f"fi auth (Karney) = {math.degrees(auth_lat_karney)}")

# invert computed authalic latitudes back to common
common_lat_snyder = auth_lat_snyder + (e_2 / 3 + 31 * e_2 ** 2 / 180 + 517 * e_2 ** 3 / 5040) * math.sin(2 * auth_lat_snyder) + (23 * e_2 ** 2 / 360 + 251 * e_2 ** 3 / 3780) * math.sin(4 * auth_lat_snyder) + (761 * e_2 **3 /45360) * math.sin(6 * auth_lat_snyder)
common_lat_rhealpix = auth_lat(auth_lat_rhealpix, e, radians=True, inverse=True)

# mat_s = [math.sin(2 * auth_lat_karney), math.sin(4 * auth_lat_karney), math.sin(6 * auth_lat_karney), math.sin(8 * auth_lat_karney), math.sin(10 * auth_lat_karney), math.sin(12 * auth_lat_karney)]
# common_lat_karney = auth_lat_karney + np.matmul(np.matmul(mat_s, c_auth_to_phi), mat_p)

common_lat_karney = auth_lat_karney + (
    (4 / 3 * n + 4 / 45 * n ** 2 - 16 / 35 * n ** 3 - 2582 / 14175 * n ** 4 + 60136 / 467775 * n ** 5 + 28112932 / 212837625 * n ** 6) * math.sin(2 * auth_lat_karney) +
    (46 / 45 * n ** 2 + 152 / 945 * n ** 3 - 11966 / 14175 * n ** 4 - 21016 / 51975 * n ** 5 + 251310128 / 638512875 * n ** 6) * math.sin(4 * auth_lat_karney) +
    (3044 / 2835 * n ** 3 + 3802 / 14175 * n ** 4 - 94388 / 66825 * n ** 5 - 8797648 / 10945935 * n ** 6) * math.sin(6 * auth_lat_karney) +
    (6059 / 4725 * n ** 4 + 41072 / 93555 * n ** 5 - 1472637812 / 638512875 * n ** 6) * math.sin(8 * auth_lat_karney) +
    (768272 / 467775 * n ** 5 + 455935736 / 638512875 * n ** 6) * math.sin(10 * auth_lat_karney) +
    (4210684958 / 1915538625 * n ** 6) * math.sin(12 * auth_lat_karney)
)

print("\nInverted authalic latitudes back to common:")
print(f"fi common (Snyder) = {math.degrees(common_lat_snyder)}")
print(f"fi common (rHEALPix) = {math.degrees(common_lat_rhealpix)}")
print(f"fi common (Karney) = {math.degrees(common_lat_karney)}")



# Calculating common latitude from authalic
fi = math.radians(45) # authalic

common_lat_snyder = fi + (e_2 / 3 + 31 * e_2 ** 2 / 180 + 517 * e_2 ** 3 / 5040) * math.sin(2 * fi) + (23 * e_2 ** 2 / 360 + 251 * e_2 ** 3 / 3780) * math.sin(4 * fi) + (761 * e_2 **3 /45360) * math.sin(6 * fi)
common_lat_rhealpix = auth_lat(fi, e, radians=True, inverse=True)

# mat_s = [math.sin(2 * fi), math.sin(4 * fi), math.sin(6 * fi), math.sin(8 * fi), math.sin(10 * fi), math.sin(12 * fi)]
# common_lat_karney = fi + np.matmul(np.matmul(mat_s, c_auth_to_phi), mat_p)

common_lat_karney = fi + (
    (4 / 3 * n + 4 / 45 * n ** 2 - 16 / 35 * n ** 3 - 2582 / 14175 * n ** 4 + 60136 / 467775 * n ** 5 + 28112932 / 212837625 * n ** 6) * math.sin(2 * fi) +
    (46 / 45 * n ** 2 + 152 / 945 * n ** 3 - 11966 / 14175 * n ** 4 - 21016 / 51975 * n ** 5 + 251310128 / 638512875 * n ** 6) * math.sin(4 * fi) +
    (3044 / 2835 * n ** 3 + 3802 / 14175 * n ** 4 - 94388 / 66825 * n ** 5 - 8797648 / 10945935 * n ** 6) * math.sin(6 * fi) +
    (6059 / 4725 * n ** 4 + 41072 / 93555 * n ** 5 - 1472637812 / 638512875 * n ** 6) * math.sin(8 * fi) +
    (768272 / 467775 * n ** 5 + 455935736 / 638512875 * n ** 6) * math.sin(10 * fi) +
    (4210684958 / 1915538625 * n ** 6) * math.sin(12 * fi)
)

print("\nCommon latitude from authalic:")
print(f"fi authalic = {math.degrees(fi)}")
print("---")
print(f"fi common (Snyder) = {math.degrees(common_lat_snyder)}")
print(f"fi common (rHEALPix) = {math.degrees(common_lat_rhealpix)}")
print(f"fi common (Karney) = {math.degrees(common_lat_karney)}")


# arrea of elipsoidal quadrangle (bounded by two lines of latitude and two lines of longitude
la_1 = math.radians(0)
la_2 = math.radians(90)
fi_1 = math.radians(0)
fi_2 = math.radians(90)
P = (
    b**2
    / 2
    * (la_2 - la_1)
    * (
        (
            math.sin(fi_1) / (1 - e_2 * (math.sin(fi_1)) ** 2)
            + 1
            / (2 * e)
            * math.log((1 + e * math.sin(fi_1)) / (1 - e * math.sin(fi_1)))
        )
        - (
            math.sin(fi_2) / (1 - e_2 * (math.sin(fi_2)) ** 2)
            + 1
            / (2 * e)
            * math.log((1 + e * math.sin(fi_2)) / (1 - e * math.sin(fi_2)))
        )
    )
)
print(f'Area: {P}')
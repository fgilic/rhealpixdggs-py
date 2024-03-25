import math
from rhealpixdggs.ellipsoids import auth_rad, auth_lat

a = 6378137
b = 6356752.3141
fi = math.radians(6)

e_2 = (a**2 - b**2) / a**2
e = math.sqrt((a**2 - b**2) / a**2)
e_crtano_2 = (a**2 - b**2) / b**2
c = a**2 / b

R_auth_moritz = c * (1 - 2 / 3 * e_crtano_2 + 26 / 45 * e_crtano_2 ** 2 - 100 / 189 * e_crtano_2 ** 3 + 7034 / 14175 * e_crtano_2 ** 4)
R_auth_snyder = a * math.sqrt((1 - (1 - e_2) / (2 * e) * math.log((1 - e) / (1 + e))) / 2)

print(f'Moritz = {R_auth_moritz}')
print(f'Snyder = {R_auth_snyder}')
print(f'rHEALPix (Snyder) = {auth_rad(a, e)}/n')

q = (1 - e_2) * (math.sin(fi) / (1 - e_2 * (math.sin(fi))**2) - 1 / (2 * e) * math.log((1 - e * math.sin(fi)) / (1 + e * math.sin(fi))))
q_p = 1 - (1 - e_2) / (2 * e) * math.log((1 - e) / (1 + e))
auth_lat_snyder = math.asin(q / q_p)
auth_lat_china = fi - 1 / 3 * (e_2 + 31 / 60 * e_2) * math.sin(2 * fi) + 17/360 * e_2**2 * math.sin(4 * fi)

print(f'fi = {math.degrees(fi)}')
print(f'fi auth (Snyder) = {math.degrees(auth_lat_snyder)}')
print(f'fi auth (China) = {math.degrees(auth_lat_china)}')
print(f'fi auth (rHEALPix) = {math.degrees(auth_lat(fi, e, radians=True))}')

la_1 = math.radians(0)
la_2 = math.radians(90)
fi_1 = math.radians(0)
fi_2 = math.radians(90)
P = b**2 / 2 * (la_2 - la_1) * ((math.sin(fi_1) / (1 - e_2 * (math.sin(fi_1))**2) + 1 / (2*e) * math.log((1 + e * math.sin(fi_1)) / (1 - e * math.sin(fi_1)))) - (math.sin(fi_2) / (1 - e_2 * (math.sin(fi_2))**2) + 1 / (2*e) * math.log((1 + e * math.sin(fi_2)) / (1 - e * math.sin(fi_2)))))
print(P)

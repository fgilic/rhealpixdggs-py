import timeit

setup1 = "import math"

code1 = """
fi = 1
e_2 = 0.00669438003551279089034
e = 0.08181919111988819270999832
q = (1 - e_2) * (
    math.sin(fi) / (1 - e_2 * (math.sin(fi)) ** 2)
    - 1 / (2 * e) * math.log((1 - e * math.sin(fi)) / (1 + e * math.sin(fi)))
)
q_p = 1 - (1 - e_2) / (2 * e) * math.log((1 - e) / (1 + e))
auth_lat_snyder = math.asin(q / q_p)
"""

exec_time1= timeit.timeit(stmt=code1, setup=setup1, number=100000)*1000

setup2 = "import math, numpy as np"

code2 = ("""
fi = 1
n = 0.0016792203978029802

c_phi_to_auth = [
    [-4 / 3, -4 / 45, 88 / 315, 538 / 4725, 20824 / 467775, -44732 / 2837835],
    [0, 34 / 45, 8 / 105, -2482 / 14175, -37192 / 467775, -12467764 / 212837625],
    [0, 0, -1532 / 2835, -898 / 14175, 54968 / 467775, 100320856 / 1915538625],
    [0, 0, 0, 6007 / 14175, 24496 / 467775, -5884124 / 70945875],
    [0, 0, 0, 0, -23356 / 66825, -839792 / 19348875],
    [0, 0, 0, 0, 0, 570284222 / 1915538625],
]

mat_s = [math.sin(2 * fi), math.sin(4 * fi), math.sin(6 * fi), math.sin(8 * fi), math.sin(10 * fi), math.sin(12 * fi)]

mat_p = np.transpose([n, n**2, n**3, n**4, n**5, n**6])

auth_lat_karney = fi + np.matmul(np.matmul(mat_s, c_phi_to_auth), mat_p)
""")

exec_time2= timeit.timeit(stmt=code2, setup=setup2, number=100000)*1000


setup3 = "import math"

code3 = ("""
fi = 1
n = 0.0016792203978029802

auth_lat_karney =  fi + (
    (- 4 / 3 * n - 4 / 45 * n ** 2 + 88 / 315 * n ** 3 + 538 / 4725 * n ** 4 + 20824 / 467775 * n ** 5 - 44732 / 2837835 * n ** 6) * math.sin(2 * fi) +
    (34 / 45 * n ** 2 + 8 / 105 * n ** 3 - 2482 / 14175 * n ** 4 - 37192 / 467775 * n ** 5 - 12467764 / 212837625 * n ** 6) * math.sin(4 * fi) + 
    (-1532 / 2835 * n ** 3 - 898 / 14175 * n ** 4 + 54968 / 467775 * n ** 5 + 100320856 / 1915538625 * n ** 6) * math.sin(6 * fi) +
    (6007 / 14175 * n ** 4 + 24496 / 467775 * n ** 5 - 5884124 / 70945875 * n ** 6) * math.sin(8 * fi) +
    (-23356 / 66825 * n ** 5 - 839792 / 19348875 * n ** 6) * math.sin(10 * fi) +
    (570284222 / 1915538625 * n ** 6) * math.sin(12 * fi)
    )
""")

exec_time3 = timeit.timeit(stmt=code3, setup=setup3, number=100000)*1000


print(f"Exec time 1 = {exec_time1} ms\nExec time 2 = {exec_time2} ms\nExec time 3 = {exec_time3} ms")
# import math
# import numpy as np
#
# def get_level0_id(point_xyz):
#     y, y, z = point_xyz
#
# def rotate_around_z_axes(point_xyz, angle = math.pi / 4, radians = True):
#     x, y, z = point_xyz
#     if not radians:
#         angle = angle * math.pi / 180.0
#     rotation_matrix = np.array([[math.cos(angle), -math.sin(angle), 0],
#                                 [math.sin(angle), math.cos(angle), 0],
#                                 [0, 0, 1]])
#     point_vector = np.array([x, y, z])
#     return np.matmul(rotation_matrix, point_vector)
#

import math
import numpy as np

def qpix_sphere(lam, phi):
    if abs(phi) <= math.pi / 4:
        if -math.pi / 4 <= lam <= math.pi / 4:


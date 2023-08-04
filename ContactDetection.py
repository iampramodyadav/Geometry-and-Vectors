import numpy as np

check_collinear = lambda p1, p2, p3: np.linalg.norm(np.cross((p2-p1), (p1-p3)))
check_collinear.__doc__ ="""
    Check if three points are collinear.

    Args:
        p1 (numpy array): First point in 3D space (x, y, z).
        p2 (numpy array): Second point in 3D space (x, y, z).
        p3 (numpy array): Third point in 3D space (x, y, z).

    Returns:
        float: The magnitude of the cross product of vectors formed by (p2 - p1) and (p1 - p3).
              Returns a value close to zero if the points are collinear.
    """

check_coplaner = lambda p1, p2, p3, p4: np.linalg.det(np.array([ [p1[0], p2[0], p3[0], p4[0]],
                                                                  [p1[1], p2[1], p3[1], p4[1]],
                                                                  [p1[2], p2[2], p3[2], p4[2]],
                                                                  [1,     1,      1,     1  ] ]))
check_coplaner.__doc__ ="""
    Check if four points are coplanar.

    Args:
        p1 (numpy array): First point in 3D space (x, y, z).
        p2 (numpy array): Second point in 3D space (x, y, z).
        p3 (numpy array): Third point in 3D space (x, y, z).
        p4 (numpy array): Fourth point in 3D space (x, y, z).

    Returns:
        float: The determinant of a 4x4 matrix formed by the given points.
              Returns a value close to zero if the points are coplanar.
    """
dcv = lambda x, y, z: np.array([x / np.sqrt(x**2 + y**2 + z**2),
                                y / np.sqrt(x**2 + y**2 + z**2),
                                z / np.sqrt(x**2 + y**2 + z**2)])
dcv.__doc__ ="""
    Calculate the direction cosines of a point.

    Args:
        x (float): x-coordinate of the point.
        y (float): y-coordinate of the point.
        z (float): z-coordinate of the point.

    Returns:
        numpy array: Array containing the direction cosines of the point.
    """
check_normal = lambda v1, v2, v3: np.abs(np.array([np.dot(v1,v2),
                                                   np.dot(v2,v3),
                                                   np.dot(v3,v1)])) < np.array([1.0e-6,
                                                                                1.0e-6,
                                                                                1.0e-6])
check_normal.__doc__ = """
    Check if three vectors are perpendicular to each other.

    Args:
        v1 (numpy array): First vector in 3D space (x, y, z).
        v2 (numpy array): Second vector in 3D space (x, y, z).
        v3 (numpy array): Third vector in 3D space (x, y, z).

    Returns:
        numpy array: Boolean array indicating whether each pair of vectors is perpendicular.
    """
def coordinate_dcm(origin,p1,p2):
    """
    Calculate the Direction Cosine Matrix (DCM) of three perpendicular vectors.

    Args:
        origin (numpy array): The centroid (origin) in 3D space (x, y, z).
        p1 (numpy array): First point in 3D space (x, y, z).
        p2 (numpy array): Second point in 3D space (x, y, z).

    Returns:
        tuple: A tuple containing the Direction Cosine Matrix (DCM) of three perpendicular vectors (v1, v4, v3).
              v1 (numpy array): Direction cosine of vector v1.
              v4 (numpy array): Direction cosine of vector v4.
              v3 (numpy array): Direction cosine of vector v3.
    """
    if abs(check_collinear(origin,p1,p2)) < 1.0e-6:
        print("Points are collinear, Select largest line and create a kp")
    else:
        v1= p1- origin                 #vector-1
        v2= p2- origin                 #vector-2
        v1_dc=dcv(v1[0],v1[1],v1[2])
        v2_dc=dcv(v2[0],v2[1],v2[2])

        v3=np.cross(v1_dc, v2_dc)      #vector-3, perp to 1,2
        v3_dc=dcv(v3[0],v3[1],v3[2])

        v4=np.cross(v3_dc, v1_dc)      #vector-4(y) perp vector-3(z) perp to vector-1(x)
        v4_dc=dcv(v4[0],v4[1],v4[2])

        return v1_dc, v4_dc, v3_dc

def dcm2angleZXY(R):
    """
    Convert Direction Cosine Matrix (DCM) to Euler rotation angles in ZXY order.

    Args:
        R (numpy array): The Direction Cosine Matrix (DCM).

    Returns:
        tuple: A tuple containing the Euler rotation angles in ZXY order (Rx, Ry, Rz).
              Rx (float): Rotation angle around the x-axis (in degrees).
              Ry (float): Rotation angle around the y-axis (in degrees).
              Rz (float): Rotation angle around the z-axis (in degrees).
    """
    Ry = np.arctan2(-R[0, 2], R[2, 2])
    Rx = np.arcsin(R[1, 2])
    Rz = np.arctan2(-R[1, 0], R[1, 1])
    return np.degrees(Rx), np.degrees(Ry) ,np.degrees(Rz)


def circle_center_radius(A, B, C):
    """
    Calculate the circumcenter, radius, and center coordinates of a circle passing through three points.

    Args:
        A (numpy array): First point in 3D space (x, y, z).
        B (numpy array): Second point in 3D space (x, y, z).
        C (numpy array): Third point in 3D space (x, y, z).

    Returns:
        tuple: A tuple containing the center coordinates (x, y, z) and the radius of the circle.
              centr (numpy array): The center coordinates (x, y, z) of the circle.
              rad (float): The radius of the circle.
    Raises:
        ValueError: If the points are collinear and no unique circle exists.
    """
    if abs(check_collinear(A,B,C)) < 1.0e-6:
        raise ValueError("Points are collinear, no unique circle exists.")
    a = np.linalg.norm(C - B)
    b = np.linalg.norm(C - A)
    c = np.linalg.norm(B - A)
    s = (a + b + c) / 2
    rad = a*b*c / 4 / np.sqrt(s * (s - a) * (s - b) * (s - c))
     #Barycentric Coordinates of circumsnter
    b1 = a*a * (b*b + c*c - a*a)
    b2 = b*b * (a*a + c*c - b*b)
    b3 = c*c * (a*a + b*b - c*c)
     #Barycentric Coordinates to cartesian Coordinates
    centr = np.column_stack((A, B, C)).dot(np.hstack((b1, b2, b3)))
    centr /= b1 + b2 + b3
#     centr = centr/(b1 + b2 + b3)
    return centr,rad

def perpendicular_vector(vector):
    """
    Return two vectors, such that all three vectors perpendicular each other

    Args:
      vector: A 3D vector (numpy array) (if Z-direction)

    Returns: 

      two perpendicular vectors(numpy array) (Y-direction, X-direction)
    """

    # Take the cross product of the given vector and any other vector
    perp_vec = np.cross(vector, np.array([1, 0, 0]))
    mod_vec = np.linalg.norm(perp_vec)

    if mod_vec<1.0e-4:
        perp_vec = np.cross(vector, np.array([0, 1, 0]))
        mod_vec = np.linalg.norm(perp_vec)

    # Normalize the resulting vector.
    perp_vec_y = perp_vec / mod_vec
    perp_vec_x = np.cross( perp_vec_y, vector)

    return perp_vec_y, perp_vec_x
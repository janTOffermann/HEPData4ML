import numpy as np

def GetRotationMatrix(alpha,beta,gamma):
    """
    Given yaw (x-y), pitch (x-z) and roll (y-z),
    determine the rotation matrix.
    """
    R_alpha = np.array([
        [np.cos(alpha), -np.sin(alpha), 0.],
        [np.sin(alpha),  np.cos(alpha), 0.],
        [           0.,             0., 1.]
        ])

    R_beta = np.array([
        [ np.cos(beta), 0., np.sin(beta)],
        [           0., 1.,           0.],
        [-np.sin(beta), 0., np.cos(beta)]
        ])

    R_gamma = np.array([
        [1.,            0.,             0.],
        [0., np.cos(gamma), -np.sin(gamma)],
        [0., np.sin(gamma),  np.cos(gamma)]
        ])

    R = np.dot(R_alpha, np.dot(R_beta,R_gamma))
    return R

def RotateVector(v, alpha,beta,gamma):
    """
    Rotate 4-vector v with the given (yaw,pitch,roll).
    Using convention v = (E, px, py, pz).
    """
    v_new = np.copy(v)
    R = GetRotationMatrix(alpha,beta,gamma)
    v_new[1:4] = np.dot(R,v_new[1:4])
    return v_new

def GetRelRotationMatrix(vec1, vec2):
    """
    Get the rotation matrix that rotates vec1 into vec2.
    """
    # Get the associated unit vectors.
    v1 = vec1 / np.linalg.norm(vec1)
    v2 = vec2 / np.linalg.norm(vec2)

    q = v2 - np.dot(v1,v2) * v1
    q /= np.linalg.norm(q)

    Finv =  np.array([v1, q, np.cross(v2,v1)]).T
    F = np.linalg.inv(Finv)

    G = np.array([
        [                  np.dot(v1,v2), -np.linalg.norm(np.cross(v1,v2)), 0.],
        [np.linalg.norm(np.cross(v1,v2)),                    np.dot(v1,v2), 0.],
        [                             0.,                               0., 1.]
    ])

    U = np.dot(Finv,np.dot(G,F))
    return U

def DetermineRotation(vec1, vec2):
    """
    Return the rotation matrix necessary to rotate the vectors
    such that vec1 is along the +x-axis, and vec2 is in the x-y
    plane (in a particular direction).
    """
    vec_i = np.copy(vec1)
    if(len(vec1) == 4): vec_i = vec_i[1:4]
    axis = np.array([1.,0.,0.])
    R1 = GetRelRotationMatrix(vec_i,axis)

    if(vec2 is None):
        return R1

    vec_j = np.copy(vec2)
    if(len(vec2) == 4): vec_j = vec_j[1:4]
    vec_j = np.dot(R1,vec_j)

    # Now, we want to rotate vec_j as to get it into the x-y plane
    # (eliminating its z-component). To do this, we rotate about x,
    # since this will not move vec_i.
    theta = np.arctan(-vec_j[2]/vec_j[1])
    if(vec_j[1] > 0.): theta += np.pi
    R2 = GetRotationMatrix(0.,0.,theta)

    R = np.dot(R2,R1)
    return R

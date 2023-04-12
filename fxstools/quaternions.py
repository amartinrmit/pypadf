"""
quaternions.py

Functions to perform 3D rotations with quaternions

This script provides a few basic tools to multiply quaternions and perform
3D rotations of vectors.

This file can not be executed independently. It is imported as a module 
and contains the following functions:

    * quaternion multiply  - multiply two quaternions
    * rotate_vector  - rotate a 3D vector using quaternions
    * rotation_matrix - generate a 3D rotation matrix using quaternions
    * random_rotation_quaternion - generate a quaternion suitable for a
        random 3D rotation with uniform probability on SO3
    * quaternion_to_axis_angle - extract the rotation axis and rotation angle
        from a quaternion
    * random_vec_angle - generate a random rotation axis and angle
"""

import numpy as np



def quaternion_multiply(quaternion1, quaternion0):
    """multiply two quaternions using quaternion multiplication
    
    Parameters
    ----------
    quaternion1 - 1D numpy array with 4 floats
    quaternion0 - 1D numpy array with 4 floats
    
    Returns
    ----------
    1D numpy array with 4 floats to store result quaternion
    """
    w0, x0, y0, z0 = quaternion0
    w1, x1, y1, z1 = quaternion1
    return np.array([-x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
                     x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
                     -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
                     x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0], dtype=np.float64)


def rotate_vector( v, axis, angle ):
    """rotate a 3D vector using quaternions
    
    Parameters
    ----------
    v - 1D numpy array (float)
        3D vector to be rotated
    
    axis - 1D numpy array (float)
        unit vector specifying axis to rotate around
        
    angle - float
        angle of rotation in radians
    
    Returns
    ----------
    qout - 1D numpy array (float)
        the rotated vector
    """

    alen = np.sqrt(np.sum(axis*axis))
    axis *= 1.0/alen

    qvect = np.zeros(4)
    qvect[1:] = v

    q = np.zeros(4)
    q[0] = np.cos(angle/2)
    q[1:] = axis*np.sin(angle/2)

    qstar = np.copy(q)
    qstar[1:] *= -1

    qout = quaternion_multiply( qvect, qstar )
    qout = quaternion_multiply( q, qout )

    return qout[1:]


def rotation_matrix( axis, angle ):
    """generate a 3D rotation matrix using a quaternion formula
    
    Parameters
    ----------
    axis - 1D numpy array (float)
        unit vector specifying axis to rotate around
        
    angle - float
        angle of rotatation in radians
    
    Returns
    ----------
    R - 2D numpy array 3x3 (float)
        the rotated matrix
    """


    alen = np.sqrt(np.sum(axis*axis))
    axis *= 1.0/alen

    q = np.zeros(4)
    q[0] = np.cos(angle/2)
    q[1:] = axis*np.sin(angle/2)


    R = np.zeros( (3,3) )

    R[0,0] = 1 - 2*q[2]*q[2] - 2*q[3]*q[3]
    R[0,1] = 2*(q[1]*q[2] - q[3]*q[0])
    R[0,2] = 2*(q[1]*q[3] + q[2]*q[0])
    
    R[1,0] = 2*(q[1]*q[2] + q[3]*q[0])
    R[1,1] = 1 - 2*q[1]*q[1] - 2*q[3]*q[3]
    R[1,2] = 2*(q[2]*q[3] - q[1]*q[0])
    
    R[2,0] = 2*(q[1]*q[3] - q[2]*q[0])
    R[2,1] = 2*(q[2]*q[3] + q[1]*q[0])
    R[2,2] = 1 - 2*q[1]*q[1] - 2*q[2]*q[2]

    return R

def random_rotation_quaternion():
    """generate a quaternion suitable for a
        random 3D rotation with uniform probability on SO3
        
    Returns
    ----------
    q- 1D numpy array with 4 elements (float)
        quaternion of a random sample of SO3
    """
    
    r = np.random.rand(3)
    s = np.sqrt(1-r[0])
    s2 = np.sqrt(r[0])

    q = np.array([ s*np.sin(2*np.pi*r[1]), s*np.cos(2*np.pi*r[1]), s2*np.sin(2*np.pi*r[2]), s2*np.cos(2*np.pi*r[2])  ])
    return q

def quaternion_to_axis_angle(q):
    """generate a 3D rotation matrix using a quaternion formula
    
    Parameters
    ----------
    q - 1D numpy array, 4 elements (float)
        quaternion to convert
        
    angle - float
        angle of rotate
    
    Returns
    ----------
    v - 1D numpy array with 3 elements (float)
        unit vector of the rotation axis
        
    angle - float
        angle of rotation
    """
    
    angle = 2*np.arccos(q[0])
    v = q[1:] / np.sin( angle/2.0)
    return v, angle*180.0/np.pi

def random_vec_angle():
    """generates a random rotation axis and angle uniformly sampled on SO3
    
    Parameters
    ----------
    None
    
    Returns
    ----------
    v - 1D numpy array with 3 elements (float)
        unit vector of the rotation axis
        
    angle - float
        angle of rotation
    """
    q = random_rotation_quaternion()
    return quaternion_to_axis_angle(q)

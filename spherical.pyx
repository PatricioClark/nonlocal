# Spherical

# Cython module containing the necessary functions to calculate the nonlocal integrals
# The C part is located in ./lib

cimport decl
import numpy as np

# Dimensions of the work domain
Nx = decl.Nx
Ny = decl.Ny
Nz = decl.Nz

def PrintNxNyNz():
    decl.PrintNxNyNz()

def IntegrateInVolume(result, dr, R, alpha, xs, ys, zs):
    """
    Integrate the field set by set_field() using a 3D L-scheme on a given set
    of points. The arrays result, xs, ys, zs must all have the same dimensions.

    Parameters
    ----------
    result : ndarray
        Array into which the results will be placed.
    dr : float
        Step size for the discretization along the radial direction.
    R : float
        Radius of the nonlocal integration.
    alpha : float
        Fractional order. Must be between 0 and 1.
    xs, ys, zs: ndarrays.
        x, y and z coordinates in which to calculate the integral.

    Returns
    -------
    result : ndarray
        Value of the nonlocal integral at each point.
    """
    # Check dims then flatten arrays
    if not type(result)==np.ndarray:
        raise TypeError("""result, xs, ys, and zs should all be arrays
                           and have the same shape""")
    for arr in [xs, ys, zs]:
        if not np.shape(result)==np.shape(xs):
            raise TypeError("""result, xs, ys, and zs should all be arrays
                               and have the same shape""")
    N = len(result.ravel())

    # Create memoryview of the arrays
    cdef double[:] result_mv = result.ravel()
    cdef double[:] xs_mv = xs.ravel()
    cdef double[:] ys_mv = ys.ravel()
    cdef double[:] zs_mv = zs.ravel()

    # Calculate integral at each point and return array
    decl.IntegrateInVolume(&result_mv[0], dr, R, alpha,
                           &xs_mv[0], &ys_mv[0], &zs_mv[0], N)
    return result

def set_field(field):
    """
    Set the field used for the integration

    Parameters
    ----------
    field : ndarray
        Field to be integrated
    """
    cdef double[:,:,:] c_field_mv = decl.field
    cdef double[:,:,:] p_field_mv = field
    c_field_mv[:] = p_field_mv

def set_domain(dx,y_domain,dz,check_boundary=True):
    """
    Set the domain parameters. x and z are assumed to be equispaced.

    Parameters
    ----------
    dx : float
        Step size in x
    y_domain : ndarray
        y coordinates
    dz : float
        Step size in z
    """
    decl.dx = dx
    decl.dz = dz
    decl.y_domain = y_domain
    if check_boundary:
        decl.checkBoundary = 1
    else:
        decl.checkBoundary = 0

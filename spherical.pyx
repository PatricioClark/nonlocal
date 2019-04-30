# Cython wrappers

cimport decl
import numpy as np

# Dimensions of the work domain
Nx = decl.Nx
Ny = decl.Ny
Nz = decl.Nz

def IntegrateInVolume(result, dr, R, alpha, xs, ys, zs):
    # Check dims then flatten arrays
    if not type(result)==np.ndarray:
        raise TypeError("""result, xs, ys, and zs should all be arrays
                           and have the same shape""")
    for arr in [xs, ys, zs]:
        if not np.shape(result)==np.shape(xs):
            raise TypeError("""result, xs, ys, and zs should all be arrays
                               and have the same shape""")
    # result = result.ravel()
    # xs = xs.ravel()
    # ys = ys.ravel()
    # zs = zs.ravel()
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
    decl.field = field

def set_domain(dx,y_domain,dz):
    decl.dx = dx
    decl.dz = dz
    decl.y_domain = y_domain

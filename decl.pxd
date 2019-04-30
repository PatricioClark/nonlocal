# Declarations wrapper
# This step is performed so I can then reuse function names in the module

# Definitions
DEF _Nx=50
DEF _Ny=400
DEF _Nz=50

cdef enum:
    Nx = _Nx
    Ny = _Ny
    Nz = _Nz

# External declarations from C part
cdef extern from "spherical.h":
    # Global variables
    double field[_Nx][_Ny][_Nz]
    double dx, dz
    double ys[_Ny]

    # Functions
    double IntegrateInVolume(double dr, double R, double alpha,
                             double x0, double y0, double z0)

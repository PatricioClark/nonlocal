# Declarations wrapper
# This step is performed so I can then reuse function names in the module

# Definitions
DEF _Nx=480
DEF _Ny=295
DEF _Nz=960
cdef enum:
    Nx = _Nx
    Ny = _Ny
    Nz = _Nz

# External declarations from C part
cdef extern from "spherical.h":
    # Global variables
    int    checkBoundary
    double field[_Nx][_Ny][_Nz]
    double dx, dz
    double y_domain[_Ny]

    # Functions
    int  SendNx()
    int  SendNy()
    int  SendNz()
    void IntegrateInVolume(double *result, double dr, double R, double alpha,
                           double *xs, double *ys, double *zs, int Npoints)
    void PrintNxNyNz()

# Declarations wrapper
# This step is performed so I can then reuse function names in the module

# Definitions
DEF _Nx=2048
DEF _Ny=512
DEF _Nz=1536
cdef enum:
    Nx = _Nx
    Ny = _Ny
    Nz = _Nz

# External declarations from C part
cdef extern from "spherical.h":
    # Global variables
    int    checkBoundary
    float field[_Nx][_Ny][_Nz]
    double dx, dz
    double y_domain[_Ny]

    # Functions
    int  SendNx()
    int  SendNy()
    int  SendNz()
    void IntegrateInVolume(float *result, double dr, double R, double alpha,
                           float *xs, float *ys, float *zs, int Npoints)

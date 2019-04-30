# Cython wrappers

cimport decl
import numpy as np

Nx = decl.Nx
Ny = decl.Ny
Nz = decl.Nz

def IntegrateInVolume(dr, R, alpha, x0, y0, z0):
    return decl.IntegrateInVolume(dr, R, alpha, x0, y0, z0)

def set_field(field):
    decl.field = field

def set_dxdydz(dx,ys,dz):
    decl.dx = dx
    decl.ys = ys
    decl.dz = dz

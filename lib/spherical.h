#ifndef SPHERICAL_H
#define SPHERICAL_H

// Includes
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

// If these values are changed, the pyx file has to be updated accordingly
# define Nx 240
# define Ny 100
# define Nz 480

int    checkBoundary;
double field[Nx][Ny][Nz];
double dx, dz;
double y_domain[Ny];

double Trilinear(double x, double y, double z);
void   PointsInSphere(int N, double radius,
                      double x0, double y0, double z0,
                      double *xsphere, double *ysphere, double *zsphere);
double IntegrateInArea(int N, double radius,
                       double x0, double y0, double z0);
void IntegrateInVolume(double *result, double dr, double R, double alpha,
                       double *xs, double *ys, double *zs, int Npoints);
#endif

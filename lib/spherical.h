#ifndef SPHERICAL_H
#define SPHERICAL_H

// Includes
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

// If these values are changed, the pyx file has to be updated accordingly
# define Nx 50
# define Ny 400
# define Nz 50

double field[Nx][Ny][Nz];
double dx, dz;
double ys[Ny];

double Trilinear(double x, double y, double z);
void   PointsInSphere(int N, double radius,
                      double x0, double y0, double z0,
                      double xs[N], double ys[N], double zs[N]);
double IntegrateInArea(int N, double radius,
                         double x0, double y0, double z0);
double IntegrateInVolume(double dr, double R, double alpha,
                         double x0, double y0, double z0);
#endif

#ifndef SPHERICAL_H
#define SPHERICAL_H

// Includes
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

// If these values are changed, the pyx file has to be updated accordingly
# define Nx 1024
# define Ny 1024
# define Nz 1024

int    checkBoundary;
float  field[Nx][Ny][Nz];
double dx, dz;
double y_domain[Ny];

int SendNx();
int SendNy();
int SendNz();

double Trilinear(double x, double y, double z);
void   PointsInSphere(int N, double radius,
                      double x0, double y0, double z0,
                      double *xsphere, double *ysphere, double *zsphere);
double IntegrateInArea(int N, double radius,
                       double x0, double y0, double z0);
void IntegrateInVolume(float *result, double dr, double R, double alpha,
                       float *xs, float *ys, float *zs, int Npoints);
#endif

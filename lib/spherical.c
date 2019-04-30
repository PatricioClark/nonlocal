/* Functions used to calculate flow and interpolations */

// Include file
#include "spherical.h"

double TrilinearInterpolation(double field[Nx][Ny][Nz], double x, double y, double z) {
  // Interpolate in 3D using trilinear method
  int    i;
  int    i0, i1, j0, j1, k0, k1;
  double x0, xd, y0, yd, z0, zd;

  double c00, c01, c10, c11;
  double c0,  c1;
  double c;

  // Check boundary conditions
  if ((y<-1.0) || (y>1.0)) 
    return 0.0;

  // Get x indexes (equispaced grid)
  i0 = floor(x/dx);
  x0 = i0*dx;
  i1 = i0 + 1;
  xd = (x-x0)/dx;

  // Get y indexes
  for (i=0; i<Ny; i++) 
    if (y<ys[i]) break;
  j0 = i-1;
  j1 = j0+1;
  y0 = ys[j0];
  yd = (y-y0)/(ys[j1]-y0);

  // Get z indexes (equispaced grid)
  k0 = floor(z/dz);
  z0 = k0*dz;
  k1 = k0 + 1;
  zd = (z-z0)/dz;

  // Interpolation in x
  c00 = field[i0][j0][k0]*(1-xd) + field[i1][j0][k0]*xd;
  c01 = field[i0][j0][k1]*(1-xd) + field[i1][j0][k1]*xd;
  c10 = field[i0][j1][k0]*(1-xd) + field[i1][j1][k0]*xd;
  c11 = field[i0][j1][k1]*(1-xd) + field[i1][j1][k1]*xd;

  // Interpolation in y
  c0 = c00*(1-yd) + c10*yd;
  c1 = c01*(1-yd) + c11*yd;

  // Interpolation in z
  c = c0*(1-zd) + c1*zd;

  return c;
}

void PointsInSphere(int N, double radius,
                    double x0, double y0, double z0,
                    double xs[N], double ys[N], double zs[N]) {
  // Taken from that paper
  int k;
  double hk, cte;
  double theta, phi, phi_prev;

  cte = 3.6/sqrt(N);
  phi_prev = 0.0;
  for (k=0; k<N; k++) {
    hk = -1.0 + (2.0*k)/(N-1.0);

    theta = acos(hk);

    if ((k==0) || (k==N-1)) {
      phi = 0.0;
    }
    else {
      phi = phi_prev + cte/sqrt(1-hk*hk);
      phi = fmod(phi, 2*M_PI);
      phi_prev = phi;
    }

    xs[k] = radius*sin(theta)*cos(phi) + x0;
    ys[k] = radius*sin(theta)*sin(phi) + y0;
    zs[k] = radius*cos(theta)          + z0;
  }
}

double IntegrateInArea(int N, double radius,
                         double x0, double y0, double z0) {
  // Integrate the strain in a sphere of a given radius centred around
  // (x0,y0,z0)

  int k;
  double value;
  double xs[N], ys[N], zs[N];

  // Get points
  PointsInSphere(N, radius, x0, y0, z0, xs, ys, zs);

  // Integrate
  value = 0.0;
  for (k=0; k<N; k++) {
    value += TrilinearInterpolation(field, xs[k], ys[k], zs[k]);
  }

  /* fprintf(stdout, "area integration: %f\n", (4*M_PI/N)*value); */
  return (4*M_PI/N)*value;
}

double IntegrateInVolume(double dr, double R,  double alpha,
                         double x0, double y0, double z0) {
  // Integrate the strain non-locally

  int i, Nr, Nsphere;
  double value, aux, radius;

  // Inner point
  value = 4*M_PI*TrilinearInterpolation(field, x0, y0, z0);

  // Spheres
  Nr = R/dr;
  for (i=1; i<Nr; i++) {
    radius  = dr*i;
    Nsphere = 4*M_PI*radius*radius/(dr*dr);
    aux = pow(i+1, 1-alpha) - pow(i, 1-alpha);
    value += aux*IntegrateInArea(Nsphere, radius, x0, y0, z0);
  }

  // Return value
  aux = 4*M_PI*tgamma(2-alpha);
  aux = pow(dr, 1-alpha)/aux;
  return aux*value;
}

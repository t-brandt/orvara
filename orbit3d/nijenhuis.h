
// A solver for Kepler's equation based on:
//
// Nijenhuis (1991)
// http://adsabs.harvard.edu/abs/1991CeMDA..51..319N
//
// and
//
// Markley (1995)
// http://adsabs.harvard.edu/abs/1995CeMDA..63..101M

#include "math.h"
#include <stdbool.h>

void sin_cos_reduc(double x, double *SnReduc, double *CsReduc);

// Calculates x - sin(x) and 1 - cos(x) to 20 significant digits for x in [0, pi)

void sin_cos_reduc(double x, double *SnReduc, double *CsReduc) {

  const double s[] = {1./6, 1./20, 1./42, 1./72, 1./110, 1./156, 1./210, 1./272, 1./342, 1./420};
  const double c[] = {0.5, 1./12, 1./30, 1./56, 1./90, 1./132, 1./182, 1./240, 1./306, 1./380};

  bool bigg = x > M_PI_2;
  double u = (bigg) ? M_PI - x : x;
  bool big = u > M_PI_2;
  double v = (big) ? M_PI_2 - u : u;
  double w = v * v;
 
  double ss = 1;
  double cc = 1;
  int i;  

  for (i = 9; i >= 1; --i) {
    ss = 1 - w * s[i] * ss;
    cc = 1 - w * c[i] * cc;
  }
  ss *= v * w * s[0];
  cc *= w * c[0];
  
  if (big) {
    *SnReduc = u - 1 + cc;
    *CsReduc = 1 - M_PI_2 + u + ss;
  } else {
    *SnReduc = ss;
    *CsReduc = cc;
  }
  if (bigg) {
    *SnReduc = 2 * x - M_PI + *SnReduc;
    *CsReduc = 2 - *CsReduc;
  }
}

//double solve_kepler(double M, double ecc, double *EA, double *sinE, double *cosE);
double solve_kepler(double M, double ecc);

//double solve_kepler(double M, double ecc, double *EA, double *sinE, double *cosE) {
double solve_kepler(double M, double ecc) {

  const double two_pi = 2 * M_PI;
  
  //double M_ref = two_pi * floor(M / two_pi);
  //M -= M_ref;

  bool high = M > M_PI;
  if (high) {
    M = two_pi - M;
  }

  double ome = 1.0 - ecc;
  
  // Get starter
  double M2 = M*M;
  double M3 = M2*M;
  double alpha = (3*M_PI + 1.6*(M_PI-fabs(M))/(1+ecc) )/(M_PI - 6/M_PI);
  double d = 3*ome + alpha*ecc;
  double r = 3*alpha*d*(d-ome)*M + M3;
  double q = 2*alpha*d*ome - M2;
  double q2 = q*q;
  double w = pow(fabs(r) + sqrt(q2*q + r*r), 2.0/3);
  double E = (2*r*w/(w*w + w*q + q2) + M) / d;

  // Approximate Mstar = E - e*sin(E) with numerically stability
  double sE, cE;
  sin_cos_reduc (E, &sE, &cE);
  
  // Refine the starter
  double f_0 = ecc * sE + E * ome - M;
  double f_1 = ecc * cE + ome;
  double f_2 = ecc * (E - sE);
  double f_3 = 1-f_1;
  double d_3 = -f_0/(f_1 - 0.5*f_0*f_2/f_1);
  double d_4 = -f_0/(f_1 + 0.5*d_3*f_2 + (d_3*d_3)*f_3/6);
  double d_42 = d_4*d_4;
  E -= f_0/(f_1 + 0.5*d_4*f_2 + d_4*d_4*f_3/6 - d_42*d_4*f_2/24);

  //sin_cos_reduc (E, &sE, &cE);
  //cE = 1 - cE;
  //sE = E - sE;
  if (high) {
    E = two_pi - E;
    //sE = -sE;
  }
  
  //*EA = E;
  //*sinE = sE;
  //*cosE = cE;
  
  return E; // + M_ref;
}

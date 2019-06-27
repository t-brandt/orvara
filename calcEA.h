/* 
 * Return the Taylor expansion for sine.  Accurate essentially to
 *  machine precision (last bit might be wrong) up to an argument of
 *  pi/4.  Faster than a straight-up call to sine.
 */

inline double shortsin(double x);

inline double shortsin(double x) {

  const double if3 = 1./6;
  const double if5 = 1./(6.*20);
  const double if7 = 1./(6.*20*42);
  const double if9 = 1./(6.*20*42*72);
  const double if11 = 1./(6.*20*42*72*110);
  const double if13 = 1./(6.*20*42*72*110*156);
  const double if15 = 1./(6.*20*42*72*110*156*210);
  
  double x2 = x*x;

  return x*(1 - x2*(if3 - x2*(if5 - x2*(if7 - x2*(if9 - x2*(if11 - x2*(if13 - x2*if15)))))));
}

/* 
 * Implements the series expansions for the eccentric anomaly starting
 * point as given in Raposo-Pulido & Pelaez (2017)
 */

inline double EAstart(double M, double ecc);

inline double EAstart(double M, double ecc) {

  double ome = 1. - ecc;
  double sqrt_ome = sqrt(ome);
    
  if (M < 1e-3*sqrt_ome*ome) {
    double xi = M/(ome*ome);
    double xi2_ome = xi*xi*ome;
    double eta = xi*(1 - xi2_ome*(ecc/6. - xi2_ome*(19*ecc - 9)/120. - xi2_ome/18.));
    return ome*eta;
  } else {
    double chi = M/(sqrt_ome*ome);
    double Lam = sqrt(8 + 9*chi*chi);
    double S = pow(Lam + 3*chi, 1./3);
    double sigma = 6*chi/(2 + S*S + 4./(S*S));
    double s2 = sigma*sigma;
    double denom = s2 + 2;
    double E = sigma*(1 + s2*ome*((s2 + 20)/(60.*denom) + s2*ome*(s2*s2*s2 + 25*s2*s2 + 340*s2 + 840)/(1400*denom*denom*denom)));
    return E*sqrt_ome;
  }
}

/* 
 * Implements the series expansions for the eccentric anomaly starting
 * point as given in Raposo-Pulido & Pelaez (2017)
 */

inline double MAmod(double M);
  
inline double MAmod(double M) {

  const double twopi = 2*3.14159265358979323846264338327950288;
    
  if (M < twopi && M >= 0)
    return M;
  
  if (M > twopi) {
    M -= twopi;
    if (M > twopi) 
      return fmod(M, twopi);
    else 
      return M;
  } else {
    M += twopi;
    if (M < 0)
      return fmod(M, twopi) + twopi;
    else
      return M;
  }
}

/* 
 * Implements the piecewise quintic polynomial given in Raposo-Pulido
 * & Pelaez (2017): calculates the coefficients.  Polynomials are
 * defined on the ranges bounds[i], bounds[i + 1] and take as their
 * argument x = M - bounds[i], where M is the mean anomaly in radians,
 * assumed to be between 0 and pi.
 */

void getbounds(double bounds[], double EA_tab[], double ecc);

void getbounds(double bounds[], double EA_tab[], double ecc) {
  
  const double pi = 3.14159265358979323846264338327950288;
  const double pi_d_12 = 3.14159265358979323846264338327950288/12;
  const double pi_d_6 = 3.14159265358979323846264338327950288/6;
  const double pi_d_4 = 3.14159265358979323846264338327950288/4;
  const double pi_d_3 = 3.14159265358979323846264338327950288/3;
  const double fivepi_d_12 = 3.14159265358979323846264338327950288*5./12;
  const double pi_d_2 = 3.14159265358979323846264338327950288/2;
  const double sevenpi_d_12 = 3.14159265358979323846264338327950288*7./12;
  const double twopi_d_3 = 3.14159265358979323846264338327950288*2./3;
  const double threepi_d_4 = 3.14159265358979323846264338327950288*3./4;
  const double fivepi_d_6 = 3.14159265358979323846264338327950288*5./6;
  const double elevenpi_d_12 = 3.14159265358979323846264338327950288*11./12;
  
  double g2s_e = 0.2588190451025207623489*ecc;
  double g3s_e = 0.5*ecc;
  double g4s_e = 0.7071067811865475244008*ecc;
  double g5s_e = 0.8660254037844386467637*ecc;
  double g6s_e = 0.9659258262890682867497*ecc;
  double g2c_e = g6s_e;
  double g3c_e = g5s_e;
  double g4c_e = g4s_e;
  double g5c_e = g3s_e;
  double g6c_e = g2s_e;
  
  bounds[0] = 0;
  bounds[1] = pi_d_12 - g2s_e;
  bounds[2] = pi_d_6 - g3s_e;
  bounds[3] = pi_d_4 - g4s_e;
  bounds[4] = pi_d_3 - g5s_e;
  bounds[5] = fivepi_d_12 - g6s_e;
  bounds[6] = pi_d_2 - ecc;
  bounds[7] = sevenpi_d_12 - g6s_e;
  bounds[8] = twopi_d_3 - g5s_e;
  bounds[9] = threepi_d_4 - g4s_e;
  bounds[10] = fivepi_d_6 - g3s_e;
  bounds[11] = elevenpi_d_12 - g2s_e;
  bounds[12] = pi;
  
  EA_tab[1] = 1/(1. - ecc);
  EA_tab[7] = 1/(1. - g2c_e);
  EA_tab[13] = 1/(1. - g3c_e);
  EA_tab[19] = 1/(1. - g4c_e);
  EA_tab[25] = 1/(1. - g5c_e);
  EA_tab[31] = 1/(1. - g6c_e);
  EA_tab[37] = 1;
  EA_tab[43] = 1/(1. + g6c_e);
  EA_tab[49] = 1/(1. + g5c_e);
  EA_tab[55] = 1/(1. + g4c_e);
  EA_tab[61] = 1/(1. + g3c_e);
  EA_tab[67] = 1/(1. + g2c_e);
  EA_tab[73] = 1/(1. + ecc);
  
  EA_tab[2] = 0;
  EA_tab[8] = -0.5*g2s_e/((1. - g2c_e)*(1. - g2c_e)*(1. - g2c_e));
  EA_tab[14] = -0.5*g3s_e/((1. - g3c_e)*(1. - g3c_e)*(1. - g3c_e));
  EA_tab[20] = -0.5*g4s_e/((1. - g4c_e)*(1. - g4c_e)*(1. - g4c_e));
  EA_tab[26] = -0.5*g5s_e/((1. - g5c_e)*(1. - g5c_e)*(1. - g5c_e));
  EA_tab[32] = -0.5*g6s_e/((1. - g6c_e)*(1. - g6c_e)*(1. - g6c_e));
  EA_tab[38] = -0.5*ecc;
  EA_tab[44] = -0.5*g6s_e/((1. + g6c_e)*(1. + g6c_e)*(1. + g6c_e));
  EA_tab[50] = -0.5*g5s_e/((1. + g5c_e)*(1. + g5c_e)*(1. + g5c_e));
  EA_tab[56] = -0.5*g4s_e/((1. + g4c_e)*(1. + g4c_e)*(1. + g4c_e));
  EA_tab[62] = -0.5*g3s_e/((1. + g3c_e)*(1. + g3c_e)*(1. + g3c_e));
  EA_tab[68] = -0.5*g2s_e/((1. + g2c_e)*(1. + g2c_e)*(1. + g2c_e));
  EA_tab[74] = 0;
  
  double B0, B1, B2, dx;
  int i, k;
  for (i = 0; i < 12; i++) {
    dx = bounds[i + 1] - bounds[i];
    k = 6*i;
    EA_tab[k] = i*pi_d_12;
    B0 = (pi_d_12 - EA_tab[k + 1]*dx - EA_tab[k + 2]*dx*dx)/(dx*dx*dx);
    B1 = (EA_tab[k + 7] - EA_tab[k + 1] - 2*EA_tab[k + 2]*dx)/(3*dx*dx);
    B2 = (EA_tab[k + 8] - EA_tab[k + 2])/(3*dx);
    EA_tab[k + 3] = 3*B2 - 12*B1 + 10*B0;
    EA_tab[k + 4] = (-6*B2 + 21*B1 - 15*B0)/dx;
    EA_tab[k + 5] = (3*B2 - 9*B1 + 6*B0)/(dx*dx);
  }
  
  return;
}

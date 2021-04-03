
void kepler_goatherd(int N_it, int n, double M[], double e, double E[]);

void kepler_goatherd(int N_it, int n, double M[], double e, double E[])
{
  
    // Solve Kepler's equation via the contour integration method of Philcox et al. (2021)
    // This uses techniques described in Ullisch (2020) to solve the `geometric goat problem'.
    // N_it specifies the number of grid-points.

  double ft_gx2, ft_gx1, this_ell, freq, zR, zI, cosC, sinC, esinRadius, ecosRadius, center;
  double fxR, fxI, ftmp, tmpcosh, tmpsinh, tmpcos, tmpsin;

  // Define sampling points (actually use one more than this)
  int N_points = N_it - 2;
  int N_fft = (N_it - 1)*2;
  
  // Define contour radius
  double radius = e/2.;
  //double M_PI = 3.141592653589793;
  //
  // Generate e^{ikx} sampling points and precompute real and imaginary parts
  double cf, sf;
  double *exp2R = malloc(sizeof(double)*N_points);
  double *exp2I = malloc(sizeof(double)*N_points);
  double *exp4R = malloc(sizeof(double)*N_points);
  double *exp4I = malloc(sizeof(double)*N_points);
  double *coshI = malloc(sizeof(double)*N_points);
  double *sinhI = malloc(sizeof(double)*N_points);
  double *ecosR = malloc(sizeof(double)*N_points);
  double *esinR = malloc(sizeof(double)*N_points);

  if (exp2R == NULL || exp2I == NULL || exp4R == NULL || exp4I == NULL ||
      coshI == NULL || sinhI == NULL || ecosR == NULL || esinR == NULL) {
    printf("Cannot allocate memory.\n");
    exit(-1);
  }

  int jj, i;
    
  for(jj = 0; jj < N_points; jj++){
    // NB: j = jj+1
      freq = 2.0*M_PI*(jj+1)/N_fft;
      cf = cos(freq);
      sf = sin(freq);
      exp2R[jj] = cf;
      exp2I[jj] = sf;
      exp4R[jj] = cf*cf-sf*sf;
      exp4I[jj] = 2.0*cf*sf;
      coshI[jj] = cosh(radius*exp2I[jj]);
      sinhI[jj] = sinh(radius*exp2I[jj]);
      ecosR[jj] = e*cos(radius*exp2R[jj]);
      esinR[jj] = e*sin(radius*exp2R[jj]);
  }
  
  // Precompute e sin(e/2) and e cos(e/2)
  esinRadius = e*sin(radius);
  ecosRadius = e*cos(radius);
  
  // Iterate over array of mean anomalies
  for(i = 0; i < n; i++){
    this_ell = M[i];
    
    // Define contour center for each ell and precompute sin(center), cos(center)
      if(this_ell<M_PI) center = this_ell+e/2;
      else center = this_ell-e/2;
      sinC = sin(center);
      cosC = cos(center);
      E[i] = center;
      //if (i == 0) printf("%.5g\n", center);
      // Accumulate Fourier coefficients
      // NB: we halve the range by symmetry, absorbing factor of 2 into ratio

      ///////////////
      // Separate out j = 0 piece, which is simpler

      // Compute z in real and imaginary parts (zI = 0 here)
      zR = center + radius;

      // Compute e*sin(zR) from precomputed quantities
      tmpsin = sinC*ecosRadius+cosC*esinRadius; // sin(zR)

      // Compute f(z(x)) in real and imaginary parts (fxI = 0)
      fxR = zR - tmpsin - this_ell;

      // Add to array, with factor of 1/2 since an edge
      ft_gx2 = 0.5/fxR;
      ft_gx1 = 0.5/fxR;

      ///////////////
      // Compute for j = 1 to N_points
      // NB: j = jj+1
      for(jj=0;jj<N_points;jj++){

        // Compute z in real and imaginary parts
        zR = center + radius*exp2R[jj];
        zI = radius*exp2I[jj];

        // Compute f(z(x)) in real and imaginary parts
        // can use precomputed cosh / sinh / cos / sin for this!
        tmpcosh = coshI[jj]; // cosh(zI)
        tmpsinh = sinhI[jj]; // sinh(zI)
        tmpsin = sinC*ecosR[jj]+cosC*esinR[jj]; // e sin(zR)
        tmpcos = cosC*ecosR[jj]-sinC*esinR[jj]; // e cos(zR)

        fxR = zR - tmpsin*tmpcosh-this_ell;
        fxI = zI - tmpcos*tmpsinh;

        // Compute 1/f(z) and append to array
        ftmp = fxR*fxR+fxI*fxI;
        fxR = fxR/ftmp;
        fxI = fxI/ftmp;

        ft_gx2 = ft_gx2 + (exp4R[jj]*fxR+exp4I[jj]*fxI);
        ft_gx1 = ft_gx1 + (exp2R[jj]*fxR+exp2I[jj]*fxI);
      }

      ///////////////
      // Separate out j = N_it piece, which is simpler

      // Compute z in real and imaginary parts (zI = 0 here)
      zR = center - radius;

      // Compute sin(zR) from precomputed quantities
      tmpsin = sinC*ecosRadius-cosC*esinRadius; // sin(zR)

      // Compute f(z(x)) in real and imaginary parts (fxI = 0 here)
      fxR = zR - tmpsin-this_ell;

      // Add to sum, with 1/2 factor for edges
      ft_gx2 = ft_gx2 + 0.5/fxR;
      ft_gx1 = ft_gx1 - 0.5/fxR;

      ///////////////
      // Compute E(ell)
    E[i] = E[i] + radius*ft_gx2/ft_gx1;
  }

  free(exp2R);
  free(exp2I);
  free(exp4R);
  free(exp4I);
  free(coshI);
  free(sinhI);
  free(ecosR);
  free(esinR);
  
  return;

}

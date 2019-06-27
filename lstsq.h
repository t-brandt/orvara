void lstsq_C(double A_in[], double b[], int m, int n, double coef[]);

void lstsq_C(double A_in[], double b[], int m, int n, double coef[]) {

  int flag, its, jj, j, i, l, k, nm, inc;
  double c, f, h, s, x, y, z, tmp, sw, eps, tsh;
  double anorm, g, scale;

  inc = 1;
  eps = 2.3e-16;
  /*
    ###############################################################
    # None of these arrays will be visible outside of this routine.
    # They are used to construct the SVD.
    ###############################################################
  */

  double *tmparr = malloc(n*sizeof(double));
  double *su = malloc(m*sizeof(double));
  double *sv = malloc(n*sizeof(double));
  double *w = malloc(n*sizeof(double));
  double *rv1 = malloc(n*sizeof(double));
  double *v = malloc(n*n*sizeof(double));
  double *A = malloc(n*m*sizeof(double));
  if (A == NULL || v == NULL || rv1 == NULL || w == NULL || sv == NULL || su == NULL || tmparr == NULL) {
    printf("Cannot allocate memory.\n");
    exit(-1);
  }
  
  for (i = 0; i < m; i++) 
    for (j = 0; j < n; j++) 
      A[i*m + j] = A_in[i*m + j];

  /*
  ###############################################################
  # The code below is largely copied from Numerical Recipes and from
  # http://www.public.iastate.edu/~dicook/JSS/paper/code/svd.c
  ###############################################################
  */

  scale = 0.;
  g = 0.;
  anorm = 0.;
  for (i = 0; i < n; i++) {
    l = i + 1;
    rv1[i] = scale*g;
    g = 0.;
    s = 0.;
    scale = 0.;
    if (i < m) {
      for (k = i; k < m; k++)
	scale = scale + fabs(A[k*m + i]);
      if (scale != 0) {
	for (k = i; k < m; k++) {
	  A[k*m + i] = A[k*m + i]/scale;
	  s = s + A[k*m + i]*A[k*m + i];
	}

	f = A[i*m + i];
	g = -1*sqrt(s);
	if (f < 0)
	  g = -1*g;
	h = f*g - s;
	A[i*m + i] = f - g;
	if (i != n - 1) {
	  for (j = l; j < n; j++) {
	    s = 0;
	    for (k = i; k < m; k++)
	      s = s + A[k*m + i]*A[k*m + j];
	    f = s/h;
	    for (k = i; k < m; k++)
	      A[k*m + j] = A[k*m + j] + f*A[k*m + i];
	  }
	}
	for (k = i; k < m; k++) 
	  A[k*m + i] = A[k*m + i]*scale;
      }
    }
    w[i] = scale*g;
    g = 0.;
    s = 0.;
    scale = 0.;
    if (i < m && i != n - 1) {
      for (k = l; k < n; k++)
	scale = scale + fabs(A[i*m + k]);
      if (scale != 0) {
	for (k = l; k < n; k++) {
	  A[i*m + k] = A[i*m + k]/scale;
	  s = s + A[i*m + k]*A[i*m + k];
	}
	f = A[i*m + l];
	g = -1*sqrt(s);
	if (f < 0) 
	  g = -1*g;
	h = f*g - s;
	A[i*m + l] = f - g;
	for (k = l; k < n; k++) 
	  rv1[k] = A[i*m + k]/h;
	if (i != m - 1) {
	  for (j = l; j < m; j++) {
	    s = 0;
	    for (k = l; k < n; k++)
		s = s + A[j*m + k]*A[i*m + k];
	    for (k = l; k < n; k++)
		A[j*m + k] = A[j*m + k] + s*rv1[k];
	  }
	}
	for (k = l; k < n; k++) 
	  A[i*m + k] = A[i*m + k]*scale;
      }
    }
    if (fabs(w[i]) + fabs(rv1[i]) > anorm)
      anorm = fabs(w[i]) + fabs(rv1[i]);
  }
  for (i = n - 1; i > -1; i--) {
    if (i < n - 1) {
      if (g != 0) {
	for (j = l; j < n; j++)
	  v[j*n + i] = A[i*m + j]/A[i*m + l]/g;
	for (j = l; j < n; j++) {
	  s = 0;
	  for (k = l; k < n; k++)
	    s = s + A[i*m + k]*v[k*n + j];
	  for (k = l; k < n; k++)
	    v[k*n + j] = v[k*n + j] + s*v[k*n + i];
	}
      }
      for (j = l; j < n; j++) {
	v[i*n + j] = 0.;
	v[j*n + i] = 0.;
      }
    }
    v[i*n + i] = 1.;
    g = rv1[i];
    l = i;
  }
  
  for (i = n - 1; i > -1; i--) {
    l = i + 1;
    g = w[i];
    if (i < n - 1) 
      for (j = l; j < n; j++) 
	A[i*m + j] = 0.;
    if (g != 0) {
      g = 1./g;
      if (i != n - 1) {
	for (j = l; j < n; j++) {
	  s = 0;
	  for (k = l; k < m; k++) 
	    s = s + A[k*m + i]*A[k*m + j];
	  f = (s/A[i*m + i])*g;
	  for (k = i; k < m; k++) 
	    A[k*m + j] = A[k*m + j] + f*A[k*m + i];
	}
      }
      for (j = i; j < m; j++) 
	A[j*m + i] = A[j*m + i]*g;
    } else {
      for (j = i; j < m; j++) 
	A[j*m + i] = 0.;
    }
    A[i*m + i] = A[i*m + i] + 1.;
  }
  for (k = n - 1; k > -1; k--) {
    for (its = 0; its < 30; its++) {
      flag = 1;
      for (l = k; l > -1; l--) {
	nm = l - 1;
	if (fabs(rv1[l]) + anorm == anorm) {
	  flag = 0;
	  break;
	}
	if (fabs(w[nm]) + anorm == anorm) 
	  break;
      }
      if (flag != 0) {
	c = 0.;
	s = 1.;
	for (i = l; i < k + 1; i++) {
	  f = s*rv1[i];
	  if (fabs(f) + anorm != anorm) {
	    g = fabs(w[i]);
	    h = sqrt(f*f + g*g);
	    w[i] = h;
	    h = 1./h;
	    c = g*h;
	    s = -1.*f*h;
	    for (j = 0; j < m; j++) {
	      y = A[j*m + nm];
	      z = A[j*m + i];
	      A[j*m + nm] = y*c + z*s;
	      A[j*m + i] = z*c - y*s;
	    }
	  }
	}
      }
      z = w[k];
      if (l == k) {
	if (z < 0.) {
	  w[k] = -1.*z;
	  for (j = 0; j < n; j++) 
	    v[j*n + k] = -1.*v[j*n + k];
	}
	break;
      }
      
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z)*(y + z) + (g - h)*(g + h))/(2.*h*y);

      g = sqrt(1. + f*f);
      tmp = g;
      if (f < 0) 
	tmp = -1*tmp;
            
      f = ((x - z)*(x + z) + h*((y/(f + tmp)) - h))/x;
            
      c = 1.;
      s = 1.;
      for (j = l; j < nm + 1; j++) {
	i = j + 1;
	g = rv1[i];
	y = w[i];
	h = s*g;
	g = c*g;
		
	z = sqrt(f*f + h*h);

	rv1[j] = z;
	c = f/z;
	s = h/z;
	f = x*c + g*s;
	g = g*c - x*s;
	h = y*s;
	y = y*c;
	for (jj = 0; jj < n; jj++) {
	  x = v[jj*n + j];
	  z = v[jj*n + i];
	  v[jj*n + j] = x*c + z*s;
	  v[jj*n + i] = z*c - x*s;
	}	    
	z = sqrt(f*f + h*h);

	w[j] = z;
	if (z != 0) {
	  z = 1./z;
	  c = f*z;
	  s = h*z;
	}
	f = c*g + s*y;
	x = c*y - s*g;
	for (jj = 0; jj < m; jj++) {
	  y = A[jj*m + j];
	  z = A[jj*m + i];
	  A[jj*m + j] = y*c + z*s;
	  A[jj*m + i] = z*c - y*s;
	}
      }
      rv1[l] = 0.;
      rv1[k] = f;
      w[k] = x;
    }
  }
  inc = 1;
  while (1) {
    inc = inc*3 + 1;
    if (inc > n)
      break;
  }
  while (1) {
    inc = inc/3;
    for (i = inc; i < n; i++) {
      sw = w[i];
      for (k = 0; k < m; k++)
	su[k] = A[k*m + i];
      for (k = 0; k < n; k++)
	sv[k] = v[k*n + i];
      j = i;

      while (w[j - inc] < sw) {
	w[j] = w[j - inc];
	for (k = 0; k < m; k++) 
	  A[k*m + j] = A[k*m + j - inc];
	for (k = 0; k < n; k++)
	  v[k*n + j] = v[k*n + j - inc];
	j = j - inc;
	if (j < inc)
	  break;
      }
      w[j] = sw;
      for (k = 0; k < m; k++)
	A[k*m + j] = su[k];
      for (k = 0; k < n; k++)
	v[k*n + j] = sv[k];
    }
    if (inc <= 1)
      break;
  }

  for (k = 0; k < n; k++) {
    jj = 0;
    for (i = 0; i < m; i++) 
      if (A[i*m + k] < 0)
	jj = jj + 1;
    for (j = 0; j < n; j++)
      if (v[j*n + k] < 0)
	jj = jj + 1;
    if (jj > (m + n)/2) {
      for (i = 0; i < m; i++)
	A[i*m + k] = -1.*A[i*m + k];
      for (j = 0; j < n; j++)
	v[j*n + k] = -1.*v[j*n + k];
    }
  }

  tsh = 0.5*sqrt(m + n + 1.)*w[0]*eps;
  
  for (j = 0; j < n; j++) {
    s = 0.;
    if (w[j] > tsh) {
      for (i = 0; i < m; i++)
	s = s + A[i*m + j]*b[i];
      s = s/w[j];
    }
    tmparr[j] = s*1.;
  }
  
  for (j = 0; j < n; j++) {
    s = 0.;
    for (jj = 0; jj < n; jj++) 
      s = s + v[j*n + jj]*tmparr[jj];
    coef[j] = s*1.;
  }
  
  free(tmparr);
  free(su);
  free(sv);
  free(w);
  free(rv1);
  free(v);
  free(A);
  return;
}

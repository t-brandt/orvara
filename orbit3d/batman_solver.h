double kepler_batman(double M, double e);

double kepler_batman(double M, double e)	//calculates the eccentric anomaly (see Seager Exoplanets book:  Murray & Correia eqn. 5 -- see section 3)
{
	double E = M, eps = 1.0e-7;
	double fe, fs;
  const int max_iter = 30;
  int i = 0;

	// modification from LK 05/07/2017:
	// add fmod to ensure convergence for diabolical inputs (following Eastman et al. 2013; Section 3.1)
	while(fmod(fabs(E - e*sin(E) - M), 2.*M_PI) > eps && i < max_iter)
	{
		fe = fmod(E - e*sin(E) - M, 2.*M_PI);
		fs = fmod(1 - e*cos(E), 2.*M_PI);
		E = E - fe/fs;
    ++i;
	}
	return E;
}

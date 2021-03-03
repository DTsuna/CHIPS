#include "flux.h"
#include "saha.h"
#include "opacity.h"
#include "constant.h"
#include <math.h>

//i must be larger than 1!!

double calc_flux_i(double r_ini, double r[], double E[], double U[], double rho[], double dt, int i, const int nsize)
{
	double Ec, rhoc, Tc;
	double T1, T2, mu;
	double dEdr;
	double kappa, lambda, R;

	saha(rho[i-1], U[2*(i-1)+1], &mu, &T1);
	saha(rho[i], U[2*(i)+1], &mu, &T2);
	Tc = (T1+T2)/2.;
	Ec = (E[2*(i-1)+1]+E[2*i+1])/2.;
	rhoc = (rho[i-1]+rho[i])/2.;
	kappa = kappa_r(rhoc, Tc);

	dEdr = (E[2*i+1]-E[2*(i-1)+1])/(r[i]-r[i-1]);
	R = fabs(dEdr)/(kappa*rhoc*Ec);
	lambda = (2.+R)/(6.+3.*R+R*R);

	return -(P_C)/(kappa*rhoc)*lambda*dEdr;
}

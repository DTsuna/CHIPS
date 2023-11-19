#include <math.h>
#include "constant.h"
#include "function.h"
#include "opacity.h"
#include "itgadve.h"
#include "saha.h"
#include "trid.h"

double func_lambda(double R)
{
	return (2.+R)/(6.+3.*R+R*R);
}

void matrix_E(double r[], double E[], double U[], double rho[], double dt, double a[], double b[], double c[], const int nsize)
{
	double lambda[nsize], R;
	double dr3[nsize], drm[nsize]; //dr3 represents r[i+1]^3-r[i]^3;
	double T_g[nsize], mu[nsize], kappa[nsize];
	double rhom[nsize];
	double Xfac[nsize];
	double k_l, k_r;
	double kappa_mid[nsize];
	int i;

	for(i = 0; i < nsize; i++){
		saha(rho[i], U[2*i+1], mu+i, T_g+i);
		rhom[i] = rho_csm(r[i]);
		kappa_mid[i] = kappa_r(rho_csm((r[i]+r[i+1])/2.), T_g[i]);
	}

	for(i = 1; i < nsize; i++){
//		kappa[i] = kappa_r(rhom[i], (T_g[i-1]+T_g[i])/2.0);
		k_l = kappa_mid[i-1];
		k_r = kappa_mid[i];
		kappa[i] = (E[2*i+1]+E[2*(i-1)+1])/(E[2*i+1]/k_r+E[2*(i-1)+1]/k_l);
		drm[i] = (r[i+1]-r[i-1])/2.;
	}

	i = 0;
	dr3[i] = (r[i+1]-r[i])*(r[i+1]*r[i+1]+r[i+1]*r[i]+r[i]*r[i]);

	for(i = 1; i < nsize; i++){
		dr3[i] = (r[i+1]-r[i])*(r[i+1]*r[i+1]+r[i+1]*r[i]+r[i]*r[i]);
		R = fabs((E[2*i+1]-E[2*(i-1)+1])/drm[i]/(kappa[i]*rhom[i]*(E[2*i+1]+E[2*(i-1)+1])/2.));
		lambda[i] = func_lambda(R);
		Xfac[i] = -(P_C)*lambda[i]/(kappa[i]*rhom[i]*drm[i]);
	}

	a[0] = 0.;
	b[0] = 1.-3.*dt/dr3[0]*r[1]*r[1]*Xfac[1];
	c[0] = 3.*dt/dr3[0]*r[1]*r[1]*Xfac[1];

	for(i = 1; i < nsize-1; i++){
		a[i] = 3.*dt/dr3[i]*r[i]*r[i]*Xfac[i];
		b[i] = 1.-3.*dt/dr3[i]*(r[i+1]*r[i+1]*Xfac[i+1]+r[i]*r[i]*Xfac[i]);
		c[i] = 3.*dt/dr3[i]*r[i+1]*r[i+1]*Xfac[i+1];
	}

	a[nsize-1] = 3.*dt/dr3[nsize-1]*r[nsize-1]*r[nsize-1]*Xfac[nsize-1];
	b[nsize-1] = 1.-3.*dt/dr3[nsize-1]*(-(P_C)*r[nsize]*r[nsize]+r[nsize-1]*r[nsize-1]*Xfac[nsize-1]);
	c[nsize-1] = 0.;
}

void init_E(double F_ini, double r[], double E[], double dt, double func[], const int nsize)
{
	double dr3 = (r[1]-r[0])*(r[1]*r[1]+r[1]*r[0]+r[0]*r[0]);
	int i;

	for(i = 0; i < nsize; i++){
		func[i] = E[2*i];
	}
	func[0] += 3.*r[0]*r[0]*F_ini/dr3*dt;
}

void itg_adv_E(double F_ini, double r[], double E[], double U[], double rho[], double dt, const int nsize)
{
	int i;
	double func[nsize];
	double a[nsize], b[nsize], c[nsize], x[nsize];

	init_E(F_ini, r, E, dt, func, nsize);
	matrix_E(r, E, U, rho, dt, a, b, c, nsize);
	trid_matrix_algorithm(a, b, c, func, x, nsize);
	for(i = 0; i < nsize; i++){
		E[2*i+1] = x[i];
	}
}

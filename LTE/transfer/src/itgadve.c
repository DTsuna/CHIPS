#include <math.h>
#include <omp.h>
#include "constant.h"
#include "opacity.h"
#include "itgadve.h"

double func_lambda(double R)
{
	return (2.+R)/(6.+3.*R+R*R);
}

void matrix_E(double r_ini, double r[], double E[], double U[], double rho[], double dt, double a[], double b[], double c[], const int nsize)
{
	double T[nsize], T_h[nsize], rho_h[nsize], mu[nsize], r_h[nsize+1];
	double kappa[nsize], R, dEdr, lambda[nsize];
	double dr[nsize], dr3;
	double x[nsize], y[nsize];
	double D0, D1;
	int i;

	for(i = 0; i < nsize; i++){
		dr[i] = r[i+1]-r[i];
	}

#pragma omp parallel for
	for(i = 0; i < nsize; i++){
		saha(rho[i], U[2*i+1], mu+i, T+i);
	}

	for(i = 1; i < nsize; i++){
		T_h[i] = (T[i]+T[i-1])/2.0;
		rho_h[i] = (rho[i]+rho[i-1])/2.0;
		r_h[i] = (r[i]+r[i-1])/2.0;
	}
	r_h[0] = r_ini;
	r_h[nsize] = r_h[nsize-1]+(r[nsize]-r[nsize-1]);

#pragma omp parallel for private(dr3, dEdr, R)
	for(i = 1; i < nsize; i++){
		dr3 = pow(r_h[i+1], 3.)-pow(r_h[i], 3.);
		x[i] = r_h[i+1]*r_h[i+1]*3.*dt/dr3;
		y[i] = r_h[i]*r_h[i]*3.*dt/dr3;
		kappa[i] = kappa_r(rho_h[i], T_h[i]);
		dEdr = (E[2*i+1]-E[2*(i-1)+1])/dr[i];
		R = fabs(dEdr)/(kappa[i]*rho_h[i]*(E[2*i+1]+E[2*(i-1)+1])/2.);
		lambda[i] = func_lambda(R);
	}
	
	dr3 = pow(r_h[1], 3.)-pow(r_h[0], 3.);
	x[0] = r_h[1]*r_h[1]*3.*dt/dr3;
	D1 = -(P_C)*lambda[1]/(kappa[1]*rho_h[1])/dr[0];

	a[0] = 0.;
	b[0] = 1.-x[0]*D1;
	c[0] = x[0]*D1;
	
	for(i = 1; i < nsize-1; i++){
		D0 = -(P_C)*lambda[i]/(kappa[i]*rho_h[i])/dr[i];
		D1 = -(P_C)*lambda[i+1]/(kappa[i+1]*rho_h[i+1])/dr[i];
		a[i] = y[i]*D0;
		b[i] = (1.-x[i]*D1-y[i]*D0);
		c[i] = x[i]*D1;
	}
	
	i = nsize-1;
	D0 = -(P_C)*lambda[i]/(kappa[i]*rho_h[i])/dr[i];
	a[i] = y[i]*D0;
	b[i] = 1.+(P_C)*x[i]*pow(r[i]/r[i+1], 2.)-y[i]*D0;
	c[i] = 0.;
}

void init_E(double r_ini, double F_ini, double r[], double E[], double dt, double func[], const int nsize)
{
	int i;
	double dr, r0, r1, dr3;
	double y;

	dr = r[1]-r[0];
	r0 = r_ini;
	r1 = r[0]+dr/2.;
	dr3 = pow(r1, 3.0)-pow(r0, 3.0);
	y = 3.*dt/dr3*r0*r0;

	for(i = 0; i < nsize; i++){
		func[i] = E[2*i];
	}
	func[0] += y*F_ini;
}

void itg_adv_E(double r_ini, double F_ini, double r[], double E[], double U[], double rho[], double dt, const int nsize)
{
	int i;
	double func[nsize];
	double a[nsize], b[nsize], c[nsize], x[nsize];

	init_E(r_ini, F_ini, r, E, dt, func, nsize);
	matrix_E(r_ini, r, E, U, rho, dt, a, b, c, nsize);
	trid_matrix_algorithm(a, b, c, func, x, nsize);
	for(i = 0; i < nsize; i++){
		E[2*i+1] = x[i];
	}
}

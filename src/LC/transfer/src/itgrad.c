#include "itgrad.h"
#include "opacity.h"
#include "function.h"
#include "constant.h"
#include "trid.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void integ_adv_eadtra(double F_ini, double E[], double T_g[], double r[], double dt, const int n, int *flag, int *loc, double *err)
{
	double a[n-1], b[n-1], c[n-1];
	double func[n-1];
	double fac = 1.;
	double x[n-1];
	double kappa[n-1];
	double eps, tol = 1.e-7;
	int i;
	int count = 0, count_max = 100;

	for(i = 0; i < n-1; i++){
		kappa[i] = kappa_r(rho_csm(r[i]), T_g[i]);
	}

	do{
		eps = 0.;
		Jacob_func_radtra(F_ini, E, T_g, kappa, r, a, b, c, func, dt, n);
		trid_matrix_algorithm(a, b, c, func, x, n-1);
		for(i = 0; i < n-1; i++){
			E[2*i+1] += fac*x[i]*E[2*i];
			if(fabs(x[i]/E[2*i+1]) > eps){
				*loc = i;
				eps = fabs(x[i]);
			}
		}
		count++;
		if(count%30 == 0){
			for(i = 0; i < n-1; i++){
				E[2*i+1] = E[2*i];
			}
			fac *= 0.7;
		}
		if(count == count_max){
			*err = eps;
			*flag = 1;
			break;
		}
	}while(eps > tol);
}

// i must not be 0, N-1
double output_R(double E[], double T_g[], double kappa[], double r[], int i)
{
	double R;
	double r_r, r_l;

	r_r = pow((r[i+1]*r[i+1]*r[i+1]+r[i]*r[i]*r[i])/2., 1./3.);
	r_l = pow((r[i-1]*r[i-1]*r[i-1]+r[i]*r[i]*r[i])/2., 1./3.);
	R = fabs(log(fabs(E[2*i+1]/E[2*(i-1)+1]))/log(r_r/r_l)/(kappa[i]*rho_csm(r[i])*r[i]));

	return fmin(R, R_MAX);
}


double Flux(double F_ini, double E[], double T_g[], double kappa[], double r[], int i, const int n)
{
	double R, lambda;
	double r_r, r_l;
	double F;

	if(i != 0 && i != n-1){
		R = output_R(E, T_g, kappa, r, i);
		lambda = (2.+R)/(6.+3.*R+R*R);
		r_r = pow((r[i+1]*r[i+1]*r[i+1]+r[i]*r[i]*r[i])/2., 1./3.);
		r_l = pow((r[i-1]*r[i-1]*r[i-1]+r[i]*r[i]*r[i])/2., 1./3.);

		F = (P_C)*lambda*(E[2*i+1]+E[2*(i-1)+1])/2.*R;
	}
	else if(i == 0){
		F = F_ini;
	}
	else if(i == n-1){
		F = (P_C)*E[2*(n-2)+1];
	}
	else{
		printf("N is out of range. stop.\n");
		exit(EXIT_FAILURE);
	}

	return F;
}

// i should be smaller than N-1
double func_radtra(double F_ini, double E[], double T_g[], double kappa[], double r[], double dt, int i, const int n)
{
	double func;
	double dr3;

	func = E[2*i+1]-E[2*i];
	dr3 = (r[i+1]-r[i])*(r[i+1]*r[i+1]+r[i+1]*r[i]+r[i]*r[i]);

	func += 3.*dt/dr3*(r[i+1]*r[i+1]*Flux(F_ini, E, T_g, kappa, r, i+1, n));
	func += -3.*dt/dr3*(r[i]*r[i]*Flux(F_ini, E, T_g, kappa, r, i, n));

	return func/E[2*i];
}

// output a, b, c
void Jacob_func_radtra(double F_ini, double E[], double T_g[], double kappa[], double r[], double a[], double b[], double c[], double func[], double dt, const int n)
{
	double dE[n-1];
	double eps = 1.e-5;
	double f4, f3, f2, f1, f0; //f4 = f(x-2dx), f3 = f(x+2dx), f2 = f(x-dx), f1 = f(x+dx), f0 = f(x).
	int i;

	for(i = 0; i < n-1; i++){
		dE[i] = eps*E[2*i];
	}

// set a[i]
	a[0] = 0.;
	for(i = 1; i < n-1; i++){
		f0 = func_radtra(F_ini, E, T_g, kappa, r, dt, i, n);
		E[2*(i-1)+1] += dE[i-1];
		f1 = func_radtra(F_ini, E, T_g, kappa, r, dt, i, n);
		E[2*(i-1)+1] -= 2.*dE[i-1];
		f2 = func_radtra(F_ini, E, T_g, kappa, r, dt, i, n);
		E[2*(i-1)+1] -= dE[i-1];
		f4 = func_radtra(F_ini, E, T_g, kappa, r, dt, i, n);
		E[2*(i-1)+1] += 4.*dE[i-1];
		f3 = func_radtra(F_ini, E, T_g, kappa, r, dt, i, n);
		E[2*(i-1)+1] -= 2.*dE[i-1];
		a[i] = (8.*(f1-f2)-(f3-f4))/(12.*dE[i-1])*E[2*(i-1)];
	}


// set b[i]
	for(i = 0; i < n-1; i++){
		f0 = func_radtra(F_ini, E, T_g, kappa, r, dt, i, n);
		func[i] = -f0;
		E[2*i+1] += dE[i];
		f1 = func_radtra(F_ini, E, T_g, kappa, r, dt, i, n);
		E[2*i+1] -= 2.*dE[i];
		f2 = func_radtra(F_ini, E, T_g, kappa, r, dt, i, n);
		E[2*i+1] -= dE[i];
		f4 = func_radtra(F_ini, E, T_g, kappa, r, dt, i, n);
		E[2*i+1] += 4.*dE[i];
		f3 = func_radtra(F_ini, E, T_g, kappa, r, dt, i, n);
		E[2*i+1] -= 2.*dE[i];
		b[i] = (8.*(f1-f2)-(f3-f4))/(12.*dE[i])*E[2*i];
	}


// set c[i]
	for(i = 0; i < n-2; i++){
		f0 = func_radtra(F_ini, E, T_g, kappa, r, dt, i, n);
		E[2*(i+1)+1] += dE[i+1];
		f1 = func_radtra(F_ini, E, T_g, kappa, r, dt, i, n);
		E[2*(i+1)+1] -= 2.*dE[i+1];
		f2 = func_radtra(F_ini, E, T_g, kappa, r, dt, i, n);
		E[2*(i+1)+1] -= dE[i+1];
		f4 = func_radtra(F_ini, E, T_g, kappa, r, dt, i, n);
		E[2*(i+1)+1] += 4.*dE[i+1];
		f3 = func_radtra(F_ini, E, T_g, kappa, r, dt, i, n);
		E[2*(i+1)+1] -= 2.*dE[i+1];
		c[i] = (8.*(f1-f2)-(f3-f4))/(12.*dE[i+1])*E[2*(i+1)];
	}
	c[n-2] = 0.;
}

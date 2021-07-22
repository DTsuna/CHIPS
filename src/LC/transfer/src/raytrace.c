#include "raytrace.h"
#include "saha.h"
#include "constant.h"
#include <stdio.h>
#include <math.h>

int imax(int a, int b)
{
	if(a > b){
		return a;
	}
	else{
		return b;
	}
}

double ds_path(double b, double r[], int j)
{
	double r0, r1;

	r0 = sqrt(r[j]*r[j]-b*b);
	r1 = sqrt(r[j+1]*r[j+1]-b*b);

	if(r[j] > b){
		return (r[j+1]+r[j])*(r[j+1]-r[j])/(r0+r1);
	}
	else{
		return r1;
	}
}

//output of (kappa+sigma)*rho
double alpha_ff_plus_beta_nu(double nu, double rho, double T)
{
	double ndens[4];
	double x = (P_H)*nu/((P_K)*T);

	get_num_density(rho, T,  ndens);
	
	return 3.692e+08*ndens[0]*(ndens[1]+ndens[2]+4.0*ndens[3])/sqrt(T)*pow(nu, -3.0)*(1.-exp(-x))+ndens[0]*(SIGMA_TH);
}

double Planck_func(double nu, double T)
{
	double fac1, fac2;

	fac1 = 2.0*(P_H)/((P_C)*(P_C))*nu*nu*nu;
	fac2 = (P_H)*nu/((P_K)*T);

	if(fac2 < 1.0e-05){
		return fac1/(fac2+fac2*(0.5*fac2+(fac2*fac2)/6.0));
	}
	else{
		return fac1/(exp(fac2)-1.0);
	}
}

/************************This subroutine returns I(b, nu)**************************/
/**********************************************************************************/
double integ_ray_tracing(double b, double nu, double r[], double rho[], double T[], double r_sh[], double rho_sh[], double T_sh[], int n, int n_sh)
{
	int i, j, jmin, jmin_sh;
	double B_nu;
	double fac;
	double ds;
	double I = 0.;
	
	jmin = jmin_func(b, r, n);

//Integrate intensity from the outermost region
	for(j = n-1; j >= imax(jmin, 0); --j){
		fac = alpha_ff_plus_beta_nu(nu, rho[j], T[j]);
		ds = ds_path(b, r, j);
		B_nu = Planck_func(nu, T[j]);
		I = (1.-0.5*fac*ds)/(1.+0.5*fac*ds)*I+fac*B_nu/(1.+0.5*fac*ds)*ds;
	}

	if(jmin == -1){
		jmin_sh = jmin_func(b, r_sh, n_sh);

		for(j = n_sh-1; j >= imax(jmin_sh, 0); --j){
			fac = alpha_ff_plus_beta_nu(nu, rho_sh[j], T_sh[j]);
			ds = ds_path(b, r_sh, j);
			B_nu = Planck_func(nu, T_sh[j]);
			I = (1.-0.5*fac*ds)/(1.+0.5*fac*ds)*I+fac*B_nu/(1.+0.5*fac*ds)*ds;
		}
		for(j = imax(jmin_sh, 0); j < n_sh; j++){
			fac = alpha_ff_plus_beta_nu(nu, rho_sh[j], T_sh[j]);
			ds = ds_path(b, r_sh, j);
			B_nu = Planck_func(nu, T_sh[j]);
			I = (1.-0.5*fac*ds)/(1.+0.5*fac*ds)*I+fac*B_nu/(1.+0.5*fac*ds)*ds;
		}
	}



	for(j = imax(jmin, 0); j < n; j++){
		fac = alpha_ff_plus_beta_nu(nu, rho[j], T[j]);
		ds = ds_path(b, r, j);
		B_nu = Planck_func(nu, T[j]);
		I = (1.-0.5*fac*ds)/(1.+0.5*fac*ds)*I+fac*B_nu/(1.+0.5*fac*ds)*ds;
	}
	if(isnan(I)){
		printf("b, nu = %e %e\n", b, nu);
	}

	return I;
}

double Lum_nu(double r_init, double r_out, double nu, double r[], double rho[], double T[], double r_sh[], double rho_sh[], double T_sh[], int n, int n_sh)
{
	int i;
	double b[NB] = {0.}, db = r_out/(double)((NB)-1);
	double I[NB], sum = 0.;
	double fac;

	b[0] = 0.;
	b[1] = r_init/100.;
	fac = pow(r_out/b[1], 1./((double)NB-2.));
	for(i = 2; i < NB; i++){
		b[i] = fac*b[i-1];
	}
	
	I[0] = integ_ray_tracing(b[0], nu, r, rho, T, r_sh, rho_sh, T_sh, n, n_sh);

	for(i = 1; i < NB-1; i++){
//		b[i] = b[i-1]+db;
		I[i] = integ_ray_tracing(b[i], nu, r, rho, T, r_sh, rho_sh, T_sh, n, n_sh);
	}
	b[NB-1] = b[NB-2]+(b[NB-1]-b[NB-2])/2.;
	I[NB-1] = integ_ray_tracing(b[NB-1], nu, r, rho, T, r_sh, rho_sh, T_sh, n, n_sh);

	for(i = 0; i < NB-1; i++){
		sum += 8.*M_PI*M_PI*0.5*(b[i]*I[i]+b[i+1]*I[i+1])*(b[i+1]-b[i]);
	}

	return sum;
}

void calc_lum(double r_init, double r_out, double r[], double rho[], double T[], double r_sh[], double rho_sh[], double T_sh[], int n, int n_sh, char *filename)
{
	FILE *fp;
	double nu[NNU], L_nu[NNU];
	double frac;
	int i;

	nu[0] = (P_C)/(3.e-05); //corresponds to lambda = 300nm;
	nu[NNU-1] = (P_C)/(1.e-04); //corresponds to lambda = 1000nm;

	frac = pow(nu[NNU-1]/nu[0], 1./((double)NNU-1.0));
	for(i = 1; i < NNU; i++){
		nu[i] = frac*nu[i-1];
	}

	printf("Started spec calculation.\n");
	for(i = 0; i < NNU; i++){
		L_nu[i] = Lum_nu(r_init, r_out, nu[i], r, rho, T, r_sh, rho_sh, T_sh, n, n_sh);
		printf("L_nu[%d] = %e\n", i, L_nu[i]);
	}
	printf("Calculation end.\n");

	fp = fopen(filename, "w");
	for(i = 0; i < NNU; i++){
		fprintf(fp, "%e %e\n", nu[i], L_nu[i]);
	}
	fclose(fp);
}

int jmin_func(double b, double r[], int n)
{
	int i, jmin = 0;

	if(b < r[0]){
		jmin = -1;
	}
	else{
		for(i = 0; i < n; i++){
			if(b >= r[i] && b < r[i+1]){
				jmin = i;
				break;
			}
		}
	}
	return jmin;
}

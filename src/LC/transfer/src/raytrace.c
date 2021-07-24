#include "raytrace.h"
#include "saha.h"
#include "constant.h"
#include <stdio.h>
#include <omp.h>
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
	double alpha;

	get_num_density(rho, T,  ndens);
	alpha = 3.692e+08*ndens[0]*(ndens[1]+ndens[2]+4.0*ndens[3])/sqrt(T)*pow(nu, -3.0)*(1.-exp(-x))+ndens[0]*(SIGMA_TH);
	if(isnan(alpha)){
		printf("nan (alpha), exp(-x) = %e, exp(x)-1 = %e\n", exp(-x), expm1(x));
	}

	return alpha;
}

double Planck_func(double nu, double T)
{
	double fac1, fac2;

	fac1 = 2.0*(P_H)/((P_C)*(P_C))*nu*nu*nu;
	fac2 = (P_H)*nu/((P_K)*T);

	if(fac2 < 15.){
		return fac1/expm1(fac2);
	}
	else{
		return fac1*exp(-fac2);
	}
}

/************************This subroutine returns I(b, nu)**************************/
/**********************************************************************************/
/*
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
*/

double integ_ray_tracing(double b, double nu, double r[], double rho[], double T[], double r_sh[], double rho_sh[], double T_sh[], int n, int n_sh)
{
	int i, j, jmin, jmin_sh;
	double B_nu;
	double fac;
	double ds, tau = 0., dtau, tau_fin;
	double sum = 0.;
	double I = 0.;
	
	jmin = jmin_func(b, r, n);

//Integrate intensity from the outermost region
	for(j = n-1; j >= imax(jmin, 0); --j){
		fac = alpha_ff_plus_beta_nu(nu, rho[j], T[j]);
		ds = ds_path(b, r, j);
		B_nu = Planck_func(nu, T[j]);
		dtau = fac*ds;
		tau += dtau;
	}

	if(jmin == -1){
		jmin_sh = jmin_func(b, r_sh, n_sh);

		for(j = n_sh-1; j >= imax(jmin_sh, 0); --j){
			fac = alpha_ff_plus_beta_nu(nu, rho_sh[j], T_sh[j]);
			ds = ds_path(b, r_sh, j);
			B_nu = Planck_func(nu, T_sh[j]);
			dtau = fac*ds;
			tau += dtau;
		}
		for(j = imax(jmin_sh, 0); j < n_sh; j++){
			fac = alpha_ff_plus_beta_nu(nu, rho_sh[j], T_sh[j]);
			ds = ds_path(b, r_sh, j);
			B_nu = Planck_func(nu, T_sh[j]);
			dtau = fac*ds;
			tau += dtau;
		}
	}



	for(j = imax(jmin, 0); j < n; j++){
		fac = alpha_ff_plus_beta_nu(nu, rho[j], T[j]);
		ds = ds_path(b, r, j);
		B_nu = Planck_func(nu, T[j]);
		dtau = fac*ds;
		tau += dtau;
	}

	tau_fin = tau;
	tau = 0.;



	for(j = n-1; j >= imax(jmin, 0); --j){
		fac = alpha_ff_plus_beta_nu(nu, rho[j], T[j]);
		ds = ds_path(b, r, j);
		B_nu = Planck_func(nu, T[j]);
		dtau = fac*ds;
		tau += dtau;
		sum += B_nu*exp(tau-tau_fin)*dtau;
	}

	if(jmin == -1){
		jmin_sh = jmin_func(b, r_sh, n_sh);

		for(j = n_sh-1; j >= imax(jmin_sh, 0); --j){
			fac = alpha_ff_plus_beta_nu(nu, rho_sh[j], T_sh[j]);
			ds = ds_path(b, r_sh, j);
			B_nu = Planck_func(nu, T_sh[j]);
			dtau = fac*ds;
			tau += dtau;
			sum += B_nu*exp(tau-tau_fin)*dtau;
		}
		for(j = imax(jmin_sh, 0); j < n_sh; j++){
			fac = alpha_ff_plus_beta_nu(nu, rho_sh[j], T_sh[j]);
			ds = ds_path(b, r_sh, j);
			B_nu = Planck_func(nu, T_sh[j]);
			dtau = fac*ds;
			tau += dtau;
			sum += B_nu*exp(tau-tau_fin)*dtau;
		}
	}



	for(j = imax(jmin, 0); j < n; j++){
		fac = alpha_ff_plus_beta_nu(nu, rho[j], T[j]);
		ds = ds_path(b, r, j);
		B_nu = Planck_func(nu, T[j]);
		dtau = fac*ds;
		tau += dtau;
		sum += B_nu*exp(tau-tau_fin)*dtau;
	}

	return sum;
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

void calc_lum(double r_init, double r_out, double r[], double rho[], double T[], 
		double r_sh[], double rho_sh[], double T_sh[], int n, int n_sh, char *filename, double abmag[])
{
	const double pc = 3.085677581e+18;
	FILE *fp, *fb;
	double nu[NNU], L_nu[NNU];
	double lam[NNU];
	double frac;
	double lam_band[100], trans_band[100], Lnu_band[100];
	double dummy[3];
	double sum1, sum2;
	int i, j, k, l, lmin;
	char bands[5][256] = {"./input/band_filters/uband.txt", "./input/band_filters/bband.txt", "./input/band_filters/vband.txt",
				"./input/band_filters/rband.txt", "./input/band_filters/iband.txt"};

	nu[0] = (P_C)/(2.9e-05); //corresponds to lambda = 290 nm;
	nu[NNU-1] = (P_C)/(1.e-04); //corresponds to lambda = 1000 nm;

	frac = pow(nu[NNU-1]/nu[0], 1./((double)NNU-1.0));
	for(i = 1; i < NNU; i++){
		nu[i] = frac*nu[i-1];
	}

	printf("Started spec calculation.\n");
#pragma omp parallel for
	for(i = 0; i < NNU; i++){
		L_nu[i] = Lum_nu(r_init, r_out, nu[i], r, rho, T, r_sh, rho_sh, T_sh, n, n_sh);
		printf("L_nu[%d] = %e\n", i, L_nu[i]);
	}
	printf("Spec calculation end.\n");

	fp = fopen(filename, "w");
	for(i = 0; i < NNU; i++){
		lam[i] = (P_C)/nu[i];
		fprintf(fp, "%e %e\n", nu[i], L_nu[i]);
	}
	fclose(fp);

	for(i = 0; i < 5; i++){
		j = 0;
		sum1 = 0.;
		sum2 = 0.;
		fb = fopen(bands[i], "r");
		while(fscanf(fb, "%lf %lf %lf %lf %lf", lam_band+j, dummy, dummy+1, dummy+2, trans_band+j) != EOF){
			lam_band[j] *= 1.0e-07;
			j++;
		}
		for(k = 0; k < j; k++){
			lmin = 0;
			for(l = lmin; l < NNU-1; l++){
				if(lam_band[k] >= lam[l] && lam_band[k] < lam[l+1]){
					lmin = l;
					Lnu_band[k] = (L_nu[l+1]-L_nu[l])/(lam[l+1]-lam[l])*(lam_band[k]-lam[l])+L_nu[l];
					break;
				}
			}
		}

		for(k = 0; k < j-1; k++){
			sum1 += (trans_band[k+1]+trans_band[k])/(lam_band[k+1]+lam_band[k])*(lam_band[k+1]-lam_band[k]);
			sum2 += (Lnu_band[k+1]+Lnu_band[k])/2.0*(trans_band[k+1]+trans_band[k])/(lam_band[k+1]+lam_band[k])*(lam_band[k+1]-lam_band[k]);
		}
		sum2 /= 4.*M_PI*100.*pc*pc;
		fclose(fb);
		abmag[i] = -2.5*log10(sum2/sum1)-48.6;
	}
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

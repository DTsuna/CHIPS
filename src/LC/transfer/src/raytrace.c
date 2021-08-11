#include "raytrace.h"
#include "dcht.h"
#include "saha.h"
#include "opacity.h"
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

double alpha_p_beta_nu(double nu, double rho, double T, opacity op)
{
	int i, j;
	double a, b, c, d;
	double kappa;
	double ndens[4];
	double alpha, beta;
	double R = rho/sqrt(T)/T;

	get_num_density(rho, T,  ndens);
	beta = ndens[0]*(SIGMA_TH);

	double logT = log10(T);
	double logR = log10(R);

	if(logT < op.T[0]){
		if(logR < op.R[0]){
			return rho*pow(10., op.kappa[0])+beta;
		}
		else if(logR > op.R[0] && logR < op.R[op.jmax-1]){
			j = dcht(logR, op.R, op.jmax);
			kappa = op.kappa[j]+(op.kappa[j+1]-op.kappa[j])/(op.R[j+1]-op.R[j])*(logR-op.R[j]);
			return rho*pow(10., kappa)+beta;
		}
		else{
			return rho*pow(10., op.kappa[op.jmax-1])+beta;
		}
	}
	else if(logR > op.R[op.jmax-1]){
		i = dcht(logT, op.T, op.imax);
		kappa = op.kappa[op.jmax*i+op.jmax-1]
		+(op.kappa[op.jmax*(i+1)+op.jmax-1]-op.kappa[op.jmax*i+op.jmax-1])/(op.T[i+1]-op.T[i])*(logT-op.T[i]);
		kappa = pow(10., kappa);
		kappa *= R/pow(10., op.R[op.jmax-1]);
		return rho*kappa+beta;
	}
	else if(logT > op.T[op.imax-1]){
		return rho*pow(10., op.kappa[op.jmax*(op.imax-1)])+beta;
	}
	else if(logT > op.T[0] && logT < op.T[op.imax-1] && logR < op.R[0]){
		i = dcht(logT, op.T, op.imax);
		kappa = op.kappa[op.jmax*i]+(op.kappa[op.jmax*(i+1)]-op.kappa[op.jmax*i])/(op.T[i+1]-op.T[i])*(logT-op.T[i]);
		return rho*pow(10., kappa)+beta;
	}
	else{
		i = dcht(logT, op.T, op.imax);
		j = dcht(logR, op.R, op.jmax);
		a = ((op.kappa[op.jmax*(i+1)+j])-(op.kappa[op.jmax*i+j]))/(op.T[i+1]-op.T[i]);
		b = ((op.kappa[op.jmax*i+j+1])-(op.kappa[op.jmax*i+j]))/(op.R[j+1]-op.R[j]);
		c = ((op.kappa[op.jmax*(i+1)+j+1])-(op.kappa[op.jmax*(i+1)+j])-(op.kappa[op.jmax*i+j+1])
			+(op.kappa[op.jmax*i+j]))/((op.T[i+1]-op.T[i])*(op.R[j+1]-op.R[j]));
		kappa = op.kappa[op.jmax*i+j]+a*(logT-op.T[i])+b*(logR-op.R[j])+c*(logT-op.T[i])*(logR-op.R[j]);
		return rho*pow(10., kappa)+beta;
	}
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
double integ_ray_tracing(double b, double nu, double r[], double rho[], double T[], double r_sh[], double rho_sh[], double T_sh[], int n, int n_sh, opacity op)
{
	FILE *fnu;
	int i, j, jmin, jmin_sh, inu = 0;
	int jminmax, jminmax_sh;
	int k = 0, l;
	double B_nu;
	double fac;
	double ds, tau = 0., dtau, tau_fin;
	double sum = 0.;
	double I = 0.;
	double opacity[8192], ds_array[8192], Planck[8192];
	
	jmin = jmin_func(b, r, n);

//Integrate intensity from the outermost region
	jminmax = imax(jmin, 0);
	for(j = n-1; j >= jminmax; --j){
		opacity[k] = alpha_p_beta_nu(nu, rho[j], T[j], op);
		fac = opacity[k];
		ds_array[k] = ds_path(b, r, j);
		Planck[k] = Planck_func(nu, T[j]);
		ds = ds_array[k];
		B_nu = Planck[k];
		k++;
		dtau = fac*ds;
		tau += dtau;
	}

	if(jmin == -1){
		jmin_sh = jmin_func(b, r_sh, n_sh);
		jminmax_sh = imax(jmin_sh, 0);

		for(j = n_sh-1; j >= jminmax_sh; --j){
			opacity[k] = alpha_p_beta_nu(nu, rho_sh[j], T_sh[j], op);
			ds_array[k] = ds_path(b, r_sh, j);
			Planck[k] = Planck_func(nu, T_sh[j]);
			fac = opacity[k];
			ds = ds_array[k];
			B_nu = Planck[k];
			k++;
			dtau = fac*ds;
			tau += dtau;
		}
		for(j = jminmax_sh; j < n_sh; j++){
			opacity[k] = alpha_p_beta_nu(nu, rho_sh[j], T_sh[j], op);
			fac = opacity[k];
			ds_array[k] = ds_path(b, r_sh, j);
			Planck[k] = Planck_func(nu, T_sh[j]);
			fac = opacity[k];
			ds = ds_array[k];
			B_nu = Planck[k];
			k++;
			dtau = fac*ds;
			tau += dtau;
		}
	}



	for(j = jminmax; j < n; j++){
		opacity[k] = alpha_p_beta_nu(nu, rho[j], T[j], op);
		fac = opacity[k];
		ds_array[k] = ds_path(b, r, j);
		Planck[k] = Planck_func(nu, T[j]);
		ds = ds_array[k];
		B_nu = Planck[k];
		k++;
		dtau = fac*ds;
		tau += dtau;
	}

	tau_fin = tau;
	tau = 0.;
	k = 0;


	for(j = n-1; j >= jminmax; --j){
		fac = opacity[k];
		ds = ds_array[k];
		B_nu = Planck[k];
		k++;
		dtau = fac*ds;
		tau += dtau;
		sum += B_nu*exp(tau-tau_fin)*dtau;
	}

	if(jmin == -1){
		jmin_sh = jmin_func(b, r_sh, n_sh);

		for(j = n_sh-1; j >= jminmax_sh; --j){
			fac = opacity[k];
			ds = ds_array[k];
			B_nu = Planck[k];
			k++;
			dtau = fac*ds;
			tau += dtau;
			sum += B_nu*exp(tau-tau_fin)*dtau;
		}
		for(j = jminmax_sh; j < n_sh; j++){
			fac = opacity[k];
			ds = ds_array[k];
			B_nu = Planck[k];
			k++;
			dtau = fac*ds;
			tau += dtau;
			sum += B_nu*exp(tau-tau_fin)*dtau;
		}
	}



	for(j = jminmax; j < n; j++){
		fac = opacity[k];
		ds = ds_array[k];
		B_nu = Planck[k];
		k++;
		dtau = fac*ds;
		tau += dtau;
		sum += B_nu*exp(tau-tau_fin)*dtau;
	}

	return sum;
}

double Lum_nu(double r_init, double r_out, double nu, double r[], double rho[], double T[], double r_sh[], double rho_sh[], double T_sh[], int n, int n_sh)
{
	int i, inu = 0, k, l;
	double b[NB] = {0.}, db = r_out/(double)((NB)-1);
	double I[NB], sum = 0.;
	double fac;
	double kappa[2000];
	double nu_opac[2000];
	FILE *fnu;
	char filename[512];

/***************************Read nu******************************/
	opacity op0, op1;
	fnu = fopen("input/TOPS_multigroup/opacity_table/frequency.txt", "r");
	while(fscanf(fnu, "%lf", &nu_opac[inu]) != EOF){
		inu++;
	}
	fclose(fnu);
	for(i = 0; i < inu-1; i++){
		if(nu >= nu_opac[i] && nu < nu_opac[i+1]){
			break;
		}
	}
	sprintf(filename, "./input/TOPS_multigroup/opacity_table/opacity_%05d.txt", i);
	set_opacity(filename, &op0);
	sprintf(filename, "./input/TOPS_multigroup/opacity_table/opacity_%05d.txt", i+1);
	set_opacity(filename, &op1);
	for(k = 0; k < op0.jmax; k++){
		for(l = 0; l < op0.imax; l++){
			kappa[l*op0.jmax+k] = ((nu_opac[i+1]-nu)*op0.kappa[l*op0.jmax+k]+(nu-nu_opac[i])*op1.kappa[l*op1.jmax+k])/(nu_opac[i+1]-nu_opac[i]);
			op0.kappa[l*op0.jmax+k] = kappa[l*op0.jmax+k];
		}
	}
/****************************************************************/

	b[0] = 0.;
	b[1] = r_init/100.;
	fac = pow(r_out/b[1], 1./((double)NB-2.));
	for(i = 2; i < NB; i++){
		b[i] = fac*b[i-1];
	}
	
	I[0] = integ_ray_tracing(b[0], nu, r, rho, T, r_sh, rho_sh, T_sh, n, n_sh, op0);

#pragma omp parallel for
	for(i = 1; i < NB-1; i++){
//		b[i] = b[i-1]+db;
		I[i] = integ_ray_tracing(b[i], nu, r, rho, T, r_sh, rho_sh, T_sh, n, n_sh, op0);
	}
	b[NB-1] = b[NB-2]+(b[NB-1]-b[NB-2])/2.;
	I[NB-1] = integ_ray_tracing(b[NB-1], nu, r, rho, T, r_sh, rho_sh, T_sh, n, n_sh, op0);

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
	double Trans[NNU];
	double frac;
	double lam_band[100], trans_band[100], Lnu_band[100];
	double dummy[3];
	double sum1, sum2;
	int i, j, k, l;
	char bands[5][256] = {"./input/band_filters/uband.txt", "./input/band_filters/bband.txt", "./input/band_filters/vband.txt",
				"./input/band_filters/rband.txt", "./input/band_filters/iband.txt"};

	nu[0] = (P_C)/(2.9e-05); //corresponds to lambda = 290 nm;
	nu[NNU-1] = (P_C)/(1.e-04); //corresponds to lambda = 1000 nm;

	frac = pow(nu[NNU-1]/nu[0], 1./((double)NNU-1.0));
	for(i = 1; i < NNU; i++){
		nu[i] = frac*nu[i-1];
	}

	printf("Started spec calculation.\n");
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

/*
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
*/

	for(i = 0; i < 5; i++){
		j = 0;
		sum1 = 0.;
		sum2 = 0.;
		fb = fopen(bands[i], "r");
		while(fscanf(fb, "%lf %lf %lf %lf %lf", lam_band+j, dummy, dummy+1, dummy+2, trans_band+j) != EOF){
			lam_band[j] *= 1.0e-07;
			j++;
		}

		for(l = 0; l < NNU-1; l++){
			if(lam[l] <= lam_band[0] || lam[l] >= lam_band[j-1]){
				Trans[l] = 0.;
			}
			else{
				for(k = 0; k < j-1; k++){
					if(lam[l] >= lam_band[k] && lam[l] < lam_band[k+1]){
						Trans[l] = (trans_band[k+1]-trans_band[k])/(lam_band[k+1]-lam_band[k])*(lam[l]-lam_band[k])+trans_band[k];
						break;
					}
				}
			}
			printf("Trans[%d] = %f\n", l, Trans[l]);
		}
		for(l = 0; l < NNU-1; l++){
			sum1 += (Trans[l]+Trans[l+1])/(lam[l]+lam[l+1])*(lam[l+1]-lam[l]);
			sum2 += (L_nu[l]+L_nu[l+1])/2.*(Trans[l]+Trans[l+1])/(lam[l]+lam[l+1])*(lam[l+1]-lam[l]);
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

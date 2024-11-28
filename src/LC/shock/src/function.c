#include <stdio.h>
#include <math.h>
#include "opacity.h"
#include "function.h"
#include "constant.h"
#include "pars.h"
#include "saha.h"

char csm[256];
extern pars pdt;
extern double X, Y;

/*Analytical solution of the thin shell model by Moriya et al. (2013)*/
double r_early(double t)
{
	double A, B, C;
	double n, s, delta, M_ej, E_ej, D;
	double r_out = 1.5e+14;
	double rho_out = rho_csm(r_out);

	n = pdt.n;
	s = pdt.s;
	delta = pdt.delta;
	M_ej = pdt.M_ej;
	E_ej = pdt.E_ej;

	D = rho_out*pow(r_out, s);

	A = pow(2.*(5.-delta)*(n-5.)*E_ej, (n-3.)/2.);
	B = pow((3.-delta)*(n-3.)*M_ej, (n-5.)/2.);
	C = pow((3.-s)*(4.-s)/(4.*M_PI*D*(n-4.)*(n-3.)*(n-delta))*(A/B), 1./(n-s));
	return C*pow(t, (n-3.)/(n-s));
}

double t_early(double r)
{
	double A, B, C;
	double n, s, delta, M_ej, E_ej, D;
	double r_out = 1.5e+14;
	double rho_out = rho_csm(r_out);

	n = pdt.n;
	s = pdt.s;
	delta = pdt.delta;
	M_ej = pdt.M_ej;
	E_ej = pdt.E_ej;

	D = rho_out*pow(r_out, s);

	A = pow(2.*(5.-delta)*(n-5.)*E_ej, (n-3.)/2.);
	B = pow((3.-delta)*(n-3.)*M_ej, (n-5.)/2.);
	C = pow((3.-s)*(4.-s)/(4.*M_PI*D*(n-4.)*(n-3.)*(n-delta))*(A/B), 1./(n-s));
	return pow(r/C, (n-s)/(n-3.));
}

/*Set CSM density profile output by the numerical simulation by Kuriyama & Shigeyama (2020)*/
double rho_csm(double r)
{
	static double r_c[1000], rho_c[1000], s;
	static int flag = 0, nsize;

	int i = 0;
	double rho;
	double dummy[5];


	if(flag == 0){
		FILE *fp;
		char filename[256];
		fp = fopen(csm, "r");
		fgets(filename, sizeof(filename), fp);
		while(fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &dummy[0], &dummy[1], &r_c[i], &dummy[2], &rho_c[i], &dummy[3], &dummy[4]) != EOF){
			i++;
		}
		nsize = i;
		flag = 1;
		fclose(fp);
		s = pdt.s;
	}
	if(r <= r_c[0]){
		rho = rho_c[0]*pow(r/r_c[0], -s);
	}
	else{
		for(i = 0; i < nsize-1; i++){
			if(r < r_c[i+1] && r >= r_c[i]){
				break;
			}
		}
		rho = (log(rho_c[i+1])-log(rho_c[i]))/(log(r_c[i+1])-log(r_c[i]))*(log(r)-log(r_c[i]))+log(rho_c[i]);
		rho = exp(rho);
	}

	return rho;
}

double set_r_ini(const char *file_csm)
{
	FILE *fp;
	char filename[256];
	double dummy[7];
	fp = fopen(file_csm, "r");
	fgets(filename, sizeof(filename), fp);
	fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &dummy[0], &dummy[1], &dummy[2], &dummy[3], &dummy[4], &dummy[5], &dummy[6]);

	fclose(fp);
	return dummy[2];
}

double set_r_diff(const char *file_csm)
{
	FILE *fp;
	char filename[256];
	double dummy[7];
	double n, s, delta, M_ej, E_ej, kappa, q;
	double A, g_to_n;
	double v_sh;
	double *rho, *r, *dr, *tau, tau_tot = 0.;
	int i = 0, j;

	rho = (double*)calloc(20000, sizeof(double));
	tau = (double*)calloc(20000, sizeof(double));
	r = (double*)calloc(20000, sizeof(double));
	dr = (double*)calloc(20000, sizeof(double));

	fp = fopen(file_csm, "r");
	fgets(filename, sizeof(filename), fp);
	while(fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &dummy[0], &dummy[1], &r[i], &dummy[3], &rho[i], &dummy[5], &dummy[6]) != EOF){
		if(i != 0){
			kappa = 0.2 * (1.+dummy[5]);
			dr[i] = r[i]-r[i-1];
			tau[i] = kappa*(rho[i]+rho[i-1])/2.*dr[i];
		}
		i++;
	}

	n = pdt.n;
	delta = pdt.delta;
	M_ej = pdt.M_ej;
	E_ej = pdt.E_ej;
	s = pdt.s;
	interp_self_similar_values(&A, dummy, dummy+1);
	q = rho[0]*pow(r[0], s);
	g_to_n = pow(2.0*(5.0-delta)*(n-5.0)*E_ej, (n-3.0)/2.0)/pow((3.0-delta)*(n-3.0)*M_ej, (n-5.0)/2.0)/((n-delta)*4.*M_PI);

	for(j = i-1; j > 0; --j){
		tau_tot += tau[j];
		v_sh = (n-3.)/(n-s)*pow(A*g_to_n/q, 1./(n-3.))*pow(r[j], -(3.-s)/(n-3.));
		if((P_C)/tau_tot < v_sh){
			break;
		}
	}

	free(rho);
	free(tau);
	free(r);
	free(dr);
	fclose(fp);
	return r[j];
}

double rho_ej(double r, double t)
{
	double A, B;
	double n, delta, M_ej, E_ej;
	n = pdt.n;
	delta = pdt.delta;
	M_ej = pdt.M_ej;
	E_ej = pdt.E_ej;
	if(t > r/pdt.v_t){
		A = pow(2.0*(5.0-delta)*(n-5.0)*E_ej, (delta-3.)/2.);
		B = pow((3.0-delta)*(n-3.0)*M_ej, (delta-5.)/2.);
		return 1.0/(4.0*M_PI*(n-delta))*(A/B)*pow(t, delta-3.0)*pow(r, -delta);
	}
	else{
		A = pow(2.0*(5.0-delta)*(n-5.0)*E_ej, (n-3.)/2.);
		B = pow((3.0-delta)*(n-3.0)*M_ej, (n-5.)/2.);
		return 1.0/(4.0*M_PI*(n-delta))*(A/B)*pow(t, n-3.0)*pow(r, -n);
	}
}

double func_M_ej(double r, double t)
{
	double A, B;
	double n, delta, M_ej, E_ej;

	n = pdt.n;
	delta = pdt.delta;
	M_ej = pdt.M_ej;
	E_ej = pdt.E_ej;
	A = 2.0*(5.0-delta)*(n-5.0);
	B = (3.0-delta)*(n-3.0);
	if(r > pdt.v_t*t){
		return pow(t/r,n-3.0)*pow(A*E_ej, (n-3.0)/2.0)/pow(B*M_ej, (n-5.0)/2.0)/((n-delta)*(n-3.0));
	}
	else{
		return pow(1./pdt.v_t, n-3.0)*pow(A*E_ej, (n-3.0)/2.0)/pow(B*M_ej, (n-5.0)/2.0)/((n-delta)*(n-3.0))
		+(pow(A*E_ej, (delta-3.)/2.)/pow(B*M_ej, (delta-5.)/2.))*(pow(pdt.v_t, -delta+3.)-pow(r/t, -delta+3.))/((n-delta)*(n-3.));
	}
}

double func_M_csm(double r, double t)
{
	static double r_c[3], rho_c[3], s;
	static int flag = 0;
	int i;
	double dammy[5];

	if(flag == 0){
		FILE *fp;
		char filename[256];
		fp = fopen(csm, "r");
		fgets(filename, sizeof(filename), fp);
		printf("%s\n", filename);
		for(i = 0; i < 3; i++){
			fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &dammy[0], &dammy[1], &r_c[i], &dammy[2], &rho_c[i], &dammy[3], &dammy[4]);
		}
		s = pdt.s;
		fclose(fp);
		flag = 1;
	}

		
	return 4.0*M_PI/(3.0-s)*r*r*r*rho_csm(r);
}

double p_tot(double y[])
{
	double Pg, Pr;
	Pg = y[1]*(P_K)*y[2]/((MU)*(MH));
	Pr = 1.0/3.0*(P_A)*pow(y[2], 4.0);

	return Pg+Pr;
}

double v_wind(double r)
{
	static double r_c[1000], v_c[1000];
	static int flag = 0, nsize;

	int i = 0;
	double v;
	double dummy[5];

	if(flag == 0){
		FILE *fp;
		char filename[256];
		fp = fopen(csm, "r");
		fgets(filename, sizeof(filename), fp);
		while(fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &dummy[0], &dummy[1], &r_c[i], &v_c[i], &dummy[2], &dummy[3], &dummy[4]) != EOF){
			i++;
		}
		nsize = i;
		flag = 1;
		fclose(fp);
	}
	if(r <= r_c[0]){
		v = v_c[0];
	}
	else{
		for(i = 0; i < nsize-1; i++){
			if(r < r_c[i+1] && r >= r_c[i]){
				break;
			}
		}
		v = (v_c[i+1]-v_c[i])/(log(r_c[i+1])-log(r_c[i]))*(log(r)-log(r_c[i]))+v_c[i];
	}

	return v;
}

void interp_self_similar_values(double *A, double *E_rev, double *E_for)
{
	FILE *fp;
	fp = fopen("LCFiles/interp_self_similar_values.txt", "r");
	fscanf(fp, "%lf %lf %lf", A, E_rev, E_for);
	fclose(fp);
}

void set_abundance(void)
{
	double r_c[2], rho_c[2];
	double dammy[5];
	FILE *fp;
	char filename[256];

	fp = fopen(csm, "r");
	fgets(filename, sizeof(filename), fp);
	fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &dammy[0], &dammy[1], &r_c[0], &dammy[2], &rho_c[0], &dammy[3], &dammy[4]);
	X = dammy[3];
	Y = dammy[4];
//	gen_opacity_sc();
	fclose(fp);
}

void gen_opacity_sc(void)
{
	FILE *fp, *fp1;
	fp = fopen("./LCFiles/opacity_sc.txt", "w");
	fp1 = fopen("./LCFiles/mean_ml_wght.txt", "w");
	double R, T;
	double op, mu;
	double x;
	int i = 0;

	R = 1.0e-008;
	T = pow(10., 3.20);
	for(x = -12.00; i < 22; x += 0.5){
		fprintf(fp, "%1.1f ", x);
		fprintf(fp1, "%1.1f ", x);
		i++;
	}
	fprintf(fp, "\n");
	fprintf(fp1, "\n");

	while(T < 1.0e+06){
		i = 0;
		fprintf(fp, "%1.2f", log10(T));
		fprintf(fp1, "%1.2f", log10(T));
		for(x = -12.0; i < 22; x += 0.5){
			R = pow(10., x);
			op = sigma_saha(R, T);
			sigma_mu_saha(R, T, &op, &mu);
			fprintf(fp, " %2.3f", log10(op));
			fprintf(fp1, " %1.5f", mu);
			i++;
		}
		fprintf(fp, "\n");
		fprintf(fp1, "\n");
		T *= pow(10.0, 0.050000);
	}

	fclose(fp);
	fclose(fp1);
}

double sigma_saha(double R, double T)
{
	double rho = R*pow(T*1.e-06, 3.);
	double mu_tmp[2] = {};
	double x;
	double n_e, n_H, n_He, n_HI, n_HII, n_HeI, n_HeII, n_HeIII;
	mu_tmp[0] = 0.5;
	mu_tmp[1] = 1.;
	while(fabs(mu_tmp[1]-mu_tmp[0]) > 1.e-15){
		mu_tmp[0] = mu_tmp[1];
		x = 2.*M_PI*P_E*P_K*T/(P_H*P_H);
		n_H = X*rho/MH;
		n_He = Y/4.*rho/MH;
		n_e = rho/(mu_tmp[0]*MH)-(X+Y/4.)*rho/MH;
		n_HII = pow(x, 1.5)*exp(-CHI_HI/(P_K*T))/(n_e+pow(x, 1.5)*exp(-CHI_HI/(P_K*T)))*n_H;
		n_HI = n_H-n_HII;
		n_HeI = pow(1.+4./n_e*pow(x, 1.5)*exp(-CHI_HeI/(P_K*T))+4./(n_e*n_e)*pow(x, 3.)*exp(-(CHI_HeI+CHI_HeII)/(P_K*T)), -1.)*n_He;
		n_HeII = 4.*n_HeI/n_e*pow(x, 1.5)*exp(-CHI_HeI/(P_K*T));
		n_HeIII = n_HeII/n_e*pow(x, 1.5)*exp(-CHI_HeII/(P_K*T));
		mu_tmp[1] = pow(X*(1.+n_HII/n_H)+Y/4.*(1.+n_HeII/n_He+2.*n_HeIII/n_He), -1.);
		mu_tmp[1] = (mu_tmp[1]+mu_tmp[0])/2.;
		mu_tmp[1] = (mu_tmp[1]+mu_tmp[0])/2.;
	}
	return n_e*SIGMA_TH/rho;
}

void sigma_mu_saha(double R, double T, double *sigma, double *mu)
{
	double rho = R*pow(T*1.e-06, 3.);
	double mu_tmp[2] = {};
	double x;
	double n_e, n_H, n_He, n_HI, n_HII, n_HeI, n_HeII, n_HeIII;
	mu_tmp[0] = 0.5;
	mu_tmp[1] = 1.;
	while(fabs(mu_tmp[1]-mu_tmp[0]) > 1.e-15){
		mu_tmp[0] = mu_tmp[1];
		x = 2.*M_PI*P_E*P_K*T/(P_H*P_H);
		n_H = X*rho/MH;
		n_He = Y/4.*rho/MH;
		n_e = rho/(mu_tmp[0]*MH)-(X+Y/4.)*rho/MH;
		n_HII = pow(x, 1.5)*exp(-CHI_HI/(P_K*T))/(n_e+pow(x, 1.5)*exp(-CHI_HI/(P_K*T)))*n_H;
		n_HI = n_H-n_HII;
		n_HeI = pow(1.+4./n_e*pow(x, 1.5)*exp(-CHI_HeI/(P_K*T))+4./(n_e*n_e)*pow(x, 3.)*exp(-(CHI_HeI+CHI_HeII)/(P_K*T)), -1.)*n_He;
		n_HeII = 4.*n_HeI/n_e*pow(x, 1.5)*exp(-CHI_HeI/(P_K*T));
		n_HeIII = n_HeII/n_e*pow(x, 1.5)*exp(-CHI_HeII/(P_K*T));
		mu_tmp[1] = pow(X*(1.+n_HII/n_H)+Y/4.*(1.+n_HeII/n_He+2.*n_HeIII/n_He), -1.);
		mu_tmp[1] = (mu_tmp[1]+mu_tmp[0])/2.;
		mu_tmp[1] = (mu_tmp[1]+mu_tmp[0])/2.;
	}
	*sigma = n_e*SIGMA_TH/rho;
	*mu = mu_tmp[1];
}

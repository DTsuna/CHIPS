#include <stdio.h>
#include <math.h>
#include "opacity.h"
#include "function.h"
#include "constant.h"
#include "pars.h"

char csm[256];
extern pars pdt;

double r_early(double t)
{
	double A, B, C;
	double n, s, delta, M_ej, E_ej, D;
	double r_out = 1.5e+14, r_in = 9.0e+13;
	double rho_out = rho_csm(r_out), rho_in = rho_csm(r_in);

	n = pdt.n;
	delta = pdt.delta;
	M_ej = pdt.M_ej;
	E_ej = pdt.E_ej;

	s = -(log(rho_out)-log(rho_in))/(log(r_out)-log(r_in));
	s = 1.5;
	D = rho_out*pow(r_out, s);
//	printf("%f %e\n", s, D);

	A = pow(2.*(5.-delta)*(n-5.)*E_ej, (n-3.)/2.);
	B = pow((3.-delta)*(n-3.)*M_ej, (n-5.)/2.);
	C = pow((3.-s)*(4.-s)/(4.*M_PI*D*(n-4.)*(n-3.)*(n-delta))*(A/B), 1./(n-s));
	return C*pow(t, (n-3.)/(n-s));
}

double t_early(double r)
{
	double A, B, C;
	double n, s, delta, M_ej, E_ej, D;
	double r_out = 1.5e+14, r_in = 9.0e+13;
	double rho_out = rho_csm(r_out), rho_in = rho_csm(r_in);

	n = pdt.n;
	delta = pdt.delta;
	M_ej = pdt.M_ej;
	E_ej = pdt.E_ej;

	s = -(log(rho_out)-log(rho_in))/(log(r_out)-log(r_in));
	s = 1.5;
	D = rho_out*pow(r_out, s);

	A = pow(2.*(5.-delta)*(n-5.)*E_ej, (n-3.)/2.);
	B = pow((3.-delta)*(n-3.)*M_ej, (n-5.)/2.);
	C = pow((3.-s)*(4.-s)/(4.*M_PI*D*(n-4.)*(n-3.)*(n-delta))*(A/B), 1./(n-s));
	return pow(r/C, (n-s)/(n-3.));
}

double rho_csm(double r)
{
	static double r_c[1000], rho_c[1000], s;
	static int flag = 0, nsize;

	int i = 0;
	double rho;
	double dammy[5];

	if(flag == 0){
		FILE *fp;
		char filename[256];
		fp = fopen(csm, "r");
		fgets(filename, 512, fp);
		while(fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &dammy[0], &dammy[1], &r_c[i], &dammy[2], &rho_c[i], &dammy[3], &dammy[4]) != EOF){
//			rho_c[i]/=3.15e+07;
			i++;
		}
		nsize = i;
		flag = 1;
		fclose(fp);
		s = (-(log10(rho_c[1])-log10(rho_c[0]))/(log10(r_c[1])-log10(r_c[0]))-(log10(rho_c[2])-log10(rho_c[1]))/(log10(r_c[2])-log10(r_c[1])))*0.5;
		s = 1.5;
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
	double dammy[7];
	fp = fopen(file_csm, "r");
	fgets(filename, 512, fp);
	fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &dammy[0], &dammy[1], &dammy[2], &dammy[3], &dammy[4], &dammy[5], &dammy[6]);

	fclose(fp);
	return dammy[2];
}

double set_r_diff(const char *file_csm)
{
	FILE *fp;
	char filename[256];
	double dammy[7];
	double n, s, delta, M_ej, E_ej, kappa, q;
	double A, g_to_n, rdiff;
	fp = fopen(file_csm, "r");
	fgets(filename, 512, fp);
	fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &dammy[0], &dammy[1], &dammy[2], &dammy[3], &dammy[4], &dammy[5], &dammy[6]);

	n = pdt.n;
	delta = pdt.delta;
	M_ej = pdt.M_ej;
	E_ej = pdt.E_ej;
	kappa = 0.2 * (1.+dammy[5]);
	A = interp_A(n);
	s = 1.5;
	q = dammy[4] * pow(dammy[2], s);
	g_to_n = pow(2.0*(5.0-delta)*(n-5.0)*E_ej, (n-3.0)/2.0)/pow((3.0-delta)*(n-3.0)*M_ej, (n-5.0)/2.0)/((n-delta)*4.*M_PI);
//	rdiff = pow(A * g_to_n * pow(q,n-4.), 2./n) * pow(2.*(n-3.)*kappa/3./(n-1.5)/P_C, 2.*(n-3.)/n);
	rdiff = pow(A * g_to_n * pow(q,n-4.), 2./n) * pow(2.*(n-3.)*kappa/21./(n-1.5)/P_C, 2.*(n-3.)/n);
	fclose(fp);
	return rdiff;
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
		fgets(filename, 256, fp);
		printf("%s\n", filename);
		for(i = 0; i < 3; i++){
			fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &dammy[0], &dammy[1], &r_c[i], &dammy[2], &rho_c[i], &dammy[3], &dammy[4]);
		}
		s = (-(log10(rho_c[1])-log10(rho_c[0]))/(log10(r_c[1])-log10(r_c[0]))-(log10(rho_c[2])-log10(rho_c[1]))/(log10(r_c[2])-log10(r_c[1])))*0.5;
		s = 1.5;
//		printf("s = %f\n", s);
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
	static double r_c[1000], v_c[1000], s;
	static int flag = 0, nsize;

	int i = 0;
	double v;
	double dammy[5];

	if(flag == 0){
		FILE *fp;
		char filename[256];
		fp = fopen(csm, "r");
		fgets(filename, 512, fp);
		while(fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &dammy[0], &dammy[1], &r_c[i], &v_c[i], &dammy[2], &dammy[3], &dammy[4]) != EOF){
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

double interp_A(double n)
{
	double array_n[20], array_A[20];
	int i = 0, j;
	FILE *fp;

	fp = fopen("./src/LC/shock/rel_n_A_gam_1.3.txt", "r");
	while(fscanf(fp, "%lf %lf", array_n+i, array_A+i) != EOF){
		i++;
	}
	for(j = 0; j < i; j++){
		if(array_n[j] <= n && n <= array_n[j+1]){
			break;
		}
	}

	fclose(fp);

	return (array_A[j+1]-array_A[j])/(array_n[j+1]-array_n[j])*(n-array_n[j])+array_A[j];
}

void interp_int_e(double n, double *E_rev, double *E_for)
{
	double array_n[20], array_E_rev[20], array_E_for[20];
	int i = 0, j;
	FILE *fp;

	fp = fopen("./src/LC/transfer/rel_n_e_gam_1.3.txt", "r");
	while(fscanf(fp, "%lf %lf %lf", array_n+i, array_E_rev+i, array_E_for+i) != EOF){
		i++;
	}

	for(j = 0; j < i; j++){
		if(array_n[j] <= n && n <= array_n[j+1]){
			break;
		}
	}

	fclose(fp);

	*E_rev = (array_E_rev[j+1]-array_E_rev[j])/(array_n[j+1]-array_n[j])*(n-array_n[j])+array_E_rev[j];
	*E_for = (array_E_for[j+1]-array_E_for[j])/(array_n[j+1]-array_n[j])*(n-array_n[j])+array_E_for[j];
}

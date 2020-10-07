#include <stdio.h>
#include <math.h>
#include <string.h>
#include "boundary.h"
#include "solver.h"
#include "pars.h"
#include "constant.h"
#include "function.h"

extern pars pdt;
extern double t_exp;
extern char csm[256];

double *calc_init_dist(double, double, double, double, double, double, const char*);
double *calc_dist(double[], double, double, double, double, double, const char*);
void init_egn(double, double[]);
void forward_egn(double[], double*, double[], double);

void init_egn(double r_ini, double egn[])
{
	egn[0] = 0.9*r_ini/t_exp;
	egn[1] = 1.1*egn[0];
	egn[3] = 1.1*r_ini;
	egn[2] = rho_csm(egn[3])*egn[1]*egn[1]*egn[1]*0.100;
}

void forward_egn(double array[], double *r_ini, double egn[], double dt)
{
	egn[2] = array[5];
	egn[1] = array[2];
	egn[0] = array[1];
	*r_ini = array[3]+egn[0]*dt;
	egn[3] = array[4]+egn[1]*dt;
}

int main(void)
{
	double int_phys[4], egn[4];
	char file_csm[] = "./inp-data/CSM_1.5.txt";
	double *array;
	double t = 86400., dt = 8640.;
	int i;
	FILE *fp;
	fp = fopen("./outp-data/CSM_1.5.txt", "w");

	t_exp = 86400.;
	egn[0] = 1.e+9;
	egn[1] = 1.1e+9;
	egn[2] = 1.e+14;
	egn[3] = 1.1e+14;

	pdt = setpars(10., 1., 1.e+51, 10.*M_SUN, 1.e+07, 86400., 1.e+12);	

	array = calc_init_dist(pdt.E_ej, pdt.M_ej, pdt.n, pdt.delta, t, 1.e+14, file_csm);
//	printf("t = %f d, %e %e %e %e %e, L = %e erg/s\n", 
//		array[0]/86400., array[1], array[2], array[3], array[4], array[5], 4.0*M_PI*array[4]*array[4]*array[5]);
		printf("t = %f d, %e %e %e %e %e, L = %e erg/s, di = %e\n", 
			array[0]/86400., array[1], array[2], array[3], array[4], array[5], 
			4.0*M_PI*array[4]*array[4]*array[5], 2.*M_PI*array[3]*array[3]*rho_csm(array[3])*pow(array[1], 3.));
		fprintf(fp, "%f %e %e %e %e %e %e\n", 
			array[0]/86400., array[1], array[2], array[3], array[4], array[5], rho_csm(array[3]));
	for(i = 0; i < 1000; i++){
		array = calc_dist(array, pdt.E_ej, pdt.M_ej, pdt.n, pdt.delta, dt, file_csm);
		printf("t = %f d, %e %e %e %e %e, L = %e erg/s, di = %e\n", 
			array[0]/86400., array[1], array[2], array[3], array[4], array[5], 
			4.0*M_PI*array[4]*array[4]*array[5], 2.*M_PI*array[3]*array[3]*rho_csm(array[3])*pow(array[1], 3.));
		fprintf(fp, "%f %e %e %e %e %e %e\n", 
			array[0]/86400., array[1], array[2], array[3], array[4], array[5], rho_csm(array[3]));
	}
	fclose(fp);
	return 0;
}

double *calc_init_dist(double E_ej, double M_ej, double n, double delta, double t_ini, double r_ini, const char *file_csm)
{
	const int nsize = 4;
	int i, j;
	double dtau = 0.50;
	double err, tol = 1.00e-09;
	double egn[4], degn[4], egnfd[4], p[4] = {};
	double phys[8], physfd[4];
	double J[16];
	static double outp_egn[6];

	t_exp = t_ini;
	strcpy(csm, file_csm);
	pdt = setpars(n, delta, E_ej, M_ej, 1.00e+07, t_exp, 0.00);
	init_egn(r_ini, egn);

	for(i = 0; i < 4; i++){
		egnfd[i] = fabs(egn[i]);
		degn[i] = 1.00e-04*egnfd[i];
	}

	physfd[0] = egn[1]; physfd[1] = rho_csm(egn[3])*egn[1]*egn[1];
	physfd[2] = egn[2]; physfd[3] = egn[3]-r_ini;

	do{
		err = 0.00;
		solver(r_ini, phys, egn, 0);
		for(i = 0; i < nsize; i++){
			egn[i] += degn[i];
			solver(r_ini, phys+nsize, egn, 0);
			for(j = 0; j < nsize; j++){
				J[nsize*j+i] = (phys[nsize+j]-phys[j])/degn[i]*egnfd[i]/physfd[j];
			}
			egn[i] -= degn[i];
		}
		for(i = 0; i < nsize; i++){
			egn[i] /= egnfd[i];
			phys[i] /= physfd[i];
		}
		get_itr_x(egn, p, J, phys, dtau, nsize);
		for(i = 0; i < nsize; i++){
			err += phys[i]*phys[i];
			egn[i] *= egnfd[i];
		}
		err = sqrt(err/(double)nsize);
	}while(err > tol);

	outp_egn[0] = t_ini; outp_egn[1] = egn[0]; outp_egn[2] = egn[1];
	outp_egn[3] = r_ini; outp_egn[4] = egn[3]; outp_egn[5] = egn[2];

	return outp_egn;
}

/*
array[0] = egn[0] = u_rs, array[1] = egn[1] = u_fs, array[2] = r_rs, array[3] = r_fs, array[4] = F_fs.
These values are old, so before calculating the distribution, r_rs, r_fs must be updated, i.e.
r_rs = r_rs + u_rs*dt, r_fs = r_fs + r_fs*dt.
*/

double *calc_dist(double array[], double E_ej, double M_ej, double n, double delta, double dt, const char *file_csm)
{
	const int nsize = 3;
	int i, j;
	double dtau = 0.50;
	double err, tol = 1.00e-08;
	double r_ini;
	double egn[4], degn[nsize], egnfd[4], p[4] = {};
	double phys[2*4], physfd[4];
	double J[16];
	static double outp_egn[6];

	t_exp = array[0]+dt;
	strcpy(csm, file_csm);
	pdt = setpars(n, delta, E_ej, M_ej, 1.00e+07, t_exp, 0.00);

	forward_egn(array, &r_ini, egn, dt);

	for(i = 0; i < nsize; i++){
		egnfd[i] = fabs(egn[i]);
		degn[i] = 1.00e-04*egnfd[i];
	}

	physfd[0] = egn[1]; physfd[1] = rho_csm(egn[3])*egn[1]*egn[1];
	physfd[2] = egn[2];

	do{
		err = 0.00;
		solver(r_ini, phys, egn, 1);
		for(i = 0; i < nsize; i++){
			egn[i] += degn[i];
			solver(r_ini, phys+nsize, egn, 1);
			for(j = 0; j < nsize; j++){
				J[nsize*j+i] = (phys[nsize+j]-phys[j])/degn[i]*egnfd[i]/physfd[j];
			}
			egn[i] -= degn[i];
		}
		for(i = 0; i < nsize; i++){
			egn[i] /= egnfd[i];
			phys[i] /= physfd[i];
		}
		get_itr_x(egn, p, J, phys, dtau, nsize);
		for(i = 0; i < nsize; i++){
			err += phys[i]*phys[i];
			egn[i] *= egnfd[i];
		}
		err = sqrt(err/(double)nsize);
	}while(err > tol);

	outp_egn[0] = t_exp; outp_egn[1] = egn[0]; outp_egn[2] = egn[1];
	outp_egn[3] = r_ini; outp_egn[4] = egn[3]; outp_egn[5] = egn[2];

	return outp_egn;
}

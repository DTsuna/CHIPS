#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "boundary.h"
#include "solver.h"
#include "pars.h"
#include "constant.h"
#include "function.h"
#include "W4.h"
#include "transfer.h"

extern pars pdt;
extern double t_exp;
extern char csm[256];

double *calc_init_dist(double, double, double, double, double, double, double, const char*);
double *calc_dist(double[], double, double, double, double, double, double, const char*, int*);
void init_egn(double, double[]);
void forward_egn(double[], double*, double[], double);
void shock_csm(double, double, double, double, double, double, const char*, const char*);

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



/*
array[0] = t_exp.
array[1] = egn[0] = u_rs, array[2] = egn[1] = u_fs, array[3] = r_rs, array[4] = r_fs, array[5] = F_fs.
*/

void shock_csm(double E_ej, double M_ej, double M_ni, double n, double delta, double s, const char *file_csm, const char *file_outp)
{
	double *array;
	double dt = 8640.;
	double t_ini, t_fin = 400.*86400.;
	double r_ini, r_ini_diff;
	int info;
	int flag = 1;
	int c = 0;
	FILE *fp;
	double egn[4], y[4]; //These arrays are used to get E_fs.


	strcpy(csm, file_csm);
	set_abundance();

	pdt = setpars(n, delta, E_ej, M_ej, M_ni, s, 1.e+07, 1.e+10);

	while(flag == 1){
		fp = fopen(file_outp, "w");
		fprintf(fp, "# day  v_rs  v_fs  r_rs  r_fs  E_(rad,fs)  rho_csm(r=r_fs)\n");
		r_ini = set_r_ini(file_csm);
		r_ini_diff = set_r_diff(file_csm);
		r_ini = fmax(r_ini, r_ini_diff);
		t_ini = t_early(r_ini)+(double)c*0.01*86400.;
		printf("t_ini = %e s\n", t_ini);
		c++;


/**************************1st step**************************/
		array = calc_init_dist(pdt.E_ej, pdt.M_ej, pdt.M_ni, pdt.n, pdt.delta, t_ini, r_ini, file_csm);
		egn[0] = array[1]; egn[1] = array[2]; egn[2] = array[5]; egn[3] = array[4];
		boundary(array[4], y, egn, 1, &info);

		printf("t = %f d, u_rs = %e cm/s, u_fs = %e cm/s, r_rs = %e cm, r_fs = %e cm, F_fs = %e erg/cm^2/s, L = %e erg/s, di = %e\n", 
			array[0]/86400., array[1], array[2], array[3], array[4], array[5], 
			4.0*M_PI*array[4]*array[4]*array[5], 2.*M_PI*array[3]*array[3]*rho_csm(array[4])*pow(array[1], 3.));
		fprintf(fp, "%f %e %e %e %e %e %e %e\n", 
			array[0]/86400., array[1], array[2], array[3], array[4], array[5], (P_A)*pow(y[2], 4.), rho_csm(array[4]));
/************************************************************/




/*************************************************************
Continue calculation when the temperature behind the forward shock drops to ~6,000 K or luminosity < 1e+40 erg/s.
*************************************************************/
		do{
			dt = 0.01*t_exp;
			if(t_exp+dt > t_fin){
				dt = t_fin-t_exp;
			}
			array = calc_dist(array, pdt.E_ej, pdt.M_ej, pdt.M_ni, pdt.n, pdt.delta, dt, file_csm, &info);
			egn[0] = array[1]; egn[1] = array[2]; egn[2] = array[5]; egn[3] = array[4];
			boundary(array[4], y, egn, 1, &info);
			if(info == 1){
				printf("t_ini += 0.01 d.\n");
				break;
			}
			if(4.*M_PI*array[4]*array[4]*array[5] < 1.e+40 || y[2] < pow(10., 3.8)){
				flag = 0;
				break;
			}

			printf("t = %f d, u_rs = %e cm/s, u_fs = %e cm/s, r_rs = %e cm, r_fs = %e cm, F_fs = %e erg/cm^2/s, L = %e erg/s, di = %e\n", 
				array[0]/86400., array[1], array[2], array[3], array[4], array[5], 
				4.0*M_PI*array[4]*array[4]*array[5], 2.*M_PI*array[3]*array[3]*rho_csm(array[3])*pow(array[1], 3.));
			fprintf(fp, "%f %e %e %e %e %e %e %e\n", 
				array[0]/86400., array[1], array[2], array[3], array[4], array[5], (P_A)*pow(y[2], 4.), rho_csm(array[3]));
		}while(t_exp < t_fin);

		if(t_exp-t_fin > -1.e-10*t_fin && t_exp-t_fin < 1.e-10*t_fin){
			flag = 0;
		}
		if(c > 10){
			flag = 0;
		}
		fclose(fp);
	}
}

double *calc_init_dist(double E_ej, double M_ej, double M_ni, double n, double delta, double t_ini, double r_ini, const char *file_csm)
{
	const int nsize = 4;
	int i, j, info = 0;
	double dtau = 0.50;
	double err, tol = 1.00e-09;
	double egn[4], degn[4], egnfd[4], p[4] = {};
	double phys[8], physfd[4];
	double J[16];
	static double outp_egn[6];

	t_exp = t_ini;
	init_egn(r_ini, egn);

	for(i = 0; i < 4; i++){
		egnfd[i] = fabs(egn[i]);
		degn[i] = 1.00e-05*egnfd[i];
	}

	physfd[0] = egn[1]; physfd[1] = rho_csm(egn[3])*egn[1]*egn[1];
	physfd[2] = egn[2]; physfd[3] = egn[3]-r_ini;

	do{
		err = 0.00;
		solver(r_ini, phys, egn, 0, &info);
		for(i = 0; i < nsize; i++){
			egn[i] += degn[i];
			solver(r_ini, phys+nsize, egn, 0, &info);
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
array[0] = t_exp, array[1] = egn[0] = u_rs, array[2] = egn[1] = u_fs, array[3] = r_rs, array[4] = r_fs, array[5] = F_fs.
These values are old, so before calculating the distribution, r_rs, r_fs must be updated, i.e.
r_rs = r_rs + u_rs*dt, r_fs = r_fs + r_fs*dt.
*/

double *calc_dist(double array[], double E_ej, double M_ej, double M_ni, double n, double delta, double dt, const char *file_csm, int *info)
{
	const int nsize = 3, cmax = 100;
	int i, j, c = 0;
	double dtau = 0.50;
	double err, tol = 1.00e-08;
	double r_ini;
	double egn_old[4], egn[4], degn[nsize], egnfd[4], p[4] = {};
	double phys[2*4], physfd[4];
	double J[16];
	static double outp_egn[6];

	memcpy(egn_old, egn, sizeof(double)*4);

	strcpy(csm, file_csm);

	*info = 0;

	while(1){
		c = 0;
		for(i = 0; i < 4; i++){
			p[i] = 0.;
//			egn[i] = (i==2) ? 0.9*egn_old[i] : egn_old[i];
			egn[i] = egn_old[i];
		}
		t_exp = array[0]+dt;
		forward_egn(array, &r_ini, egn, dt);
	
		for(i = 0; i < nsize; i++){
			egnfd[i] = fabs(egn[i]);
			degn[i] = 1.00e-06*egnfd[i];
		}
	
		physfd[0] = egn[1]; physfd[1] = rho_csm(egn[3])*egn[1]*egn[1];
		physfd[2] = egn[2];
	
		do{
			err = 0.00;
			solver(r_ini, phys, egn, 1, info);
			if(*info == 1){
				break;
			}
			for(i = 0; i < nsize; i++){
				egn[i] += degn[i];
				solver(r_ini, phys+nsize, egn, 1, info);
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
			c++;
			if(c == cmax){
				printf("count max (iteration failure)\n");
				break;
			}
		}while(err > tol);
		if(isnan(egn[0]) || isnan(egn[1]) || isnan(egn[2]) || isnan(egn[3]) || c == cmax){
			dt *= 0.1;
		}
		else{
			break;
		}
/*Calculation failure (NaN output)*/
		if(*info == 1){
			break;
		}
		if(dt < 0.01){
			break;
		}
	}

	if(dt < 0.01){
		*info = 1;
	}
	
	outp_egn[0] = t_exp; outp_egn[1] = egn[0]; outp_egn[2] = egn[1];
	outp_egn[3] = r_ini; outp_egn[4] = egn[3]; outp_egn[5] = egn[2];

	return outp_egn;
}

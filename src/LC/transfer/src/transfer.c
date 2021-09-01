#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "function.h"
#include "constant.h"
#include "itgadvu.h"
#include "itgadve.h"
#include "srctrm.h"
#include "rhocsm.h"
#include "raytrace.h"
#include "flux.h"
#include "opacity.h"
#include "saha.h"
#include "pars.h"

#ifdef _OPENMP
	#include <omp.h>
#endif

extern char csm[256];

void rad_transfer_csm(double, double, double, double, double, const char*, const char*, const char*, const char*, const char*);
void init_E_U(double, double, double[], double[], double[], double[], double[], const int);
void read_shockprofiles(FILE*, double[], double[], double[], int*);


void rad_transfer_csm(double Eexp, double Mej, double nej, double delta, double r_out, 
	const char *file_csm, const char *file_inp, const char *file_outp, const char *file_outp_band, const char *dir_shockprofiles)
{
	FILE *fp, *fl, *fw, *fnu_time, *fsh;
	double F_max = 0., F_out = 0.;
	double r_eff, T_eff, r_eff_interp, T_color;
	double E[2*NSIZE], U[2*NSIZE], r[NSIZE+1], E_old[NSIZE], rho[NSIZE], v_w[NSIZE], E0[NSIZE], U0[NSIZE];
	double T_g[NSIZE], mu[NSIZE], F[NSIZE];
	double kappa_s, kappa_a;
	double r_ini, F_ini, u_ini, E_ini;
	double t, dt = 4.;
	double err = 0., tol = 1.e-06;
	double rho_ed[2];
	double tf[20000], rf[20000], Ff[20000], Ef[20000], uf[20000];
	int i = 0, ii = 0, j = 0, k, l, n = NSIZE, fsize, flag = 0;
	int F_neg_flag = 0;
	int count = 0;
	int c = 0, cmax = 100;
	int L_outp_flag = 0, count_nu = 0;
	int FLNU;
	double dummy[8];
	double dr;
	double CFL = 1.0000000000000000;
	double tau[NSIZE], tau_eff[NSIZE];
	double tau_tot, tau_eff_tot;
	char filename[1024], profiles[1024];
	double g, A, q;
	double gam = 4./3.;
	double E_rev, E_for;
	double t_diff;
	double E_to_T;
	double F_mean;
	double r_sh[2000], T_sh[2000], rho_sh[2000];
	double abmag[5];
	char *cat_str = "mag_";
	int n_sh;
	int interval = 50;
	int outp_date_int = 0, outp_date_min;
	pars pdt;

	int num_of_threads, i_threads;



	pdt = setpars(nej, delta, Eexp, Mej, 1.e+7, 1.e+12);

	/* OpenMP parameters */
	if(getenv("OMP_NUM_THREADS") == NULL){
		num_of_threads = 1;
	}
	else if(atoi(getenv("OMP_NUM_THREADS")) == 0){
		num_of_threads = 1;
	}
	else{
		num_of_threads = atoi(getenv("OMP_NUM_THREADS"));
	}

	/* multi-band flag based on whether the multi-band output file string is empty or not*/
	if(strlen(file_outp_band) == 0) {
		FLNU = 0;
	} else {
		FLNU = 1;
	}

#ifdef _OPENMP
	printf("The number of used threads = %d.\n", num_of_threads);
#endif

	sprintf(csm, "%s", file_csm);
	sprintf(filename, "%s", file_inp);
	fp = fopen(filename, "r");
	set_abundance();

	if(FLNU == 1) {
		sprintf(filename, "%s", file_outp_band);
		fnu_time = fopen(filename, "w");
		fprintf(fnu_time, "%s %s %s %s %s %s\n", "time [day]", "U", "B", "V", "R", "I");
	}

	sprintf(filename, "%s", file_outp);
	fl = fopen(filename, "w");
	while(fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf", &dummy[0], &dummy[1], &dummy[2], &dummy[3], &dummy[4], &dummy[5], &dummy[6], &dummy[7]) != EOF){
		tf[i] = dummy[0]*86400.000;
		rf[i] = dummy[4];
		uf[i] = dummy[2];
		Ef[i] = dummy[6];
		if (F_neg_flag==0 && dummy[5]>0.0) {
			Ff[i] = dummy[5];
		} else {
			F_neg_flag = 1;
			Ff[i] = 0.0;
		}
		i++;
	}
	fsize = i;

	for(outp_date_int = 0; outp_date_int < 100; outp_date_int++){
		if(tf[0]/86400. >= (double)outp_date_int && tf[0]/86400. < (double)(outp_date_int+1)){
			break;
		}
	}
	outp_date_min = outp_date_int+1;
	outp_date_int++;

	double last_dt = tf[fsize-1] - tf[fsize-2];
	double slope_dlogr_dlogt = (tf[fsize-1]/rf[fsize-1]) * (rf[fsize-1]-rf[fsize-11])/(tf[fsize-1]-tf[fsize-11]);
	double slope_dlogu_dlogt = (tf[fsize-1]/rf[fsize-1]) * (uf[fsize-1]-uf[fsize-11])/(tf[fsize-1]-tf[fsize-11]);
	double slope_dlogF_dlogt = (tf[fsize-1]/Ff[fsize-1]) * (Ff[fsize-1]-Ff[fsize-11])/(tf[fsize-1]-tf[fsize-11]);
	double slope_dlogE_dlogt = (tf[fsize-1]/Ef[fsize-1]) * (Ef[fsize-1]-Ef[fsize-11])/(tf[fsize-1]-tf[fsize-11]);
	for(i = fsize; i < 2*fsize; i++){
		tf[i] = tf[fsize-1] + (double) (i-fsize+1) * last_dt;
		rf[i] = rf[fsize-1] * pow(tf[i]/tf[fsize-1], slope_dlogr_dlogt);
		Ff[i] = Ff[fsize-1] * pow(tf[i]/tf[fsize-1], slope_dlogF_dlogt);
		Ef[i] = Ef[fsize-1] * pow(tf[i]/tf[fsize-1], slope_dlogE_dlogt);
		uf[i] = Ef[fsize-1] * pow(tf[i]/tf[fsize-1], slope_dlogu_dlogt);
	}

	t = tf[0];
	r_ini = rf[0];
	F_ini = Ff[0];
	u_ini = uf[0];
	E_ini = 0.;



/********************Calculate internal energy deposited in shocked region********************/
	A = interp_A(nej);
	interp_int_e(nej, &E_rev, &E_for);
	g = 1./(4.*M_PI*(nej-delta))*pow(2.*(5.-delta)*(nej-5.), (nej-3.)/2.)/pow((3.-delta)*(nej-3.), (nej-5.)/2.);
	g *= pow(Eexp/Mej, (nej-3.)/2.);
	g *= Mej;
//g -> g^nej
	q = rho_csm(rf[0])*pow(rf[0], 1.5);
	E_rev *= 4.*M_PI*(nej-3.)/((gam-1.)*(nej-1.5))*g*pow(A*g/q, (5.-nej)/(nej-1.5))*pow(tf[0], 1.5*(nej-5.)/(nej-1.5));
	E_for *= 4.*M_PI*(nej-3.)/((gam-1.)*(nej-1.5))*q*pow(A*g/q, 3.5/(nej-1.5))*pow(tf[0], 1.5*(nej-5.)/(nej-1.5));
	t_diff = rf[0]/uf[0];
/*********************************************************************************************/




	init_E_U(r_ini, r_out, r, rho, v_w, E, U, NSIZE);
	dr = r[1]-r[0];

/*
If r_fs(t) > r[0]-dr/3, r[0] = r[1], r[1] = r[2], ... .
Here, dr = r[N]-r[N-1] must be fixed.
for i = 0, 1, ..., n-1, 
X[i] = X[i+1], where X is physical quantity, i.e. E, U, rho.
*/
	while(t < tf[2*fsize-1]){

/*
dt does not necesarrily satisfy CFL condition.
*/
		if(flag == 0){
			dt = CFL*dr/(P_C);
		}
		else{
			flag = 0;
		}

/*
Identify the position of forward shock, and estimate by linear interpolation.
*/
		for(i = j; i < 2*fsize; i++){
			if(t >= tf[i] && t < tf[i+1]){
				j = i;
			}
		}

		if(FLNU == 1 && t+dt > tf[j+1]){
			L_outp_flag = 1;
			dt = tf[j+1]-t;
		}

/*Interpolation of r, E, u_fs, F*/
		r_ini = rf[j]*exp(log(rf[j+1]/rf[j])/log(tf[j+1]/tf[j])*log((t+dt)/tf[j]));
		E_to_T = Ef[j]*exp(log(Ef[j+1]/Ef[j])/log(tf[j+1]/tf[j])*log((t+dt)/tf[j]));
		u_ini = uf[j]*exp(log(uf[j+1]/uf[j])/log(tf[j+1]/tf[j])*log((t+dt)/tf[j]));
		if(Ff[j] > 0. && Ff[j+1] > 0.){
			F_ini = Ff[j]*exp(log(Ff[j+1]/Ff[j])/log(tf[j+1]/tf[j])*log((t+dt)/tf[j]));
		}
		else{
			printf("Negative Flux value.\n");
			F_ini = 1.e+04;
		}

		r[0] = r_ini;
		if(r_ini > r[1]-dr/4.){
			E_ini = E[0];
			n--;
			for(i = 0; i < n; i++){
				if(i == 0){
					E[0] = E[0]*(r[1]*r[1]*r[1]-r[0]*r[0]*r[0])/(r[2]*r[2]*r[2]-r[0]*r[0]*r[0])
						+E[2]*(r[2]*r[2]*r[2]-r[1]*r[1]*r[1])/(r[2]*r[2]*r[2]-r[0]*r[0]*r[0]);
					U[0] = U[0]*(r[1]*r[1]*r[1]-r[0]*r[0]*r[0])/(r[2]*r[2]*r[2]-r[0]*r[0]*r[0])
						+U[2]*(r[2]*r[2]*r[2]-r[1]*r[1]*r[1])/(r[2]*r[2]*r[2]-r[0]*r[0]*r[0]);
				}
				else{
					E[2*i] = E[2*(i+1)];
					U[2*i] = U[2*(i+1)];
				}
				E[2*i+1] = E[2*(i+1)+1];
				U[2*i+1] = U[2*(i+1)+1];
				if(i != 0){
					r[i] = r[i+1];
				}
				rho[i] = rho[i+1];
				v_w[i] = v_w[i+1];
			}
			r[n] = r[n+1];
		}

		for(i = 0; i < n; i++){
			E0[i] = E[2*i];
			U0[i] = U[2*i];
//			printf("E = %e U = %e\n", E0[i], U0[i]);
		}

#ifdef EADD
		F_ini += E_ini*u_ini+(E_rev+E_for)/(4.*M_PI*r_ini*r_ini*t_diff) * exp(1.0 - t/t_diff);
#endif


/************************************************************************************************************
In the following, integrate radiative transfer equation and energy eqation using operator splitting.
i)   Integrate source term dE/dt = kappa*rho*(acT^4-E), dU/dt = -kappa*rho(acT^4-E), by implicit Euler method.
ii)  Integrate energy equation implicitly.
iii) Integrate 0th moment equation implicitly. Iteration is needed to complete the calculation.
************************************************************************************************************/


/************************************************************
At first, intergrate source term using implicit Euler method.
************************************************************/
//#pragma omp parallel for shared(E, U), private(i)
		for(i_threads = 0; i_threads < num_of_threads; i_threads++){
			for(i = i_threads; i < n; i += num_of_threads){
				itg_src(&E[2*i], &U[2*i], rho[i], dt, tol);
			}
		}

		for(i = 0; i < n; i++){
			if(isnan(E[2*i+1]) != 0){
				flag = 1;
				break;
			}
		}
		
		if(flag == 1){
			for(k = 0; k <= i; k++){
				E[2*k+1] = E[2*k];
				U[2*k+1] = U[2*k];
			}
			dt *= 0.5;
			fprintf(stderr, "iteration failure (source term). Too large time step.");
			continue;
		}
		else{
			for(i = 0; i < n; i++){
				E[2*i] = E[2*i+1];
				U[2*i] = U[2*i+1];
			}
			flag = 0;
		}

/*
Integrate energy equation impilicitly.
*/
		rho_ed[0] = rho_csm(r[0]-(r[1]-r[0])/2.);
		rho_ed[1] = rho_csm(r[n]+dr/2.);
		itg_adv_U(r, U, rho, rho_ed, v_w, dt, n);


	
/*
Integrate 0th moment equation implicitly.
Here, iteration is needed to complete the calculation.
*/



/*
E_old[n] must keep values of E[2*i+1] before iteration, so that error is estimated.
*/
		for(i = 0; i < n; i++){
			E_old[i] = E[2*i];
		}
		

		do{
			err = 0.;
			itg_adv_E(F_ini, r, E, U, rho, dt, n);
			for(i = 0; i < n; i++){
				if(isnan(E[2*i+1])){
					flag = 1;
					break;
				}
			}
			if(flag == 1){
				fprintf(stderr, "iteration failure. time step is too large.\n");
				break;
			}
			for(i = 0; i < n; i++){
				err += pow(1.-E[2*i+1]/E_old[i], 2.);
				E_old[i] = E[2*i+1];
			}
	
			err = sqrt(err/(double)n);
		}while(err > tol);

		if(flag == 1){
			for(i = 0; i < n; i++){
				E[2*i] = E0[i];
				E[2*i+1] = E[2*i];
				U[2*i] = U0[i];
				U[2*i+1] = U[2*i];
			}
			dt *= 0.5;
			c++;
			if(c == cmax){
				fprintf(stderr, "too small dt. exit.\n");
				exit(EXIT_FAILURE);
			}
			continue;
		}
		else{
			c = 0;
			flag = 0;
		}

		t += dt;

		for(i = 0; i < n; i++){
			E[2*i] = E[2*i+1];
			U[2*i] = U[2*i+1];
		}

		F_out = (P_C)*E[2*(n-1)];
		if(F_max < F_out){
			F_max = F_out;
		}

		if(F_out < F_max*0.1 && 4.*M_PI*r[n-1]*r[n-1]*F_out < 1.e+41){
			break;
		}


		for(i = 0; i < n; i++){
			saha(rho[i], U[2*i], mu+i, T_g+i);
			kappa_s = kappa_r(rho[i], T_g[i]);
			kappa_a = kappa_p(rho[i], T_g[i]);
			if(i != 0){
				F[i-1] = calc_flux_i(r_ini, r, E, U, rho, dt, i, n);
			}
			tau_eff[i] = sqrt(kappa_a*kappa_s)*rho[i]*(r[i+1]-r[i]);
			tau[i] = kappa_s*rho[i]*(r[i+1]-r[i]);
		}

		tau_tot = 0.;
		tau_eff_tot = 0.;
		for(i = n-1; i > 0; --i){
			tau_tot += tau[i];
			if(tau_tot < 1. && tau_tot+tau[i-1] > 1.){
				break;
			}
		}
		for(ii = n-1; ii > 0; --ii){
			tau_eff_tot += tau_eff[ii];
			if(tau_eff_tot < 1. && tau_eff_tot+tau_eff[ii-1] > 1.){
				break;
			}
		}


/********************Calculate color temperature*********************/
		if(ii != 1 && ii != 0){
			T_color = (T_g[ii-2]-T_g[ii])/2./tau_eff[ii-1]*(1.-tau_eff_tot)+(T_g[ii-1]+T_g[ii])/2.;
		}
		else{
			r_eff_interp = r_eff;
			F_mean = F_ini;
			T_color = pow(E_to_T/(P_A), 0.25);
			if(fabs(T_color) < 1.e-20){
				T_color = NAN;
			}
		}
/********************************************************************/

/***********************Calculate photosphere************************/
		if(i != 1 && i != 0){
			r_eff = (r[i-2]-r[i])/2./tau[i-1]*(1.-tau_tot)+(r[i-1]+r[i])/2.;
		}
		else{
			r_eff = r_ini;
		}
/********************************************************************/

/******************Calculate L_nu if requested************************/
		if(FLNU == 1 && L_outp_flag == 1){
			if(t/86400. > (double)outp_date_int && j < fsize-1){
				printf("/*********************************************************/\n");
				n_sh = 0;
				printf("%e\n", tf[j+1]);
				sprintf(filename, "%s/profiles%08d.txt", dir_shockprofiles, j+1);
				printf("Open shock profile: %s/profiles%08d.txt\n", dir_shockprofiles, j+1);
				fsh = fopen(filename, "r");
				read_shockprofiles(fsh, r_sh, rho_sh, T_sh, &n_sh);
				sprintf(filename, "%s/Lnu%08d", dir_shockprofiles, outp_date_int-outp_date_min);
				calc_lum(t, r[0], r_out, F_ini, r, rho, T_g, r_sh, rho_sh, T_sh, n, n_sh, filename, abmag, pdt);
				fprintf(fnu_time, "%e %e %e %e %e %e\n", tf[j+1]/86400., abmag[0], abmag[1], abmag[2], abmag[3], abmag[4]);
				printf("%e %e %e %e %e %e\n", tf[j+1]/86400., abmag[0], abmag[1], abmag[2], abmag[3], abmag[4]);
				outp_date_int++;
				fclose(fsh);
				printf("/*********************************************************/\n");
			}
			count_nu++;
			L_outp_flag = 0;
		}
/********************************************************************/




/**********************************************************************************
Output of temperature, radiation energy density, flux as functions of radius.
**********************************************************************************/
#if 0
		sprintf(filename, "LCFiles/dist%08d.txt", count);
		count++;
		fw = fopen(filename, "w");
		for(l = 0; l < n; l++){
			fprintf(fw, "%e %e %e %e %e\n", r[l], F[l], E[l], U[l], T_g[l]);
		}
		fclose(fw);
#endif
/*********************************************************************************/


	
		fprintf(fl, "%f %e %e %e\n", t/86400., 4.*M_PI*r[n-1]*r[n-1]*(P_C)*E[2*(n-1)], r_eff, T_color);
		printf("j = %d, n = %d, %f %e %e %e\n", j, n, t/86400., 4.*M_PI*r[n-1]*r[n-1]*(P_C)*E[2*(n-1)], r_eff, T_color);
	}
	fclose(fp);
	fclose(fl);
	if(FLNU == 1){
		fclose(fnu_time);
	}
}

void init_E_U(double r_ini, double r_out, double r[], double rho[], double v_w[], double E[], double U[], const int nsize)
{
	int i;
	double dr = (r_out-r_ini)/((double)NSIZE);
	double T;

	for(i = 0; i < nsize; i++){
		r[i] = r_ini+dr*(double)i;
		E[2*i] = (P_A)*pow(1.e+03, 4.)*pow(r[0]/r[i], 2.);
		T = pow(E[2*i]/(P_A), 0.25);
		E[2*i+1] = E[2*i];
		rho[i] = rho_csm(r[i]+dr/2.);
		v_w[i] = v_wind(r[i]);
		U[2*i] = 1.5*rho[i]*(P_K)*T/(1.3*(MH));
		U[2*i+1] = U[2*i];
	}
	r[nsize] = r[nsize-1]+dr;
}

void read_shockprofiles(FILE *f_sh, double r_sh[], double rho_sh[], double T_sh[], int *n_sh)
{
	int i = 0, j, positive = 0;
	double dummy[2];
	while(fscanf(f_sh, "%lf %lf %lf %lf %lf", r_sh+i, dummy, rho_sh+i, T_sh+i, dummy+1) != EOF){
		i++;
		if(dummy[0] > 0.){
			positive++;
		}
	}
	*n_sh = i-1;
	for(j = 0; j < *n_sh+1; j++){
	}

	for(j = 0; j < positive-1; j++){
		rho_sh[j] = (rho_sh[j]+rho_sh[j+1])/2.;
		T_sh[j] = (T_sh[j]+T_sh[j+1])/2.;
	}
	for(j = positive-1; j < i-1; j++){
		r_sh[j] = r_sh[j+1];
		rho_sh[j] = rho_sh[j+1];
		T_sh[j] = T_sh[j+1];
	}
	for(j = positive-1; j < i-2; j++){
		rho_sh[j] = (rho_sh[j]+rho_sh[j+1])/2.;
		T_sh[j] = (T_sh[j]+T_sh[j+1])/2.;
	}
	*n_sh = *n_sh-1;
}

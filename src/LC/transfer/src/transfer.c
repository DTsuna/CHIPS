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
#include "raytrace.h"
#include "flux.h"
#include "opacity.h"
#include "saha.h"
#include "pars.h"
#include "popov.h"

#ifdef _OPENMP
	#include <omp.h>
#endif

extern char csm[256];
extern double X;

void rad_transfer_csm(double, double, double, double, double, double, double, double, const char*, const char*, const char*, const char*, const char*);
void init_E_U(double, double, double[], double[], double[], double[], double[], const int);


void rad_transfer_csm(double Eej, double Mej, double Mni, double nej, double delta, double r_out, double out_interval, double plateau_flag,
	const char *file_csm, const char *file_inp, const char *file_outp, const char *file_outp_band, const char *dir_shockprofiles)
{
	FILE *fp, *fl, *fnu_time;
	double F_max = 0., F_out = 0.;
	double r_eff, v_eff, r_eff_interp, T_color;
	double E[2*NSIZE], U[2*NSIZE], r[NSIZE+1], E_old[NSIZE], rho[NSIZE], v_w[NSIZE], E0[NSIZE], U0[NSIZE];
	double T_g[NSIZE], mu[NSIZE], F[NSIZE];
	double kappa_s, kappa_a;
	double r_ini, F_ini, u_ini, E_ini;
	double t, dt = 4.;
	double err = 0., tol = 1.e-08;
	double rho_ed[2];
	double tf[20000], rf[20000], Ff[20000], Ef[20000], uf[20000];
	int count_e;
	int i = 0, ii = 0, j = 0, k, n = NSIZE, fsize, flag = 0, jjj;
	int F_neg_flag = 0;
	int output_flag = 0;
	int c = 0, cmax = 100;
	int FLNU;
	double dummy[8];
	double dr;
	double CFL = 0.9000000000000000;
	double tau[NSIZE], tau_eff[NSIZE];
	double tau_tot, tau_eff_tot;
	char buf[256], filename[1024];
	double g, A, q;
	double gam = 4./3.;
	double E_rev, E_for;
	double t_diff;
	double E_to_T;
	double F_mean;
	double F_from_shocked;
	double when_out[10000], time_interval = 86400.*out_interval;
	pars pdt;

	int num_of_threads, i_threads;


	FILE *fw;
	int l, count = 0;


	pdt = setpars(nej, delta, Eej, Mej, Mni, 1.e+7, 1.e+12);

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
	fgets(buf, 256, fp);
	set_abundance();


	if(FLNU == 1) {
		sprintf(filename, "%s", file_outp_band);
		fnu_time = fopen(filename, "w");
		fprintf(fnu_time, "%s %s %s %s %s %s\n", "time [day]", "U", "B", "V", "R", "I");
	}

	sprintf(filename, "%s", file_outp);
	fl = fopen(filename, "w");
	fprintf(fl, "# day  luminosity  r_eff  v_eff  T_color\n");
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

//************************ set output date**********************************
	when_out[0] = 86400.*(floor(tf[0]/86400.)+1.);
	for(jjj = 1; jjj < 10000; jjj++){
		when_out[jjj] = when_out[jjj-1]+time_interval;
	}
//**************************************************************************

	double last_dt = tf[fsize-1] - tf[fsize-2];
	double slope_dlogr_dlogt = (tf[fsize-1]/rf[fsize-1]) * (rf[fsize-1]-rf[fsize-11])/(tf[fsize-1]-tf[fsize-11]);
	double slope_dlogu_dlogt = (tf[fsize-1]/rf[fsize-1]) * (uf[fsize-1]-uf[fsize-11])/(tf[fsize-1]-tf[fsize-11]);
	double slope_dlogF_dlogt = (tf[fsize-1]/Ff[fsize-1]) * (Ff[fsize-1]-Ff[fsize-11])/(tf[fsize-1]-tf[fsize-11]);
	double slope_dlogE_dlogt = (tf[fsize-1]/Ef[fsize-1]) * (Ef[fsize-1]-Ef[fsize-11])/(tf[fsize-1]-tf[fsize-11]);
	for(i = fsize; i < 5*fsize; i++){
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
	g *= pow(pdt.E_ej/pdt.M_ej, (nej-3.)/2.);
	g *= pdt.M_ej;
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
	while(t < tf[5*fsize-1] + r_out/P_C){

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
		for(i = j; i < 5*fsize; i++){
			if(t >= tf[i] && t < tf[i+1]){
				j = i;
			}
		}

		if(FLNU == 1 && t+dt > tf[j+1]){
			dt = tf[j+1]-t;
		}

		if(t+dt > when_out[count]){
			dt = when_out[count]-t;
			output_flag = 1;
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
		F_from_shocked = (E_rev+E_for)/(4.*M_PI*r_ini*r_ini*t_diff) * exp(-t/t_diff - 0.5*t*t/t_diff/t_diff) * sqrt(2./M_E/M_PI) / (1.-erf(1./sqrt(2.)));
		F_ini += E_ini*u_ini+F_from_shocked;
		F_ini += plateau_flag*L_popov(t, pdt.M_ej, pdt.E_ej, 1.e+14, X)/(4.*M_PI*r_ini*r_ini);
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
			if(dt < 1.e-05 || isnan(dt)){
				fprintf(stderr, "Too small time interval. end.\n");
				exit(EXIT_FAILURE);
			}
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

		for(i = 0; i < n; i++){
			saha(rho[i], U[2*i+1], mu+i, T_g+i);
		}

	
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
		
		count_e = 0;
		do{
			count_e++;
			err = 0.;
			itg_adv_E(F_ini, r, E, U, rho, dt, n);
			for(i = 0; i < n; i++){
				if(isnan(E[2*i+1])){
					flag = 1;
					break;
				}
			}
			if(flag == 1){
				fprintf(stderr, "iteration failure. time step is too large. i = %d.\n", i);
				break;
			}
			for(i = 0; i < n; i++){
				err += pow(1.-E[2*i+1]/E_old[i], 2.);
				E_old[i] = E[2*i+1];
			}
			err = sqrt(err/(double)n);
			if(count_e == 300){
				printf("count max (rad tra iter), err = %e\n", err);
				break;
			}
	
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

		/* exit when luminosity is less than 1e41 ergs, or forward shock reached outer edge */
		if((F_out < F_max*0.1 && 4.*M_PI*r[n-1]*r[n-1]*F_out < 1.e+41) || r_ini > 0.99*r_out){
			break;
		}

		for(i = 0; i < n; i++){
			saha(rho[i], U[2*i], mu+i, T_g+i);
			kappa_s = kappa_r(rho[i], T_g[i]);
			kappa_a = kappa_p(rho[i], T_g[i]);
			if(i != 0){
				F[i] = calc_flux_i(r_ini, r, E, U, rho, dt, i, n);
				v_w[i] += dt*(kappa_r(rho[i], T_g[i])+kappa_r(rho[i-1], T_g[i-1]))/2.*F[i]/(P_C);
			}
			else{
				F[i] = F_ini;
				v_w[i] += dt*(kappa_r(rho[i], T_g[i]))*F[i]/(P_C);
			}

			tau_eff[i] = sqrt(kappa_a*kappa_s)*rho[i]*(r[i+1]-r[i]);
			tau[i] = kappa_s*rho[i]*(r[i+1]-r[i]);
		}


		for(i = 0; i < n; i++){
			saha(rho[i], U[2*i], mu+i, T_g+i);
			kappa_s = kappa_r(rho[i], T_g[i]);
			kappa_a = kappa_p(rho[i], T_g[i]);
			if(i != 0){
				F[i] = calc_flux_i(r_ini, r, E, U, rho, dt, i, n);
				v_w[i] += dt*(kappa_r(rho[i], T_g[i])+kappa_r(rho[i-1], T_g[i-1]))/2.*F[i]/(P_C);
			}
			else{
				F[i] = F_ini;
				v_w[i] += dt*(kappa_r(rho[i], T_g[i]))*F[i]/(P_C);
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
			v_eff = (v_w[i-2]-v_w[i])/2./tau[i-1]*(1.-tau_tot)+(v_w[i-1]+v_w[i])/2.;
		}
		else{
			r_eff = r_ini;
		}
/********************************************************************/

/**********************************************************************************
Output of temperature, radiation energy density, flux as functions of radius.
**********************************************************************************/
		if(output_flag == 1){
			sprintf(filename, "LCFiles/transfer_output/dist%05d.txt", count);
			count++;
			fw = fopen(filename, "w");
			fprintf(fw, "#r, F, E, U, T_g, rho, v_w\n");
			fprintf(fw, "#t = %f days.\n", t/86400.);
			fprintf(fw, "#Temp = %e K.\n", pow(4.*F_ini/(P_A)/(P_C), 0.25));
			fprintf(fw, "#input vel. = %e cm/s\n", u_ini*6./7.+v_wind(r_ini)/7.);
			for(l = 0; l < n; l++){
				fprintf(fw, "%e %e %e %e %e %e %e\n", r[l], F[l], E[2*l], U[2*l], T_g[l], rho[l], v_w[l]);
			}
			fclose(fw);
			printf("t = %f d, L = %e erg/s, r_eff = %e cm, T_col = %e K, F_ini = %e erg/s/cm2\n", t/86400., 4.*M_PI*r[n-1]*r[n-1]*(P_C)*E[2*(n-1)], r_eff, T_color, F_ini);
			output_flag = 0;
		}
/*********************************************************************************/


	
		fprintf(fl, "%f %e %e %e %e\n", t/86400., 4.*M_PI*r[n-1]*r[n-1]*(P_C)*E[2*(n-1)], r_eff, v_eff, T_color);
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

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
#include "flux.h"
#include "opacity.h"
#include "saha.h"

extern char csm[256];

void rad_transfer_csm(double, double, double, double, double, const char*, const char*, const char*);
void init_E_U(double, double, double[], double[], double[], double[], double[], const int);

//int main(void)
//{
//	double r_out = 9.9e+15;
//	char *file_csm = "./inp-data/CSM_1.5.txt";
//	char *file_inp = "./inp-data/CSM_1.5_profile.txt";
//
//
//	rad_transfer_csm(r_out, file_csm, file_inp);
//
//	return 0;
//}

void rad_transfer_csm(double Eexp, double Mej, double nej, double delta, double r_out, const char *file_csm, const char *file_inp, const char *file_outp)
{
	FILE *fp, *fl, *fw;
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
	double dummy[8];
	double dr;
	double CFL = 1.00000000000000;
	double tau[NSIZE], tau_eff[NSIZE];
	double tau_tot, tau_eff_tot;
	char filename[256];
	double g, A, q;
	double gam = 4./3.;
	double E_rev, E_for;
	double t_diff;
	double E_to_T;
	double F_mean;

	sprintf(csm, "%s", file_csm);
	sprintf(filename, "%s", file_inp);
	fp = fopen(filename, "r");


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
	double last_dt = tf[fsize-1] - tf[fsize-2];
	double slope_dlogr_dlogt = (tf[fsize-1]/rf[fsize-1]) * (rf[fsize-1]-rf[fsize-11])/(tf[fsize-1]-tf[fsize-11]);
	double slope_dlogF_dlogt = (tf[fsize-1]/Ff[fsize-1]) * (Ff[fsize-1]-Ff[fsize-11])/(tf[fsize-1]-tf[fsize-11]);
	double slope_dlogE_dlogt = (tf[fsize-1]/Ef[fsize-1]) * (Ef[fsize-1]-Ef[fsize-11])/(tf[fsize-1]-tf[fsize-11]);
	for(i = fsize; i < 2*fsize; i++){
		tf[i] = tf[fsize-1] + (double) (i-fsize+1) * last_dt;
		rf[i] = rf[fsize-1] * pow(tf[i]/tf[fsize-1], slope_dlogr_dlogt);
		Ff[i] = Ff[fsize-1] * pow(tf[i]/tf[fsize-1], slope_dlogF_dlogt);
		Ef[i] = Ef[fsize-1] * pow(tf[i]/tf[fsize-1], slope_dlogE_dlogt);
		uf[i] = 0.0;
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
		for(i = j; i < fsize; i++){
			if(t >= tf[i] && t < tf[i+1]){
				j = i;
			}
		}

/*Interpolation of r, E, u_fs, F*/
		r_ini = rf[j]*exp(log(rf[j+1]/rf[j])/log(tf[j+1]/tf[j])*log(t/tf[j]));
		E_to_T = Ef[j]*exp(log(Ef[j+1]/Ef[j])/log(tf[j+1]/tf[j])*log(t/tf[j]));
		u_ini = uf[j]*exp(log(uf[j+1]/uf[j])/log(tf[j+1]/tf[j])*log(t/tf[j]));
		if(Ff[j] > 0. && Ff[j+1] > 0.){
			F_ini = Ff[j]*exp(log(Ff[j+1]/Ff[j])/log(tf[j+1]/tf[j])*log(t/tf[j]));
		}
		else{
			printf("Negative Flux value.\n");
			F_ini = 1.e+04;
		}

		if(r_ini > r[0]-dr/4.){
			E_ini = E[0];
			n--;
			for(i = 0; i < n; i++){
				E[2*i] = E[2*(i+1)];
				U[2*i] = U[2*(i+1)];
				E[2*i+1] = E[2*(i+1)+1];
				U[2*i+1] = U[2*(i+1)+1];
				r[i] = r[i+1];
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
		//if(t < tf[0]+t_diff){
		//	F_ini += E_ini*u_ini+(E_rev+E_for)/(4.*M_PI*r_ini*r_ini*t_diff);
		//}
		//else{
		//	F_ini += E_ini*u_ini;
		//}
		F_ini += E_ini*u_ini+(E_rev+E_for)/(4.*M_PI*r_ini*r_ini*t_diff) * exp(1.0 - t/t_diff);

#endif


/*
In the following, integrate radiative transfer equation and energy eqation using operator splitting.
i)   Integrate source term dE/dt = kappa*rho*(acT^4-E), dU/dt = -kappa*rho(acT^4-E), by implicit Euler method.
ii)  Integrate energy equation implicitly.
iii) Integrate 0th moment equation implicitly. Iteration is needed to complete the calculation.
*/




/*
At first, intergrate source term using implicit Euler method.
*/
		for(i = 0; i < n; i++){
//			printf("i = %d r[%d] = %e cm rho[%d] = %e\n", i, i, r[i], i, rho[i]);
			itg_src(E+2*i, U+2*i, rho[i], dt, tol);
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
		rho_ed[0] = rho_csm(r[0]-(r[1]-r[0]));
		rho_ed[1] = rho_csm(r[n]);
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
			itg_adv_E(r_ini, F_ini, r, E, U, rho, dt, n);
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
//			printf("err = %e\n", err);
		}while(err > tol);

		if(flag == 1){
			for(i = 0; i < n; i++){
				E[2*i] = E0[i];
				E[2*i+1] = E[2*i];
				U[2*i] = U0[i];
				U[2*i+1] = U[2*i];
			}
			dt *= 0.5;
//			printf("back to (?)\n");
			continue;
		}
		else{
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
//			r_eff = (r[j-1]+r[j-2])/2.;
//			T_eff = pow(2.*(F[i-1]+F[i-2])/((P_A)*(P_C)), 0.25);
//			T_color = T_g[j-1];
//			F_mean = (F[j-1]+F[j-2])/2.;
//			F_mean = (F[j-2]-F[j-1])/tau[j-1]*(1.-tau_eff_tot)+F[j-1];
			T_color = (T_g[ii-2]-T_g[ii])/2./tau_eff[ii-1]*(1.-tau_eff_tot)+(T_g[ii-1]+T_g[ii])/2.;
//			r_eff_interp = (r[i-2]-r[i])/2./tau[i-1]*(1.-tau_tot)+(r[i-1]+r[i])/2.;
//			r_eff_interp = sqrt(4.*M_PI*r_eff_interp*r_eff_interp*F_mean/(4.*M_PI*(P_A)*(P_C)/4.*pow(T_eff, 4.)));
		}
		else{
//			r_eff = r_ini;
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
		printf("%f %e %e %e\n", t/86400., 4.*M_PI*r[n-1]*r[n-1]*(P_C)*E[2*(n-1)], r_eff, T_color);
//		printf("%f %e\n", t/86400., 4.*M_PI*r[n-1]*r[n-1]*(P_C)*E[2*(n-1)]);
	}
	fclose(fp);
	fclose(fl);
}

void init_E_U(double r_ini, double r_out, double r[], double rho[], double v_w[], double E[], double U[], const int nsize)
{
	int i;
	double dr = (r_out-r_ini)/((double)(nsize)+1./2.)/2.;
	double T;

	for(i = 0; i < nsize; i++){
		r[i] = r_ini-dr+2.*(double)(i+1)*dr;
		E[2*i] = (P_A)*pow(1.e+03, 4.)*pow(r[0]/r[i], 2.);
		T = pow(E[2*i]/(P_A), 0.25);
		E[2*i+1] = E[2*i];
		rho[i] = rho_csm(r[i]);
		v_w[i] = v_wind(r[i]);
//		U[2*i] = 1.5*rho[0]*(P_K)*T/(0.62*(MH))*pow(rho[i]/rho[0], 5./3.);
		U[2*i] = 1.5*rho[i]*(P_K)*T/(1.3*(MH));
		U[2*i+1] = U[2*i];
	}
	r[nsize] = r[nsize-1]+2.*dr;
}

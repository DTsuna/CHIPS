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

extern char csm[256];

void rad_transfer_csm(double, const char*, const char*, const char*);
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

void rad_transfer_csm(double r_out, const char *file_csm, const char *file_inp, const char *file_outp)
{
	FILE *fp, *fl;
	double F_max = 0., F_out = 0.;
	double E[2*NSIZE], U[2*NSIZE], r[NSIZE+1], E_old[NSIZE], rho[NSIZE], v_w[NSIZE], E0[NSIZE], U0[NSIZE];
	double r_ini, F_ini, u_ini, E_ini;
	double t, dt = 4.;
	double err = 0., tol = 1.e-06;
	double rho_ed[2];
	double tf[7000], rf[7000], Ff[7000], Ef[7000], uf[7000];
	int i = 0, j = 0, k, n = NSIZE, fsize, flag = 0;
	int F_neg_flag = 0;
	double dummy[8];
	double dr;
	double CFL = 1.00000000000000;
	char filename[256];

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

	t = tf[0];
	r_ini = rf[0];
	F_ini = Ff[0];
	u_ini = uf[0];
	E_ini = 0.;


	init_E_U(r_ini, r_out, r, rho, v_w, E, U, NSIZE);
	dr = r[1]-r[0];

/*
If r_fs(t) > r[0]-dr/3, r[0] = r[1], r[1] = r[2], ... .
Here, dr = r[N]-r[N-1] must be fixed.
for i = 0, 1, ..., n-1, 
X[i] = X[i+1], where X is physical quantity, i.e. E, U, rho.
*/
	while(t < tf[fsize-1]){

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
//		E_ini = Ef[j]*exp(log(Ef[j+1]/Ef[j])/log(tf[j+1]/tf[j])*log(t/tf[j]));
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
		F_ini += E_ini*u_ini;
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

	
		fprintf(fl, "%f %e\n", t/86400., 4.*M_PI*r[n-1]*r[n-1]*(P_C)*E[2*(n-1)]);
		printf("%f %e\n", t/86400., 4.*M_PI*r[n-1]*r[n-1]*(P_C)*E[2*(n-1)]);
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

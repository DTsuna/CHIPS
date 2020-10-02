#include <stdio.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include "constant.h"
#include "itgadvu.h"
#include "itgadve.h"
#include "srctrm.h"
#include "rhocsm.h"

#define N 3000

extern char csm[256];

void init_E_U(double, double, double[], double[], double[], double[], const int);

int main(void)
{
	FILE *fp, *fw, *fl;
	double E[2*N], U[2*N], r[N+1], E_old[N], rho[N];
	double T_g[N], T_r[N], mu[N];
	double r_ini, F_ini;
	double t, dt = 4.;
	double err = 0., tol = 1.e-06;
	double rho_ed[2], v_w = 1.e+07;
	double tf[2000], rf[2000], Ff[2000];
	int i = 0, j = 0, n = N, fsize, count = 0;
	double dammy[7];
	double dr;
	double time1, cpu_time;


	time1 = omp_get_wtime();
	sprintf(csm, "./inp-data/%s.txt", "CSM_1.5");

	fp = fopen("./inp-data/CSM_1.5_profile.txt", "r");
	fw = fopen("distribution_1.5.txt", "w");
	fl = fopen("luminosity_1.5.txt", "w");
	while(fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &dammy[0], &dammy[1], &dammy[2], &dammy[3], &dammy[4], &dammy[5], &dammy[6]) != EOF){
		tf[i] = dammy[0]*86400.000;
		rf[i] = dammy[4];
		Ff[i] = dammy[5];
		i++;
	}
	fsize = i;

	t = tf[0];
	r_ini = rf[0];
	F_ini = Ff[0];


	init_E_U(r_ini, 9.e+15, r, rho, E, U, N);
	dr = r[1]-r[0];

/*
If r_fs(t) > r[0]-dr/3, r[0] = r[1], r[1] = r[2], ... .
Here, dr = r[N]-r[N-1] must be fixed.
for i = 0, 1, ..., n-1, 
X[i] = X[i+1], where X is physical quantity, i.e. E, U, rho.
*/
	cpu_time = omp_get_wtime()-time1;
	cpu_time /= 60.;
while(cpu_time < 300. && count < 10000 && t < tf[fsize-1]){
	dt = t*0.001;
	for(i = 0; i < fsize; i++){
		if(t >= tf[i] && t < tf[i+1]){
			j = i;
		}
	}
	r_ini = rf[j]*exp(log(rf[j+1]/rf[j])/log(tf[j+1]/tf[j])*log(t/tf[j]));
	F_ini = Ff[j]*exp(log(Ff[j+1]/Ff[j])/log(tf[j+1]/tf[j])*log(t/tf[j]));
	if(r_ini > r[0]-dr/4.){
		n--;
		for(i = 0; i < n; i++){
			E[2*i] = E[2*(i+1)];
			U[2*i] = U[2*(i+1)];
			E[2*i+1] = E[2*(i+1)+1];
			U[2*i+1] = U[2*(i+1)+1];
			r[i] = r[i+1];
			rho[i] = rho[i+1];
		}
		r[n] = r[n+1];
	}
	printf("count = %d, %d %f %e %e %f\n", count, n, t/86400., r_ini, F_ini, cpu_time);

/*
At first, intergrate source term using implicit Euler method.
*/
	for(i = 0; i < n; i++){
		itg_src(E+2*i, U+2*i, rho[i], dt, tol);
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
	for(i = 0; i < n; i++){
		E_old[i] = E[2*i];
	}

	do{
		err = 0.;
		itg_adv_E(r_ini, F_ini, r, E, U, rho, dt, n);
		for(i = 0; i < n; i++){
			err += pow(1.-E[2*i+1]/E_old[i], 2.);
			E_old[i] = E[2*i+1];
		}

		err = sqrt(err/(double)n);
//		printf("err = %e\n", err);
	}while(err > tol);

	for(i = 0; i < n; i++){
		saha(rho[i], U[2*i+1], mu+i, T_g+i);
		T_r[i] = pow(E[2*i+1]/(P_A), 0.25);
	}

	t += dt;
	count++;
	cpu_time = omp_get_wtime()-time1;
	cpu_time /= 60.;
	for(i = 0; i < n; i++){
		if(count%100 == 0){
			fprintf(fw, "%e %e %e %e %e\n", r[i], E[2*i+1], U[2*i+1], T_r[i], T_g[i]);
		}
		E[2*i] = E[2*i+1];
		U[2*i] = U[2*i+1];
	}
	if(count%100 == 0){
		fprintf(fw, "\n\n");
	}

	fprintf(fl, "%f %e\n", t/86400., 4.*M_PI*r[n-1]*r[n-1]*(P_C)*E[2*(n-1)]);
}

	return 0;
}

void init_E_U(double r_ini, double r_out, double r[], double rho[], double E[], double U[], const int nsize)
{
	int i;
	double dr = (r_out-r_ini)/((double)(nsize)+1./2.)/2.;
	double T;

	printf("r_ini = %e cm\n", r_ini);
	for(i = 0; i < nsize; i++){
		r[i] = r_ini-dr+2.*(double)(i+1)*dr;
		E[2*i] = (P_A)*pow(1.e+03, 4.)*pow(r[0]/r[i], 2.);
		T = pow(E[2*i]/(P_A), 0.25);
		E[2*i+1] = E[2*i];
		rho[i] = rho_csm(r[i]);
		U[2*i] = 1.5*rho[0]*(P_K)*T/(0.62*(MH))*pow(rho[i]/rho[0], 5./3.);
		U[2*i+1] = U[2*i];
		printf("r[%d] = %e cm, E = %e, U = %e, rho = %e\n", i, r[i], E[2*i], U[2*i], rho[i]);
	}
	r[nsize] = r[nsize-1]+2.*dr;
	printf("r_out = %e cm\n", r[nsize]);
}

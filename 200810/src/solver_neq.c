#include <stdio.h>
#include <omp.h>
#include "solver_neq.h"
#include "boundary.h"
#include "diff_eq.h"
#include "irk.h"
#include "opacity.h"	

extern double t_exp;
extern double xfd, yrev[N], yfor[N];
extern pars pdt;

void solver_rev(double x_ini, double int_phys[], double egn[], FILE *fprev)
{
	double kappa, kappa_planck, R;
	int count = 0, k = 0, i, mem_count;
	double y[N], yout[N], yl[N], ylout[N], dydx[N];
	double T_g = 1.e+08, T_r = 1.e+04;
	double Mass = 0., dMass = 0., Mass_rev;
	double tmp_m = 0.;
	double err = 0.1, tol = 1.e-05;
	double dx, dx_min, dx_lim, mem_dx;
	double x;
	double source;
	double rho[2];
	double dydx0;

	Mass_rev = func_M_ej(x_ini, t_exp);
//	printf("Mass_rev = %e M_sun\n", Mass_rev/(M_SUN));
	boundary_rev(&x_ini, y, egn);
	memcpy(yout, y, sizeof(double)*N);

	x = 1.;
	dx_min = 1.e-8;
	dx_lim = 1.e-15;
	dx = dx_min;

	diff_eq_rev(x, y, dydx);
	dydx0 = fabs(dydx[1]);
	mem_dx = 1.;

//保存量を見て、ずれすぎてたら刻み幅を制御する
	do{
		if(count == 0 && fabs(T_g-T_r) > 1.){
			diff_eq_rev(x, y, dydx);
//		printf("dlnrho/dx = %e\n", dydx[1]);
			dx = dx_min*fabs(dydx0/dydx[1]);
		}
		if((fabs(T_g-T_r) < 500. && dx < 1.e-07)){
			break;
		}
		irk24(y, N, dx, x, yout, tol, &err, diff_eq_rev);
		dimless2dim(yout, yl, yrev);
		T_g = T_gas(yl);
		T_r = T_rad(yl);
		if(isnan(yout[0]) || yout[2] < 0.0 || T_g < T_r){
			count++;
		//	dx = dx_min*pow(1.5, -((double)count));
			dx *= 0.9;
			continue;
		}
		else{
			if(fabs(dx) > mem_dx){
				break;
			}
			mem_count = count;
			mem_dx = fabs(dx);
			count = 0;
		}

		rho[0] = exp(y[1]*yrev[1]);
		rho[1] = exp(yout[1]*yrev[1]);
		dMass = 4.*M_PI*(x+dx/2.0)*(x+dx/2.)*(rho[0]+rho[1])*0.5*dx*xfd*xfd*xfd;

		if(Mass+dMass > Mass_rev){
			break;
		}

		memcpy(y, yout, sizeof(double)*N);
		x += dx;
		Mass += dMass;
		for(i = 0; i < N; i++){
			yl[i] = y[i]*yrev[i];
		}
		T_r = T_rad(yl);
		T_g = T_gas(yl);
//			printf("%d %.15e %.10e %.10e %.10e %.10e %.10e %.10e %.8f\n", 
//			mem_count, x*xfd, exp(yl[0]), exp(yl[1]), yl[2], T_g, T_r, dx, Mass/Mass_rev);
		if(fprev != NULL){
			fprintf(fprev, "%.15e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.8f\n", 
			x*xfd, exp(yl[0]), exp(yl[1]), yl[2], yl[3], yl[4], func_Erad(yl), T_g, T_r, dx, Mass/Mass_rev);
		}
	}while(Mass < Mass_rev);
	if(Mass > Mass_rev){
		dx = set_final_dx(x, y, yout, Mass, Mass_rev, dMass, yrev);
		sub_Mass_shocked(x, dx, y, yout, yrev, int_phys, egn[0], diff_eq_rev);
	}
	if(fprev != NULL){
		fprintf(fprev, "\n\n");
		printf("Mass/Mass_rev = %f\n", Mass/Mass_rev);
	}
	else{
		int_phys[0] = exp(yl[0])+egn[0];
		int_phys[1] = yl[2]+yl[3];
		int_phys[2] = yl[4];
		int_phys[3] = x*xfd;
	}
}

void solver_for(double int_phys[], double ext_phys[], double egn[], FILE *fpfor)
{
	double kappa, kappa_planck, R;
	int count = 0, k = 0, i, l = 0;
	double y[N], yout[N], yl[N], ylout[N], dydx[N];
	double T_g = 1.e+08, T_r = 1.e+04;
	double Mass = 0., dMass = 0., Mass_for;
	double tmp_m = 0.;
	double err = 0.1, tol = 1.e-05;
	double x_ini;
	double dx, dx_min, dx_lim, mem_dx = 1.;
	double x;
	double source;
	double rho[2];
	double dydx0;

	//Mass_for = 4.*M_PI*pdt.D/(3.-pdt.s)*(pow(egn[3], 3.-pdt.s)-pow(pdt.R_p, 3.-pdt.s));
	Mass_for = 4.*M_PI*pdt.D/(3.-pdt.s)*(pow(egn[3], 3.-pdt.s));
//	printf("Mass_for = %e M_sun\n", Mass_for/(M_SUN));
	boundary_for(&x_ini, y, egn);
	memcpy(yout, y, sizeof(double)*N);

	x = egn[3]/xfd;
	dx_min = -1.e-03;
	dx_lim = 1.e-15;
	dx = dx_min;

	diff_eq_for(x, y, dydx);
	dydx0 = fabs(dydx[1]);

	do{
		l++;
		if(l == 10000){
			printf("for\n");
			l = 0;
		}
		if(count == 0 && fabs(T_g-T_r) > 50.){
			diff_eq_for(x, y, dydx);
//		printf("dlnrho/dx = %e\n", dydx[1]);
			dx = -fabs(dx_min*dydx0/dydx[1]);
		}
		if(fabs(T_g-T_r) < 100. && fabs(dx) < 1.e-08){
			break;
			dx *= 1.1;
		}
		irk24(y, N, dx, x, yout, tol, &err, diff_eq_for);
		dimless2dim(yout, yl, yfor);
		T_r = T_rad(yl);
		T_g = T_gas(yl);
		if(isnan(yout[0]) || T_r > T_g){
			count = 1;
		//	dx = dx_min*pow(1.5, -((double)count));
			dx *= 0.1;
			continue;
		}
		else{
			if(fabs(dx) > mem_dx){
				break;
			}
			mem_dx = fabs(dx);
			count = 0;
		}

		rho[0] = exp(y[1]*yfor[1]);
		rho[1] = exp(yout[1]*yfor[1]);
		dMass = -4.*M_PI*(x+dx/2.0)*(x+dx/2.)*(rho[0]+rho[1])*0.5*dx*xfd*xfd*xfd;

		if(Mass+dMass > Mass_for){
			break;
		}

		memcpy(y, yout, sizeof(double)*N);
		x += dx;
		Mass += dMass;
		for(i = 0; i < N; i++){
			yl[i] = y[i]*yfor[i];
		}
		T_r = T_rad(yl);
		T_g = T_gas(yl);
//		printf("%.15e %.10e %.10e %.10e %.10e %.10e %.10e %.8f mass conserv. %.15e momentum conserv. %.15e\n", 
//		x*xfd, exp(yl[0]), exp(yl[1]), yl[2], T_g, T_r, dx, Mass/Mass_for, x*x*xfd*xfd*exp(yl[0])*exp(yl[1]), exp(yl[1]+2.*yl[0])+yl[2]+yl[3]);
		if(fpfor != NULL){
			kappa_planck = kappa_p(exp(yl[1]), T_g);
			fprintf(fpfor, "%.15e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %e %.8f\n", 
			x*xfd, exp(yl[0]), exp(yl[1]), yl[2], yl[3], yl[4], func_Erad(yl), T_g, T_r, dx, sqrt(5./3.*yl[2]/exp(yl[1])), Mass/Mass_for);
		}
	}while(Mass < Mass_for);

	if(fpfor != NULL){
		printf("Mass/Mass_for = %f\n", Mass/Mass_for);
		fprintf(fpfor, "\n\n");
	}
	if(Mass+dMass > Mass_for){
		dx = -set_final_dx(x, y, yout, Mass, Mass_for, dMass, yfor);
		sub_Mass_shocked(x, dx, y, yout, yfor, ext_phys, -egn[1], diff_eq_for);
		ext_phys[0] *= -1.0;
	}
	else{
		ext_phys[0] = -exp(yl[0])+egn[1];
		ext_phys[1] = yl[2]+yl[3];
		ext_phys[2] = yl[4];
		ext_phys[3] = x*xfd;
	}
}

void dimless2dim(double y[], double yl[], double yfd[])
{
	int i;
	for(i = 0; i < N; i++){
		yl[i] = y[i]*yfd[i];
	}
}

void set_phys_cd(double x, double y[], double yfd[], double u_s, double phys[])
{
	double yl[N];
	double yy, p_r;
	dimless2dim(y, yl, yfd);
	phys[0] = exp(yl[0])+u_s;
	phys[1] = yl[2]+yl[3];
	phys[2] = yl[4];
	phys[3] = x*xfd;
}

void sub_Mass_shocked(double x, double dx, double y[], double yout[],  double yfd[], double phys[], double u_s, void (*derivs)(double, double[], double[]))
{
	int i;
	double err = 0.1, tol = 1.e-05;
	irk24(y, N, dx, x, yout, tol, &err, derivs);
	set_phys_cd(x, yout, yfd, u_s, phys);
}

double set_final_dx(double x, double y[], double yout[], double Mass, double Mass_shocked, double dMass, double yfd[])
{
	double dx, dx2;
	yout[1] = yout[1]*(Mass_shocked-Mass)/dMass+y[1]*(Mass+dMass-Mass_shocked)/dMass;
	dx2 = (Mass-Mass_shocked)/((exp(yfd[1]*y[1])+exp(yfd[1]*yout[1]))*M_PI*x)/(xfd*xfd*xfd);
	dx = dx2/(x+sqrt(x*x-dx2));
	return fabs(dx);
}


void solver(double x_ini, double egn[], double phys[], FILE *fprev, FILE *fpfor)
{
	int i = 0;
	double int_phys[4], ext_phys[4];
		solver_rev(x_ini, int_phys, egn, fprev);
		solver_for(int_phys, ext_phys, egn, fpfor);
//	printf("int_v = %e, ext_v = %e\n", int_phys[0], ext_phys[0]);
	for(i = 0; i < 4; i++){
		phys[i] = ext_phys[i]-int_phys[i];
	}
}

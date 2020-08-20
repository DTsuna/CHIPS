#include "solver.h"
#include "diff_eq.h"
#include "rk.h"
#include "main.h"
#include "calc_array.h"

extern double t_exp;
extern double xfd, yrev[N], yfor[N];
extern pars pdt;

//x_ini, int_phys, egnは全て次元を持つ量だが、中で行われている計算は全て無次元で行われている

void solver_rev(double x_ini, double int_phys[], double egn[])
{
	double y[2*N+1], yout[N], yl[N], ylout[N], dydx[N];
	double Mass = 0., dMass = 0., Mass_rev;
	double tmp_m = 0.;
	double err = 0.1;
	double dx, x;
	int i;
	double yy, R, p_r;
	double tol = 1.e-05;
	double chi, p_rad;
	double j, T_g, T_r, kappa;
	double dTdr;

	Mass_rev = func_M_ej(x_ini, t_exp);
	boundary_rev(&x_ini, y, egn);
	for(i = 0; i < N; i++){
		yl[i] = y[i]*yrev[i];
	}
	printf("%e %e %e %e %e\n", yl[0], yl[1], yl[2], yl[3], yl[4]);

	for(i = 0; i < N; i++){
		yout[i] = y[i];
	}
	x = x_ini;
	dx = 1.e-07;

	do{
		if(Mass+dMass <= Mass_rev){
			Mass += dMass;
		}
		else{
			break;
		}
//		dx = 1.e-9;
		irk24(y, N, dx, x, yout, tol, &err, diff_eq_rev);
		for(i = 0; i < N; i++){
			yl[i] = y[i]*yrev[i];
			ylout[i] = yout[i]*yrev[i];
		}
		T_g = T_gas(ylout);
		T_r = pow(ylout[3]/(P_A), 0.25);
		if(T_gas(ylout) < pow(ylout[3]/(P_A), 0.25)-1. || isnan(yout[0])){
			printf("aaa %.15e %e %e %e %e %e T_g = %e T_r = %e\n", x*xfd, ylout[0], ylout[1], ylout[2], ylout[3], ylout[4], T_g, T_r);
			diff_eq_rev(x, y, dydx);
			for(i = 0; i < N; i++){
				dydx[i] *= yrev[i]/xfd;
			}
			T_g = T_gas(yl);
			T_r = pow(yl[3]/(P_A), 0.25);
			dTdr = (P_K)/(MU)*(dydx[2]/yl[1]-yl[2]/yl[1]/yl[1]*dydx[1]);
			printf("dTdr = %e, dx = %e\n", dTdr, (T_r-T_g)/dTdr);
			dx = (T_r-T_g)/dTdr;
			continue;
		}
		memcpy(y, yout, sizeof(double)*N);
		dMass = 4.*M_PI*(x+dx/2.)*(x+dx/2.)*(y[1]+yout[1])/2.*dx*pow(xfd, 3.)*yrev[1];
		x += dx;
		for(i = 0; i < N; i++){
			yl[i] = y[i]*yrev[i];
		}
		T_g = T_gas(yl);
		T_r = pow(yl[3]/(P_A), 0.25);
		if(fabs(T_g-T_r) < 10.){
			dx = 1.e-06;
		}
		printf("%.15e %e %e %e %e %e T_g = %e T_r = %e %e\n", x*xfd, yl[0], yl[1], yl[2], yl[3], yl[4], T_g, T_r, Mass/Mass_rev);
	}while(Mass <= Mass_rev);
	x -= dx;
	yout[1] = yout[1]*Mass_rev/(Mass+dMass);
	tmp_m = (Mass_rev-Mass)*(1.-(yout[1]-y[1])/(2.*y[1]))/(4.*M_PI*x*y[1])/(pow(xfd, 3.)*yrev[1]);
	dx = tmp_m/(sqrt(tmp_m+x*x*0.25)+x*0.5);
	irk24(y, N, dx, x, yout, tol, &err, diff_eq_rev);
	for(i = 0; i < N; i++){
		yout[i] = yout[i]*yrev[i];
	}
	yy = yout[4]/((P_C)*yout[3]);
	chi = (2.-3.*yy+sqrt(4.+12.*yy-15.*yy*yy))/12.+yy*yy;
	p_r = chi*yout[3];
	int_phys[0] = yout[0]+egn[0];
	int_phys[1] = yout[2]+p_r;
	int_phys[2] = yout[4];
	int_phys[3] = (x+dx)*xfd;
}

void solver_rev_outp(double x_ini, double int_phys[], double egn[], FILE *fp)
{
	double y[2*N+1], yout[N], yl[N], ylout[N], dydx[N];
	double Mass = 0., dMass = 0., Mass_rev;
	double tmp_m = 0.;
	double err = 0.1;
	double dx, x;
	int i;
	double yy, R, p_r;
	double tol = 1.e-05;
	double chi, p_rad;
	double j, T_g, T_r, kappa;
	double dTdr;

	Mass_rev = func_M_ej(x_ini, t_exp);
	boundary_rev(&x_ini, y, egn);
	for(i = 0; i < N; i++){
		yl[i] = y[i]*yrev[i];
	}
	printf("%e %e %e %e %e\n", yl[0], yl[1], yl[2], yl[3], yl[4]);

	for(i = 0; i < N; i++){
		yout[i] = y[i];
	}
	x = x_ini;
	dx = 1.e-07;

	do{
		if(Mass+dMass <= Mass_rev){
			Mass += dMass;
		}
		else{
			break;
		}
//		dx = 1.e-9;
		irk24(y, N, dx, x, yout, tol, &err, diff_eq_rev);
		for(i = 0; i < N; i++){
			yl[i] = y[i]*yrev[i];
			ylout[i] = yout[i]*yrev[i];
		}
		T_g = T_gas(ylout);
		T_r = pow(ylout[3]/(P_A), 0.25);
		if(T_gas(ylout) < pow(ylout[3]/(P_A), 0.25)-1. || isnan(yout[0])){
			printf("aaa %.15e %e %e %e %e %e T_g = %e T_r = %e\n", x*xfd, ylout[0], ylout[1], ylout[2], ylout[3], ylout[4], T_g, T_r);
			diff_eq_rev(x, y, dydx);
			for(i = 0; i < N; i++){
				dydx[i] *= yrev[i]/xfd;
			}
			T_g = T_gas(yl);
			T_r = pow(yl[3]/(P_A), 0.25);
			dTdr = (P_K)/(MU)*(dydx[2]/yl[1]-yl[2]/yl[1]/yl[1]*dydx[1]);
			printf("dTdr = %e, dx = %e\n", dTdr, (T_r-T_g)/dTdr/xfd);
			dx = (T_r-T_g)/dTdr;
			continue;
		}
		memcpy(y, yout, sizeof(double)*N);
		dMass = 4.*M_PI*(x+dx/2.)*(x+dx/2.)*(y[1]+yout[1])/2.*dx*pow(xfd, 3.)*yrev[1];
		x += dx;
		for(i = 0; i < N; i++){
			yl[i] = y[i]*yrev[i];
		}
		T_g = T_gas(yl);
		T_r = pow(yl[3]/(P_A), 0.25);
		if(fabs(T_g-T_r) < 10.){
			dx = 1.e-06;
		}

		fprintf(fp, "%.15e %e %e %e %e %e %e %e %e\n", x*xfd, yl[0], yl[1], yl[2], yl[3], yl[4], T_g, T_r, Mass/Mass_rev);
	}while(Mass <= Mass_rev);
	x -= dx;
	yout[1] = yout[1]*Mass_rev/(Mass+dMass);
	tmp_m = (Mass_rev-Mass)*(1.-(yout[1]-y[1])/(2.*y[1]))/(4.*M_PI*x*y[1])/(pow(xfd, 3.)*yrev[1]);
	dx = tmp_m/(sqrt(tmp_m+x*x*0.25)+x*0.5);
	irk24(y, N, dx, x, yout, tol, &err, diff_eq_rev);
	for(i = 0; i < N; i++){
		yout[i] = yout[i]*yrev[i];
	}
	T_g = T_gas(yout);
	T_r = pow(yout[3]/(P_A), 0.25);
	fprintf(fp, "%.15e %e %e %e %e %e %e %e %e\n", (x+dx)*xfd, yout[0], yout[1], yout[2], yout[3], yout[4], T_g, T_r, Mass/Mass_rev);
	yy = yout[4]/((P_C)*yout[3]);
	chi = (2.-3.*yy+sqrt(4.+12.*yy-15.*yy*yy))/12.+yy*yy;
	p_r = chi*yout[3];
	int_phys[0] = yout[0]+egn[0];
	int_phys[1] = yout[2]+p_r;
	int_phys[2] = yout[4];
	int_phys[3] = (x+dx)*xfd;
}

void solver_for(double int_phys[], double ext_phys[], double egn[])
{
	double y[2*N+1], yout[N], yl[N], ylout[N], dydx[N];
	double Mass = 0., dMass = 0., Mass_for;
	double tmp_m = 0.;
	double err = 0.1;
	double dx, x;
	int i;
	double yy, R, p_r;
	double tol = 1.e-05;
	double chi, p_rad;
	double j, T_g, T_r, kappa;
	double dTdr;

	Mass_for = 4.*M_PI*pdt.D/(3.-pdt.s)*(pow(egn[3], 3.-pdt.s)-pow(pdt.R_p, 3.-pdt.s));
	boundary_for(&x, y, egn);
	for(i = 0; i < N; i++){
		yl[i] = y[i]*yfor[i];
	}
	printf("%e %e %e %e %e\n", yl[0], yl[1], yl[2], yl[3], yl[4]);

	for(i = 0; i < N; i++){
		yout[i] = y[i];
	}

	printf("x_ini = %e\n", x);
	dx = 1.e-07;

	do{
		if(Mass+dMass <= Mass_for){
			Mass += dMass;
		}
		else{
			break;
		}
//		dx = 1.e-9;
		irk24(y, N, -dx, x, yout, tol, &err, diff_eq_for);
		for(i = 0; i < N; i++){
			yl[i] = y[i]*yfor[i];
			ylout[i] = yout[i]*yfor[i];
		}
		T_g = T_gas(ylout);
		T_r = pow(ylout[3]/(P_A), 0.25);
		if(T_gas(ylout) < pow(ylout[3]/(P_A), 0.25)-1. || isnan(yout[0])){
			printf("aaa %.15e %e %e %e %e %e T_g = %e T_r = %e\n", x*xfd, ylout[0], ylout[1], ylout[2], ylout[3], ylout[4], T_g, T_r);
			diff_eq_rev(x, y, dydx);
			for(i = 0; i < N; i++){
				dydx[i] *= yfor[i]/xfd;
			}
			T_g = T_gas(yl);
			T_r = pow(yl[3]/(P_A), 0.25);
			dx = fabs(((MU)/(P_K)*yl[2]/T_r-yl[1])/dydx[1])/xfd;
			printf("dx = %e\n", dx);
		}
			continue;
		memcpy(y, yout, sizeof(double)*N);
		dMass = 4.*M_PI*(x+dx/2.)*(x+dx/2.)*(y[1]+yout[1])/2.*dx*pow(xfd, 3.)*yfor[1];
		x -= dx;
		for(i = 0; i < N; i++){
			yl[i] = y[i]*yfor[i];
		}
		T_g = T_gas(yl);
		T_r = pow(yl[3]/(P_A), 0.25);
		if(fabs(T_g-T_r) < 10.){
			dx = 1.e-06;
		}
		printf("%.15e %e %e %e %e %e T_g = %e T_r = %e %e\n", x*xfd, yl[0], yl[1], yl[2], yl[3], yl[4], T_g, T_r, Mass/Mass_for);
	}while(Mass <= Mass_for);
	x += dx;
	yout[1] = yout[1]*Mass_for/(Mass+dMass);
	tmp_m = (Mass_for-Mass)*(1.-(yout[1]-y[1])/(2.*y[1]))/(4.*M_PI*x*y[1])/(xfd*xfd*xfd*yfor[1]);
	dx = tmp_m/(sqrt(-tmp_m+x*x*0.25)+x*0.5);
	irk24(y, N, -dx, x, yout, tol, &err, diff_eq_for);

	for(i = 0; i < N; i++){
	        yout[i] = yout[i]*fabs(yfor[i]);
	}
	yy = yout[4]/((P_C)*yout[3]);
	chi = (2.-3.*yy+sqrt(4.+12.*yy-15.*yy*yy))/12.+yy*yy;
	p_r = chi*yout[3];
	ext_phys[0] = yout[0]+egn[1];
	ext_phys[1] = yout[2]+p_r;
	ext_phys[2] = yout[4];
	ext_phys[3] = (x-dx)*xfd;
}

void solver_for_outp(double int_phys[], double ext_phys[], double egn[], FILE *fp)
{
	double y[2*N+1], yout[N], yl[N], ylout[N], dydx[N];
	double Mass = 0., dMass = 0., Mass_for;
	double tmp_m = 0.;
	double err = 0.1;
	double dx, x;
	int i, count = 0, flag = 0;
	double yy, R, p_r;
	double tol = 1.e-05;
	double chi, p_rad;
	double j, T_g, T_r, kappa;
	double dTdr;
	double c_s;

	Mass_for = 4.*M_PI*pdt.D/(3.-pdt.s)*(pow(egn[3], 3.-pdt.s)-pow(pdt.R_p, 3.-pdt.s));
	boundary_for(&x, y, egn);
	for(i = 0; i < N; i++){
		yl[i] = y[i]*yfor[i];
	}

	for(i = 0; i < N; i++){
		yout[i] = y[i];
	}

	dx = 0.2e-06;

	do{
		if(Mass+dMass <= Mass_for){
			Mass += dMass;
		}
		else{
			break;
		}
//		dx = 0.2e-06;
//		dx = 1.e-9;
//		rk4(y, N, -dx, x, yout, diff_eq_for);
		irk24(y, N, -dx, x, yout, tol, &err, diff_eq_for);
		diff_eq_for(x, y, dydx);
		for(i = 0; i < N; i++){
			yl[i] = y[i]*yfor[i];
			ylout[i] = yout[i]*yfor[i];
			dydx[i] *= yfor[i]/xfd;
		}
		T_g = T_gas(ylout);
		T_r = pow(ylout[3]/(P_A), 0.25);
		if(isnan(yout[0])){
			printf("aaa %.15e %e %e %e %e %e T_g = %e T_r = %e\n", x*xfd, ylout[0], ylout[1], ylout[2], ylout[3], ylout[4], T_g, T_r);
			T_g = T_gas(yl);
			T_r = pow(yl[3]/(P_A), 0.25);
//			dx = 0.5*fabs(((MU)/(P_K)*yl[2]/T_r-yl[1])/dydx[1])/xfd;
			if(flag == 0){
				dx = 0.9*fabs(((MU)/(P_K)*yl[2]/T_r-yl[1])/dydx[1])/xfd;
				dx = fmin(dx, 1.e-10);
				printf("dx = %e, y[0] = %e y[1] = %e y[2] = %e y[3] = %e T_g = %e T_r = %e\n", dx, yl[0], yl[1], yl[2], yl[3], T_g, T_r);
				printf("dx = %e, %e %e %e %e %e\n", dx, dydx[0], dydx[1], dydx[2], dydx[3], dydx[4]);
				flag = 1;
			}
			else{
				dx *= 0.9;
			}
			if(isinf(dx)){
				c_s = sqrt(5./3.*yl[2]/yl[1]);
				dx = 0.9*fabs(((MU)/(P_K)*yl[2]/T_r-yl[1])/dydx[1])/xfd;
				printf("dx = %e, y[0] = %e y[1] = %e y[2] = %e y[3] = %e\n", dx, yl[0], yl[1], yl[2], yl[3]);
				printf("dydx[0] = %e dydx[1] = %e dydx[2] = %e dydx[3] = %e\n", dydx[0], dydx[1], dydx[2], dydx[3]);
				printf("%.15e %e %e %e %e %e %e c_s = %e\n", x*xfd, yl[0], yl[1], yl[2], T_g, T_r, Mass/Mass_for, c_s);
				
				break;
			}
			printf("dx = %e\n", dx);
			continue;
		}
		flag = 0;
		memcpy(y, yout, sizeof(double)*N);
		dMass = 4.*M_PI*(x+dx/2.)*(x+dx/2.)*(y[1]+yout[1])/2.*dx*pow(xfd, 3.)*yfor[1];
		x -= dx;
		for(i = 0; i < N; i++){
			yl[i] = y[i]*yfor[i];
		}
		T_g = T_gas(yl);
		T_r = pow(yl[3]/(P_A), 0.25);
		if(fabs(T_g-T_r) < 10.){
			dx = 1.e-06;
		}
		//dx = 1.e-06;
		count++;
		fprintf(fp, "%.15e %e %e %e %e %e %e %e %e %e\n", x*xfd, yl[0], yl[1], yl[2], yl[3], yl[4], T_g, T_r, Mass/Mass_for, dx);
		if(count%100 == 0){
			c_s = sqrt(5./3.*ylout[2]/ylout[1]);
			printf("%.15e %e %e %e %e %e %e c_s = %e\n", x*xfd, yl[0], yl[1], yl[2], T_g, T_r, Mass/Mass_for, c_s);
				printf("dydx[0] = %e dydx[1] = %e dydx[2] = %e dydx[3] = %e\n", dydx[0], dydx[1], dydx[2], dydx[3]);
		}
	}while(Mass <= Mass_for);
	x += dx;
	yout[1] = yout[1]*Mass_for/(Mass+dMass);
	tmp_m = (Mass_for-Mass)*(1.-(yout[1]-y[1])/(2.*y[1]))/(4.*M_PI*x*y[1])/(xfd*xfd*xfd*yfor[1]);
	dx = tmp_m/(sqrt(-tmp_m+x*x*0.25)+x*0.5);
	irk24(y, N, -dx, x, yout, tol, &err, diff_eq_for);

	for(i = 0; i < N; i++){
	        yout[i] = yout[i]*fabs(yfor[i]);
	}
	yy = yout[4]/((P_C)*yout[3]);
	chi = (2.-3.*yy+sqrt(4.+12.*yy-15.*yy*yy))/12.+yy*yy;
	p_r = chi*yout[3];
	ext_phys[0] = yout[0]+egn[1];
	ext_phys[1] = yout[2]+p_r;
	ext_phys[2] = yout[4];
	ext_phys[3] = (x-dx)*xfd;
}

void dimless2dim(double y[], double yl[], double yfd[])
{
	int i;
	for(i = 0; i < N; i++){
		yl[i] = y[i]*yfd[i];
	}
}

void dim2dimless(double y[], double yl[], double yfd[])
{
	int i;
	for(i = 0; i < N; i++){
		y[i] = yl[i]/yfd[i];
	}
}

void jump_calc(double y[], double yout[], double yfd[])
{
	double T_g, T_r;
	double yl[N];
	double rho;
	int i;

	dimless2dim(y, yl, yfd);

	T_g = T_gas(yl);
	T_r = pow(yl[3]/(P_A), 0.25);
	yout[2] = yl[2];
	yout[3] = yl[3];
	yout[4] = yl[4];
	rho = yl[1];
	yout[1] = rho;

	do{
		rho = yout[1];
		yout[1] = mmw(rho, T_r)*(MH)/(P_K)*yout[2]/T_r;
	}while(fabs(1.-yout[1]/rho) > 1.e-10);
	
	yout[0] = yl[0]*yl[1]/yout[1];
	dim2dimless(y, yout, yfd);
	memcpy(yout, y, sizeof(double)*N);
}

void solver_f(double int_phys[], double out_phys[], double egn[], FILE *fp)
{
	double y[N], yout[N], yl[N], ylout[N], dydx[N];
	double Mass = 0., dMass = 0., Mass_for;
	double tmp_m = 0.;
	double err = 0.1;
	double dx, dx_min, x, x0;
	int i, count = 0, flag = 0;
	double tol = 1.e-05;
	double j, T_g, T_r, kappa, p_r, yy;
	double dTdr;
	double c_s;
	double C[N];
	double (*y_ana[])(double, double, double[]) = {ana_sol_kconst_v, ana_sol_kconst_rho, ana_sol_kconst_pgas, ana_sol_kconst_erad, ana_sol_kconst_flux};

	Mass_for = 4.*M_PI*pdt.D/(3.-pdt.s)*(pow(egn[3], 3.-pdt.s)-pow(pdt.R_p, 3.-pdt.s));
	boundary_for(&x, y, egn);
	memcpy(yout, y, sizeof(double)*N);
	dimless2dim(y, yl, yfor);

	printf("%e\n", Mass_for);
//dx must be an negative value so that the same subroutine can be shared.
	dx_min = -2.e-06;
	dx = dx_min;
	do{
		T_g = T_gas(yl);
		T_r = T_rad(yl);
		if(fabs(dx) < 2.5e-16){
			break;
		}
		if(Mass < Mass_for){
			Mass += dMass;
		}
		else{
			break;
		}
		if(fabs(T_g-T_r) < 1.){
			dx = -1.e-9;
		}
		irk24(y, N, dx, x, yout, tol, &err, diff_eq_for);
		diff_eq_for(x, y, dydx);
		
		if(isnan(yout[0]) || isnan(yout[1])){
			dx = -dx_estim(x, y, yfor, 0.9, diff_eq_for);
			dx = fmax(dx, dx_min*pow(2., -(double)count));
			count++;
			continue;
		}
		memcpy(y, yout, sizeof(double)*N);
		dimless2dim(y, yl, yfor);
		dMass = 4.*M_PI*(x+dx/2.)*(x+dx/2.)*(y[1]+yout[1])/2.*dx*pow(xfd, 3.)*yfor[1];
		dMass = fabs(dMass);
		x += dx;
		yy = fabs(yl[4]/((P_C)*yl[3]));
		p_r = chi_f(yy)*yl[3];
		fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %e %e %e\n", 
			x*xfd, yl[0], yl[1], yl[2], yl[3], yl[4], p_r, yl[2]+p_r, T_gas(yl), T_rad(yl), kappa_r(yl[1], T_gas(yl)), kappa_p(yl[1], T_gas(yl)), err);
		diff_eq_for(x, y, dydx);
		printf("%e %e %e %e %e Mass/Mass_for = %e\n", dydx[0], dydx[1], dydx[2], dydx[3], dydx[4], Mass/Mass_for);
		printf("%e %e %e %e %e Mass/Mass_for = %e\n", y[0], y[1], y[2], y[3], y[4], Mass/Mass_for);
		printf("%e %e %e %e %e Mass/Mass_for = %e\n", yfor[0], yfor[1], yfor[2], yfor[3], yfor[4], Mass/Mass_for);
	}while(Mass <= Mass_for);
	dimless2dim(y, yl, yfor);
	c_s = sqrt(5./3.*yl[2]/yl[1]);
	printf("v*v-c_s*c_s = %e\n", yl[0]*yl[0]-c_s*c_s);
	jump_calc(y, yout, yfor);
	dimless2dim(yout, yl, yfor);
	c_s = sqrt(5./3.*yl[2]/yl[1]);
	printf("v*v-c_s*c_s = %e\n", yl[0]*yl[0]-c_s*c_s);
		fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %e %e\n", 
			x*xfd, yl[0], yl[1], yl[2], yl[3], yl[4], T_gas(yl), T_rad(yl), kappa_r(yl[1], T_gas(yl)), kappa_p(yl[1], T_gas(yl)));
//	printf("%e %e %e %e %e %e %e\n", yl[0], yl[1], yl[2], yl[3], yl[4], T_gas(yl), T_rad(yl));
	dimless2dim(yout, yl, yfor);
	memcpy(yfor, yl, sizeof(double)*N);
	for(i = 0; i < N; i++){
		y[i] = 1.;
	}
	memcpy(yout, y, sizeof(double)*N);
	dx *= 0.1;
	dx_min = dx;
	count = 0;
	x0 = x;

//	for(i = 0; i < 1000; i++){	
//	fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e\n", x*xfd, 
//	y_ana[0](x0*xfd, (x)*xfd, yfor), y_ana[1](x0*xfd, (x)*xfd, yfor),
//	y_ana[2](x0*xfd, (x)*xfd, yfor), y_ana[3](x0*xfd, (x)*xfd, yfor), y_ana[4](x0*xfd, (x)*xfd, yfor));
//	x += dx_min*100000.;
//	}
//	fprintf(fp, "\n\n");
//
//	x = x0;
//	
//	y_ana[0] = ana_sol_kprho_v;
//	y_ana[1] = ana_sol_kprho_rho;
//	y_ana[2] = ana_sol_kprho_pgas;
//	y_ana[3] = ana_sol_kprho_erad;
//	y_ana[4] = ana_sol_kprho_flux;
//	
//	for(i = 0; i < 1000; i++){	
//	fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e\n", x*xfd, 
//	y_ana[0](x0*xfd, (x)*xfd, yfor), y_ana[1](x0*xfd, (x)*xfd, yfor),
//	y_ana[2](x0*xfd, (x)*xfd, yfor), y_ana[3](x0*xfd, (x)*xfd, yfor), y_ana[4](x0*xfd, (x)*xfd, yfor));
//	x += dx_min*100000.;
//	}
//	fprintf(fp, "\n\n");
//
//	x = x0;

	for(count = 0; count < 100000; count++){
		irk24(y, N, dx, x, yout, tol, &err, diff_eq_for);
		if(isnan(yout[0])){
			dx *= 0.9;
			continue;
		}
		dMass = 4.*M_PI*(x+dx/2.)*(x+dx/2.)*(y[1]+yout[1])/2.*dx*pow(xfd, 3.)*yfor[1];
		dMass = fabs(dMass);
		Mass += dMass;
		dimless2dim(yout, yl, yfor);
		x += dx;
		if(count%100==0){
		fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %e %e\n", 
			x*xfd, yl[0], yl[1], yl[2], yl[3], yl[4], T_gas(yl), T_rad(yl), kappa_r(yl[1], T_gas(yl)), kappa_p(yl[1], T_gas(yl)));}
		diff_eq_for(x, yout, dydx);
		c_s = sqrt(5./3.*yl[2]/yl[1]);
//		printf("count = %d %e %e %e %e %e %e %e (v/c_s)**2 = %e\n", count, yl[0], yl[1], yl[2], yl[3], yl[4], Mass/Mass_for, err, yl[0]*yl[0]/c_s/c_s);
		dimless2dim(y, yl, yfor);
		memcpy(y, yout, sizeof(double)*N);
		dimless2dim(yout, yl, yfor);
		diff_eq_for(x, y, dydx);
		for(i = 0; i < N; i++){
			dydx[i] *= yfor[i]/xfd;
		}
		devide_array(yl, dydx, C, N);
		for(i = 0; i < N; i++){
			C[i] = fabs(C[i])/xfd;
		}

		dx = fabs(min_array(C, N));
		dx = -0.2*dx;
		if(count > 1000 && count < 2000){
			dx *= 20.;
		}
		else if(count >= 2000){
			dx *= 50.;
		}
//		printf("%e %e %e %e %e %e %e\n", x*xfd, dx*xfd, dydx[0]*xfd/yfor[0], dydx[1]*xfd/yfor[1], dydx[2]*xfd/yfor[2], dydx[3]*xfd/yfor[3], dydx[4]*xfd/yfor[4]);
		if(count%100 == 0){
		printf("count = %d %.14e %e %e %e %e %e %e Mass/Mass_for = %e\n", 
		count, x*xfd, dx*xfd, dydx[0], dydx[1], dydx[2], dydx[3], dydx[4], Mass/Mass_for);
		}
	}
	printf("Mass/Mass_for = %e\n", Mass/Mass_for);
}

double dx_estim(double x, double y[], double yfd[], double alpha, void (*derivs)(double, double[], double[]))
//NOTE that this function outputs a positive value, so be careful when you use this subroutine.
{
	double yl[N], dydx[N], T_r, dx;
	(*derivs)(x, y, dydx);
	dimless2dim(y, yl, yfd);
	T_r = T_rad(yl);
	dx = alpha*fabs(((MU)/(P_K)*yl[2]/T_r-yl[1])/(dydx[1]*yl[1]));
	return dx;
}

double ana_sol_kconst_rho(double r0, double r, double y0[]) //次元あり
{
	double kappa = kappa_r(y0[1], T_gas(y0)), c_s0 = sqrt(5./3.*y0[2]/y0[1]);
	double a = kappa*y0[4]*r0/(c_s0*c_s0*(P_C));
	return pow(1.+2./3.*a*(1.-r0/r), 2./3.)*y0[1];
}

double ana_sol_kconst_pgas(double r0, double r, double y0[])
{
	double rho = ana_sol_kconst_rho(r0, r, y0);
	return y0[2]/pow(y0[1]/rho, 5./3.);
}

double ana_sol_kconst_v(double r0, double r, double y0[])
{
	double rho = ana_sol_kconst_rho(r0, r, y0);
	return r0*r0*y0[0]*y0[1]/(r*r*rho);
}

double ana_sol_kconst_erad(double r0, double r, double y0[])
{
	double y[N] = {}, T_g;
	y[1] = ana_sol_kconst_rho(r0, r, y0);
	y[2] = ana_sol_kconst_pgas(r0, r, y0);
	T_g = T_gas(y);
	return (P_A)*pow(T_g, 4.);
}

double ana_sol_kconst_flux(double r0, double r, double y0[])
{
	return r0*r0*y0[4]/(r*r);
}

double ana_sol_kprho_rho(double r0, double r, double y0[]) //次元あり
{
	double kappa = kappa_r(y0[1], T_gas(y0)), c_s0 = sqrt(5./3.*y0[2]/y0[1]);
	double a = kappa*y0[4]*r0/(c_s0*c_s0*(P_C));
	return pow(1.-1./3.*a*(1.-r0/r), -3.)*y0[1];
}

double ana_sol_kprho_pgas(double r0, double r, double y0[])
{
	double rho = ana_sol_kprho_rho(r0, r, y0);
	return y0[2]/pow(y0[1]/rho, 5./3.);
}

double ana_sol_kprho_v(double r0, double r, double y0[])
{
	double rho = ana_sol_kprho_rho(r0, r, y0);
	return r0*r0*y0[0]*y0[1]/(r*r*rho);
}

double ana_sol_kprho_erad(double r0, double r, double y0[])
{
	double y[N] = {}, T_g;
	y[1] = ana_sol_kprho_rho(r0, r, y0);
	y[2] = ana_sol_kprho_pgas(r0, r, y0);
	T_g = T_gas(y);
	return (P_A)*pow(T_g, 4.);
}

double ana_sol_kprho_flux(double r0, double r, double y0[])
{
	return r0*r0*y0[4]/(r*r);
}

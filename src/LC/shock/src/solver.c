#include <stdio.h>
#include <math.h>
#include <string.h>
#include "solver.h"
#include "boundary.h"
#include "constant.h"
#include "rk.h"
#include "diff_eq.h"
#include "function.h"

extern double t_exp;

double func_dM(double x, double dx, double rho, double rho_out)
{
	double x_m = x+0.5*dx, rho_m = (rho+rho_out)*0.5;
	return 4.0*M_PI*x_m*x_m*rho_m*fabs(dx);
}

void solver_rev(double x_ini, double int_phys[], double egn[], int *info)
{
	double x = x_ini, y[N], yout[N];
	double dx = 0.0;
	double M = 0.0, dM = 0.0, M_rev = func_M_ej(x_ini, t_exp);
	int bflag = 0;

	boundary(x, yout, egn, bflag, info);
	
	do{
		x += dx;
		M += dM;
		dx = 1.0e-04*x_ini;
		memcpy(y, yout, sizeof(double)*N);
		rk4(y, N, dx, x, yout, diff_eq);
		dM = func_dM(x, dx, y[1], yout[1]);
	}while(M+dM < M_rev);

	dx = dx*(M_rev-M)/dM;
	rk4(y, N, dx, x, yout, diff_eq);

	set_phys(x, dx, yout, int_phys, egn, 0);
}

void solver_for(double int_phys[], double ext_phys[], double egn[], int flag, int *info)
{
	double x = egn[3], y[N], yout[N];
	double dx = 0.0;
	double M = 0.0, dM = 0.0;
	double M_for = func_M_csm(x, t_exp);
	int bflag = 1;

	boundary(x, yout, egn, bflag, info);

	if(flag == 0){
		do{
			x += dx;
			M += dM;
			dx = -1.0e-04*egn[3];
			memcpy(y, yout, sizeof(double)*N);
			rk4(y, N, dx, x, yout, diff_eq);
			dM = func_dM(x, dx, y[1], yout[1]);
		}while(M+dM < M_for);
		dx = dx*(M_for-M)/dM;
		rk4(y, N, dx, x, yout, diff_eq);
	}
	else if(flag == 1){
		do{
			x += dx;
			M += dM;
			dx = -1.0e-04*egn[3];
			memcpy(y, yout, sizeof(double)*N);
			rk4(y, N, dx, x, yout, diff_eq);
			dM = func_dM(x, dx, y[1], yout[1]);
		}while(x+dx > int_phys[3]);
		dx = -fabs(int_phys[3]-x);
		rk4(y, N, dx, x, yout, diff_eq);
	}
	set_phys(x, dx, yout, ext_phys, egn, 1);
}

void solver(double x_ini, double phys[], double egn[], int flag, int *info)
{
	double int_phys[4], ext_phys[4];
	int i;

	solver_rev(x_ini, int_phys, egn, info);
	solver_for(int_phys, ext_phys, egn, flag, info);
	for(i = 0; i < 4; i++){
		phys[i] = ext_phys[i]-int_phys[i];
	}
}

void set_phys(double x, double dx, double y[], double phys[], double egn[], int flag)
{
	phys[0] = y[0]+egn[flag];
	phys[1] = p_tot(y);
	phys[2] = y[3];
	phys[3] = x+dx;
}

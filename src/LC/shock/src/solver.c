#include <stdio.h>
#include <math.h>
#include <string.h>
#include "solver.h"
#include "boundary.h"
#include "constant.h"
#include "rk.h"
#include "diff_eq.h"
#include "function.h"
#include "pars.h"
#include "opacity.h"

extern pars pdt;
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


/*Profile output routine.*/
void solver_rev_outp(double x_ini, double int_phys[], double egn[], int *info, char *filename)
{
	double x = x_ini, y[N], yout[N];
	double dx = 0.0;
	double M = 0.0, dM = 0.0, M_rev = func_M_ej(x_ini, t_exp);
	int bflag = 0;
	FILE *fp;

	fp = fopen(filename, "w");

	boundary(x, yout, egn, bflag, info);
	fprintf(fp, "#radius, velocity, density, temperature, flux, rosseland mean\n");
	fprintf(fp, "# t = %f days\n", t_exp/86400.);
	fprintf(fp, "%e %e %e %e %e %e\n", x, yout[0]+egn[0], yout[1], yout[2], yout[3], kappa_r(yout[1], yout[2]));
	
	do{
		x += dx;
		M += dM;
		dx = 5.0e-04*x_ini;
		memcpy(y, yout, sizeof(double)*N);
		rk4(y, N, dx, x, yout, diff_eq);
		dM = func_dM(x, dx, y[1], yout[1]);
		if(M+dM < M_rev){
			fprintf(fp, "%e %e %e %e %e %e\n", x+dx, yout[0]+egn[0], yout[1], yout[2], yout[3], kappa_r(yout[1], yout[2]));
		}
	}while(M+dM < M_rev);

	dx = dx*(M_rev-M)/dM;
	rk4(y, N, dx, x, yout, diff_eq);
	fprintf(fp, "%e %e %e %e %e %e\n", x+dx, yout[0]+egn[0], yout[1], yout[2], yout[3], kappa_r(yout[1], yout[2]));
	fprintf(fp, "\n");
	fclose(fp);

	set_phys(x, dx, yout, int_phys, egn, 0);
}

//outp_array[0]=M, outp_array[1]=rho_ave, outp_array[2]=T_ave, outp_array[3]=v_ave, outp_array[4]=tau
void solver_for_outp(double int_phys[], double ext_phys[], double egn[], int flag, int *info, double outp_array[], char *filename)
{
	static int count = 0;
	FILE *fp;

	double xrec[2000], vrec[2000], rhorec[2000], Trec[2000], Frec[2000];
	double x = egn[3], y[N], yout[N];
	double dx = 0.0;
	double M = 0.0, dM = 0.0;
	double M_for = func_M_csm(x, t_exp);
	double rho_ave, T_ave, v_ave = 0., U_ave = 0., tau = 0., mom_ave = 0.;
	double dV = 4./3.*M_PI*(egn[3]*egn[3]*egn[3]-int_phys[3]*int_phys[3]*int_phys[3]);
	double kappa;
	int bflag = 1, i = 0, j;

	fp = fopen(filename, "a+");

	boundary(x, yout, egn, bflag, info);
	xrec[0] = x;
	vrec[0] = yout[0];
	rhorec[0] = yout[1];
	Trec[0] = yout[2];
	Frec[0] = yout[3];

	if(flag == 0){
		do{
			x += dx;
			M += dM;
			dx = -5.0e-04*egn[3];
			memcpy(y, yout, sizeof(double)*N);
			rk4(y, N, dx, x, yout, diff_eq);
			dM = func_dM(x, dx, y[1], yout[1]);
			i++;
			xrec[i] = x+dx;
			vrec[i] = yout[0];
			rhorec[i] = yout[1];
			Trec[i] = yout[2];
			Frec[i] = yout[3];
		}while(M+dM < M_for);
		dx = dx*(M_for-M)/dM;
		rk4(y, N, dx, x, yout, diff_eq);
		i++;
		xrec[i] = x+dx;
		vrec[i] = yout[0];
		rhorec[i] = yout[1];
		Trec[i] = yout[2];
		Frec[i] = yout[3];
	}
	else if(flag == 1){
		do{
			x += dx;
			M += dM;
			dx = -5.0e-04*egn[3];
			memcpy(y, yout, sizeof(double)*N);
			rk4(y, N, dx, x, yout, diff_eq);
			dM = func_dM(x, dx, y[1], yout[1]);
			i++;
			xrec[i] = x+dx;
			vrec[i] = yout[0];
			rhorec[i] = yout[1];
			Trec[i] = yout[2];
			Frec[i] = yout[3];
		}while(x+dx > int_phys[3]);
		dx = -fabs(int_phys[3]-x);
		rk4(y, N, dx, x, yout, diff_eq);
		xrec[i] = x+dx;
		vrec[i] = yout[0];
		rhorec[i] = yout[1];
		Trec[i] = yout[2];
		Frec[i] = yout[3];
	}
	set_phys(x, dx, yout, ext_phys, egn, 1);
	M = 0.;
	for(j = i; j >= 0; --j){
		kappa = kappa_r(rhorec[j], Trec[j]);
		if(j > 0){
			tau += kappa*rhorec[j]*fabs(xrec[j]-xrec[j-1]);
			M += 4.*M_PI*rhorec[j]*xrec[j]*xrec[j]*fabs(xrec[j]-xrec[j-1]);
			U_ave += (3./2.*rhorec[j]/((MH)*(MU))*(P_K)*Trec[j]+(P_A)*pow(Trec[j], 4.))*(4.*M_PI*xrec[j]*xrec[j]*fabs(xrec[j]-xrec[j-1]));
			mom_ave += (rhorec[j]*(vrec[j]+egn[1]))*4.*M_PI*xrec[j]*xrec[j]*fabs(xrec[j]-xrec[j-1]);
		}
		fprintf(fp, "%e %e %e %e %e %e\n", xrec[j], vrec[j]+egn[1], rhorec[j], Trec[j], Frec[j], kappa_r(rhorec[j], Trec[j]));
	}
	U_ave = U_ave/dV;
	v_ave = mom_ave/M;
	rho_ave = M/dV;
	T_ave = get_T_from_U(U_ave, rho_ave);


	outp_array[0] = M;
	outp_array[1] = rho_ave;
	outp_array[2] = T_ave;
	outp_array[3] = v_ave;
	outp_array[4] = tau;

	fclose(fp);
	count++;

}

void solver_outp(double x_ini, double phys[], double egn[], int flag, int *info, double outp_array[], char *filename)
{
	double int_phys[4], ext_phys[4];
	int i;

	solver_rev_outp(x_ini, int_phys, egn, info, filename);
	solver_for_outp(int_phys, ext_phys, egn, flag, info, outp_array, filename);
	for(i = 0; i < 4; i++){
		phys[i] = ext_phys[i]-int_phys[i];
	}
}

double get_T_from_U(double U, double rho)
{
	double dT, T = pow(U/(P_A), 0.25);
	double tol = 1.e-10;
	int count = 0;

	do{
		dT = -((P_A)*pow(T, 4.)+1.5*rho/((MU)*(MH))*(P_K)*T-U)/(4.*(P_A)*T*T*T+1.5*rho/((MU)*(MH))*(P_K));
		T += dT;
		count++;
		if(count > 100){
			printf("failed to derive T. abnormal end.\n");
			exit(EXIT_FAILURE);
		}
	}while(fabs(dT/T) > tol);
	
//	printf("count = %d, T = %e K.\n", count, T);
	return T;
}

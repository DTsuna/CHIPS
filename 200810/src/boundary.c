#include "boundary.h"
#include "main.h"

extern double t_exp;
extern pars pdt;
extern double xfd, yrev[N], yfor[N];

//インプットの*xは次元ありだが、最後に無次元量に変換される。

void boundary_rev(double *x, double y[], double egn[])
{
	double y_up[2];
	double gamma_g = 5./3.;
	double Mach = 300.;
	double jump = (gamma_g+1.)/(gamma_g-1.);
	double F0, e_r0, e_r1;
	double f0, f1, chi0, chi1;
	double tau;
	double xx, yy;
	int i;

	e_r0 = 3.*pdt.alpha*pdt.E_ej/(4.*M_PI*pdt.R_p*pdt.R_p*pdt.R_p)*pow(pdt.R_p/(*x), 4.);
	jump = (gamma_g+1.)*Mach*Mach/((gamma_g-1.)*Mach*Mach+2.);

	y_up[0] = (*x)/t_exp-egn[0];
	y_up[1] = rho_ej((*x), t_exp);
	y[0] = y_up[0]/jump;
	y[1] = y_up[1]*jump;
	y[2] = y_up[0]*y_up[1]*(y_up[0]-y[0]);
	tau = tau_ej((*x), t_exp);
	f0 = 1.0/tau;
	chi0 = func_chi_f(f0);
	yy = jump*jump*(1.0+chi0)*(1.0+chi0)/f0/f0;
	chi1 = (3.*yy+4.)/(5.*yy-4.+4.*sqrt(yy*(yy-4.)));
	f1 = sqrt((-3.0+10.*chi1-3.0*chi1*chi1)/4.0);
	F0 = (P_C)*e_r0/tau;
	e_r1 = F0/((P_C)*f1);
	y[3] = chi1*e_r1;

	xx = 3.*jump;
	y[4] = F0;
	
//	printf("%.10e %.10e %.10e %.10e %.10e %e\n", y[0], y[1], y[2], y[3], y[4], y[4]/((P_C)*y[3]));
	
	y[0] = log(y[0]);
	y[1] = log(y[1]);

	for(i = 0; i < N; i++){
		y[i] /= yrev[i];
//		y[i] = 1.;
	}
//	xfd = *x;
//	*x = 1.;
}

void boundary_for(double *x, double y[], double egn[])
{
	double y_up[2];
	double gamma_g = 5./3.;
	double Mach = 300.;
	double jump = (gamma_g+1.)/(gamma_g-1.);
	double F0, e_r0;
	double xx, yy;
	double f1, chi1;
	double f0, chi0;
	double tau;
	int i;

	jump = (gamma_g+1.)*Mach*Mach/((gamma_g-1.)*Mach*Mach+2.);

	y_up[0] = pdt.v_w-egn[1];
	y_up[1] = rho_csm(egn[3]);
	y[0] = y_up[0]/jump;
	y[1] = y_up[1]*jump;
	y[2] = y_up[0]*y_up[1]*(y_up[0]-y[0]);
	
	tau = func_tau_csm(*x);
	f0 = 1.0/fmax(tau, 1.);
	chi0 = func_chi_f(f0);
	yy = jump*jump*(1.0+chi0)*(1.0+chi0)/f0/f0;
	chi1 = (5.*yy-4.0-4.0*sqrt(yy*(yy-4.0)))/(3.0*yy+4.0);
	f1 = sqrt((-3.0+10.*chi1-3.0*chi1*chi1)/4.0);
	y[4] = egn[2];
	y[3] = chi1*(y[4]/((P_C)*f1));



	y[0] = log(-y[0]);
	y[1] = log(y[1]);

	for(i = 0; i < N; i++){
		y[i] /= yfor[i];
//		y[i] = 1.;
	}
//	*x = egn[3]/xfd;
}

/*
void boundary_for(double *x, double y[], double egn[])
{
	double chi, yy, chi0;
	double y_up[2];
	double gamma_g = 5./3.;
	double Mach = 300.;
	double jump = (gamma_g+1.)/(gamma_g-1.);
	double e_r0, e_r[2], F0;
	double tau;
	int i;

	jump = (gamma_g+1.)*Mach*Mach/((gamma_g-1.)*Mach*Mach+2.);

	y_up[0] = pdt.v_w-egn[1];
	y_up[1] = rho_csm(egn[3]);
	y[0] = y_up[0]/jump;
	y[1] = y_up[1]*jump;
	y[2] = y_up[0]*y_up[1]*(y_up[0]-y[0]);
	F0 = egn[2];
	e_r0 = F0/(P_C);
	yy = fabs(F0)/(P_C)/e_r0;
	chi0 = 1.;

	e_r[0] = 2.*e_r0*jump/(1.+1./3.);
	e_r[1] = e_r[0];
	do{
		e_r[0] = e_r[1];
		yy = fabs(F0)/(P_C)/e_r[1];
		chi = chi_f(yy);
		e_r[1] = (1.+chi0)*e_r0*jump/(1.+chi);
	}while(fabs(1.-e_r[0]/e_r[1]) > 1.e-14);
	y[0] = log(-y[0]);
	y[1] = log(y[1]);
	y[3] = e_r[0];
	y[4] = F0;
//	printf("%e %e %e %e %e T_gas = %e, T_rad = %e\n", y[0], y[1], y[2], y[3], y[4], y[2]/y[1]/(P_K)*(MU), pow(y[3]/(P_A), 0.25));
	for(i = 0; i < N; i++){
		yfor[i] = y[i];
		y[i] = 1.;
	}
	*x = egn[3]/xfd;
}
*/

void set_yrev(double *x, double y[], double egn[])
{
	double y_up[2];
	double gamma_g = 5./3.;
	double Mach = 300.;
	double jump = (gamma_g+1.)/(gamma_g-1.);
	double F0, e_r0, e_r1;
	double f0, f1, chi0, chi1;
	double tau;
	double xx, yy;
	int i;

	e_r0 = 3.*pdt.alpha*pdt.E_ej/(4.*M_PI*pdt.R_p*pdt.R_p*pdt.R_p)*pow(pdt.R_p/(*x), 4.);
	jump = (gamma_g+1.)*Mach*Mach/((gamma_g-1.)*Mach*Mach+2.);

	y_up[0] = (*x)/t_exp-egn[0];
	y_up[1] = rho_ej((*x), t_exp);
	y[0] = y_up[0]/jump;
	y[1] = y_up[1]*jump;
	y[2] = y_up[0]*y_up[1]*(y_up[0]-y[0]);
	tau = tau_ej((*x), t_exp);
	f0 = 1.0/tau;
	chi0 = func_chi_f(f0);
	yy = jump*jump*(1.0+chi0)*(1.0*chi0)/f0/f0;
	chi1 = (5.*yy-4.0-4.0*sqrt(yy*(yy-4.0)))/(3.0*yy+4.0);
	f1 = sqrt((-3.0+10.*chi1-3.0*chi1*chi1)/4.0);
	F0 = (P_C)*e_r0/tau;
	e_r1 = F0/((P_C)*f1);
	y[3] = chi1*e_r1;

	xx = 3.*jump;
	y[4] = F0;
	
//	printf("%.10e %.10e %.10e %.10e %.10e %e\n", y[0], y[1], y[2], y[3], y[4], y[4]/((P_C)*y[3]));
	
	y[0] = log(y[0]);
	y[1] = log(y[1]);

	for(i = 0; i < N; i++){
		yrev[i] = yrev[i];
//		y[i] = 1.;
	}
}

void set_yfor(double *x, double y[], double egn[])
{
	double y_up[2];
	double gamma_g = 5./3.;
	double Mach = 300.;
	double jump = (gamma_g+1.)/(gamma_g-1.);
	double F0, e_r0;
	double tau;
	double xx, yy;
	double f1, chi1;
	int i;

	jump = (gamma_g+1.)*Mach*Mach/((gamma_g-1.)*Mach*Mach+2.);

	y_up[0] = pdt.v_w-egn[1];
	y_up[1] = rho_csm(egn[3]);
	y[0] = y_up[0]/jump;
	y[1] = y_up[1]*jump;
	y[2] = y_up[0]*y_up[1]*(y_up[0]-y[0]);

	yy = 4.0*jump*jump;
	chi1 = (5.*yy-4.0-4.0*sqrt(yy*(yy-4.0)))/(3.0*yy+4.0);
	f1 = sqrt((-3.0+10.*chi1-3.0*chi1*chi1)/4.0);
	y[4] = egn[2];
	y[3] = chi1*(y[4]/((P_C)*f1));

	y[0] = log(-y[0]);
	y[1] = log(y[1]);

	for(i = 0; i < N; i++){
		yfor[i] = y[i];
	}
//	*x = egn[3]/xfd;
}

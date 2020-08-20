#include "diff_eq.h"
#include "main.h"

extern double t_exp;
extern opacity rosse, sct, pla;
extern double xfd, yrev[N], yfor[N];

//インプット、アウトプット共に全て無次元量

/*
void diff_eq_rev(double x, double y[], double dydx[])
{
	double j;
	double kappa_s, kappa, kappa_planck;
	double yl[N];
	double xl;
	double xx, yy;
	double T_g, T[2];
	double mu;
	double c_s;
	double R, p_r;
	double source;
	double rho, v;
	int i;
	for(i = 0; i < N; i++){
		yl[i] = yrev[i]*y[i];
	}
	v = exp(yl[0]);
	rho = exp(yl[1]);

	xl = x*xfd;
	T_g = T_gas(yl);
	yy = yl[4]/((P_C)*yl[3]);
	R = 12.*fabs(yy)/(2.-3.*fabs(yy)+sqrt(4.+12.*yy-15.*yy*yy));
	p_r = ((2.-3.*fabs(yy)+sqrt(4.+12.*fabs(yy)-15.*yy*yy))/12.+yy*yy)*yl[3];

	xx = p_r/yl[3];
	kappa = kappa_r(rho, T_g);
	kappa_s = sigma_sc(rho, T_g);
	kappa_planck = kappa_p(rho, T_g);
	j = j_ff(rho, T_g);
	j = kappa_planck*(P_C)*(P_A)*pow(T_g, 4.);
//	kappa_planck = j/((P_A)*(P_C)*pow(T_g, 4.));
	c_s = sqrt(5./3.*yl[2]/rho);
	source = j-(P_C)*kappa_planck*yl[3];
	
	dydx[0] = (2./3.*source+kappa*v*(yl[4])/(P_C)+2./xl*v*c_s*c_s)/(v*v-c_s*c_s)/v;
	dydx[1] = -2./xl*rho-dydx[0]*rho;
	dydx[1] /= rho;
	dydx[2] = kappa*rho*(yl[4])/(P_C)-rho*v*v*dydx[0];
	dydx[3] = -kappa*rho*yl[3]*yy/(xx-yy*yy);
	dydx[4] = -2./xl*yl[4]+rho*source;

	for(i = 0; i < N; i++){
		dydx[i] *= xfd/yrev[i];
	}
}
*/


void diff_eq_rev(double x, double y[], double dydx[])//y[0] = ln(v), y[1] = ln(rho), y[2] = p_g, y[3] = p_r, y[4] = F.
{
	double j;
	double kappa_s, kappa, kappa_planck;
	double yl[N];
	double xl;
	double T_g, T[2];
	double mu;
	double c_s;
	double source;
	double rho, v;
	double Erad;
	int i;
	for(i = 0; i < N; i++){
		yl[i] = yrev[i]*y[i];
	}
	Erad = func_Erad(yl);
	v = exp(yl[0]);
	rho = exp(yl[1]);

	xl = x*xfd;
	T_g = T_gas(yl);

	kappa = kappa_r(rho, T_g);
	kappa_s = sigma_sc(rho, T_g);
	kappa_planck = kappa_p(rho, T_g);
	j = j_ff(rho, T_g);
	j = kappa_planck*(P_C)*(P_A)*pow(T_g, 4.);
//	kappa_planck = j/((P_A)*(P_C)*pow(T_g, 4.));
	c_s = sqrt(5./3.*yl[2]/rho);
	source = j-(P_C)*kappa_planck*Erad;
	
	dydx[0] = (2./3.*source+kappa*v*(yl[4])/(P_C)+2./xl*v*c_s*c_s)/(v*v-c_s*c_s)/v;
	dydx[1] = -2./xl*rho-dydx[0]*rho;
	dydx[1] /= rho;
	dydx[2] = kappa*rho*(yl[4])/(P_C)-rho*v*v*dydx[0];
	dydx[3] = -kappa*rho*(yl[4])/(P_C)-(3.*yl[3]-Erad)/xl;
	dydx[4] = -2./xl*yl[4]+rho*source;

	for(i = 0; i < N; i++){
		dydx[i] *= xfd/yrev[i];
	}
}

void diff_eq_for(double x, double y[], double dydx[])//y[0] = ln(v), y[1] = ln(rho), y[2] = p_g, y[3] = p_r, y[4] = F.
{
	double j;
	double kappa_s, kappa, kappa_planck;
	double yl[N];
	double xl;
	double T_g, T[2];
	double mu;
	double c_s;
	double source;
	double rho, v;
	double Erad;
	int i;
	for(i = 0; i < N; i++){
		yl[i] = yfor[i]*y[i];
	}
	Erad = func_Erad(yl);
	v = -exp(yl[0]);
	rho = exp(yl[1]);

	xl = x*xfd;
	T_g = T_gas(yl);

	kappa = kappa_r(rho, T_g);
	kappa_s = sigma_sc(rho, T_g);
	kappa_planck = kappa_p(rho, T_g);
	j = j_ff(rho, T_g);
	j = kappa_planck*(P_C)*(P_A)*pow(T_g, 4.);
//	kappa_planck = j/((P_A)*(P_C)*pow(T_g, 4.));
	c_s = sqrt(5./3.*yl[2]/rho);
	source = j-(P_C)*kappa_planck*Erad;
	
	dydx[0] = (2./3.*source+kappa*v*(yl[4])/(P_C)+2./xl*v*c_s*c_s)/(v*v-c_s*c_s)/v;
	dydx[1] = -2./xl*rho-dydx[0]*rho;
	dydx[1] /= rho;
	dydx[2] = kappa*rho*(yl[4])/(P_C)-rho*v*v*dydx[0];
	dydx[3] = -kappa*rho*(yl[4])/(P_C)-(3.*yl[3]-Erad)/xl;
	dydx[4] = -2./xl*yl[4]+rho*source;

	for(i = 0; i < N; i++){
		dydx[i] *= xfd/yfor[i];
	}
}

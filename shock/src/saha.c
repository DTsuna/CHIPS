#include "saha.h"
#include "opacity.h"
#include "constant.h"

void saha(double rho, double U, double *mu, double *T)
{
	double GAMMA = 1.6666666666666666;
        double mu_tmp[2] = {};
        double x;
        double T_tmp = 2000.;
        double n_e, n_H, n_He, n_HI, n_HII, n_HeI, n_HeII, n_HeIII;
        mu_tmp[0] = 0.5;
        mu_tmp[1] = 1.;
        T_tmp = (GAMMA-1.)*mu_tmp[0]*U/P_K/rho;
	if(U > rho*3000./((GAMMA)-1.)*(P_K)/(1.3*(MH))){
        	while(fabs(mu_tmp[0]-mu_tmp[1]) > 1.e-13){
                	mu_tmp[0] = mu_tmp[1];
                	T_tmp = (GAMMA-1.)*mu_tmp[0]*MH*U/P_K/rho;
                	x = 2.*M_PI*P_E*P_K*T_tmp/(P_H*P_H);
                	n_H = X*rho/MH;
                	n_He = Y/4.*rho/MH;
                	n_e = rho/(mu_tmp[0]*MH)-(X+Y/4.)*rho/MH;
                	n_HII = pow(x, 1.5)*exp(-CHI_HI/(P_K*T_tmp))/(n_e+pow(x, 1.5)*exp(-CHI_HI/(P_K*T_tmp)))*n_H;
                	n_HI = n_H-n_HII;
                	n_HeI = pow(1.+4./n_e*pow(x, 1.5)*exp(-CHI_HeI/(P_K*T_tmp))+4./(n_e*n_e)*pow(x, 3.)*exp(-(CHI_HeI+CHI_HeII)/(P_K*T_tmp)), -1.)*n_He;
                	n_HeII = 4.*n_HeI/n_e*pow(x, 1.5)*exp(-CHI_HeI/(P_K*T_tmp));
                	n_HeIII = n_HeII/n_e*pow(x, 1.5)*exp(-CHI_HeII/(P_K*T_tmp));
                	mu_tmp[1] = pow(X*(1.+n_HII/n_H)+Y/4.*(1.+n_HeII/n_He+2.*n_HeIII/n_He), -1.);
                	mu_tmp[1] = (mu_tmp[1]+mu_tmp[0])/2.;
                	mu_tmp[1] = (mu_tmp[1]+mu_tmp[0])/2.;
                	mu_tmp[1] = (mu_tmp[1]+mu_tmp[0])/2.;
                	mu_tmp[1] = (mu_tmp[1]+mu_tmp[0])/2.;
        	}
		*mu = mu_tmp[0];
		*T = T_tmp;
	}
	else{
		*mu = pow(X*1.+(Y)/4., -1.);
		*T = (GAMMA-1.)*(*mu)*(MH)*U/(P_K)/rho;
	}
}

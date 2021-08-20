#include <stdio.h>
#include <math.h>
#include "opacity.h"
#include "diff_eq.h"
#include "constant.h"

void diff_eq(double x, double y[], double dydx[])
{
	double kappa;
	double dTdr[2] = {};
	double eps = 1.e-10;
	int count = 0, count_max = 100;
	kappa = kappa_r(y[1], y[2]);

	do{
		count++;
		if(count == count_max){
			printf("count max (diff_eq).\n");
			break;
		}
		dTdr[0] = dTdr[1];
		dTdr[1] = func_dTdr(y, kappa, dTdr[0]);
	}while(fabs(1.-dTdr[1]/dTdr[0]) > eps);

	dydx[2] = func_dTdr(y, kappa, dTdr[1]);
	dydx[1] = (-(4.0/3.0*(P_A)*pow(y[2], 3.0)+y[1]*(P_K)/((MU)*(MH)))*dydx[2]+2.0*y[1]*y[0]*y[0]/x)/(P_K*y[2]/((MU)*(MH))-y[0]*y[0]);
	dydx[0] = -2.0*y[0]/x-(y[0]/y[1])*dydx[1];
	dydx[3] = -2.0/x*(y[3]+0.5*y[1]*pow(y[0], 3.)+5.0/2.0*y[1]/((MU)*(MH))*(P_K)*y[2]*y[0]
		+(4.0/3.0)*(P_A)*pow(y[2], 4.)*y[0])-(1.5*y[1]*y[0]*y[0]+5.0/2.0*y[1]/((MU)*(MH))*(P_K)*y[2]
		+(4.0/3.0)*(P_A)*pow(y[2], 4.0))*dydx[0]-(0.5*pow(y[0], 3.0)+5.0/2.0*(P_K)/((MU)*(MH))*y[2]/y[0])*dydx[1]
		-(5.0/2.0*y[1]/((MU)*(MH))*(P_K)*y[0]+(16.0/3.0)*(P_A)*pow(y[2], 3.0)*y[0])*dydx[2];
}

double func_dTdr(double y[], double kappa, double dTdr)
{
	double R = 4.*fabs(dTdr/(kappa*y[1]*y[2]));
	double lambda = (2.+R)/(6.+3.*R+R*R);

	return -kappa*y[1]/(4.*(P_A)*(P_C)*pow(y[2], 3.)*lambda)*y[3];
}

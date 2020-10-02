#include <stdio.h>
#include <math.h>
#include "opacity.h"
#include "diff_eq.h"
#include "constant.h"

void diff_eq(double x, double y[], double dydx[])
{
	double kappa;
	kappa = kappa_r(y[1], y[2]);

	dydx[2] = -3.0*kappa*y[1]/(4.0*(P_A)*(P_C)*y[2]*y[2]*y[2])*y[3];
	dydx[1] = (-(4.0/3.0*(P_A)*pow(y[2], 3.0)+y[1]*(P_K)/((MU)*(MH)))*dydx[2]+2.0*y[1]*y[0]*y[0]/x)/(P_K*y[2]/((MU)*(MH))-y[0]*y[0]);
	dydx[0] = -2.0*y[0]/x-(y[0]/y[1])*dydx[1];
	dydx[3] = -2.0/x*(y[3]+0.5*y[1]*pow(y[0], 3.)+5.0/2.0*y[1]/((MU)*(MH))*(P_K)*y[2]*y[0]
		+(4.0/3.0)*(P_A)*pow(y[2], 4.)*y[0])-(1.5*y[1]*y[0]*y[0]+5.0/2.0*y[1]/((MU)*(MH))*(P_K)*y[2]
		+(4.0/3.0)*(P_A)*pow(y[2], 4.0))*dydx[0]-(0.5*pow(y[0], 3.0)+5.0/2.0*(P_K)/((MU)*(MH))*y[2]/y[0])*dydx[1]
		-(5.0/2.0*y[1]/((MU)*(MH))*(P_K)*y[0]+(16.0/3.0)*(P_A)*pow(y[2], 3.0)*y[0])*dydx[2];
}

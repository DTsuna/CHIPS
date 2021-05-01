#include <stdio.h>
#include <math.h>
#include "ranhugo.h"
#include "constant.h"

double func_mass(double y[])
{
	return y[0]*y[1];
}

double func_momentum(double y[])
{
	double Pg, Pr;
	Pg = y[1]*(P_K)*y[2]/((MU)*(MH));
	Pr = 1.0/3.0*(P_A)*pow(y[2], 4.0);

	return y[1]*y[0]*y[0]+Pg+Pr;
}

double func_energy(double y[])
{
	double Pg, Pr;
	Pg = y[1]*(P_K)*y[2]/((MU)*(MH));
	Pr = 1.0/3.0*(P_A)*pow(y[2], 4.0);

	return (0.5*y[1]*y[0]*y[0]+5.0/2.0*Pg+4.0*Pr)*y[0];
}

void set_func(double y_up[], double y_down[], double func[])
{
	int i;
	double (*func_rh[3])(double[]) = {func_mass, func_momentum, func_energy};

	for(i = 0; i < 3; i++){
		func[i] = (*func_rh[i])(y_down)-(*func_rh[i])(y_up);
	}
}

void jacobian_rh(double y[], double J[])
{
	double Pg, Pr;
	Pg = y[1]*(P_K)*y[2]/((MU)*(MH));
	Pr = 1.0/3.0*(P_A)*pow(y[2], 4.0);

	J[0] = y[1];
	J[1] = y[0];
	J[2] = 0.0;

	J[3] = 2.0*y[1]*y[0];
	J[4] = y[0]*y[0]+(P_K)*y[2]/((MU)*(MH));
	J[5] = y[1]*(P_K)/((MU)*(MH))+4.0/3.0*(P_A)*y[2]*y[2]*y[2];

	J[6] = (0.5*y[1]*y[0]*y[0]+5.0/2.0*Pg+4.0*Pr)+y[1]*y[0]*y[0];
	J[7] = (0.5*y[0]*y[0]+5.0/2.0*(P_K)*y[2]/((MU)*(MH)))*y[0];
	J[8] = (5.0/2.0*y[1]*(P_K)/((MU)*(MH))+16.*(P_A)*y[2]*y[2]*y[2])*y[0];
}

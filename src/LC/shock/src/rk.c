#include "rk.h"
#include <string.h>

/*This subroutine solves an ordinary differential equation by the 4th order (Runge-Kutta method)*/
void rk4(double y[], int m, double h, double x, double yout[], void(*dervs)(double, double[], double[]))
{
	double hh = 0.5*h, xh = x+hh, h6 = h/6.;
	double dym[2*m], dyt[2*m], yt[2*m+1];
	int i;
	memcpy(yt+m, y+m, sizeof(double)*m);
	yt[2*m] = x*y[2*m];
	(*dervs)(x, y, dym);
	for(i = 0; i < m; i++){
		yt[i] = y[i]+hh*dym[i];
	}
	(*dervs)(xh, yt, dyt);
	for(i = 0; i < m; i++){
		dym[i] += 2.*dyt[i];
		yt[i] = y[i]+hh*dyt[i];
	}
	(*dervs)(xh, yt, dyt);
	for(i = 0; i < m; i++){
		dym[i] += 2.*dyt[i];
		yt[i] = y[i]+h*dyt[i];
	}
	(*dervs)(x+h, yt, yout);
	for(i = 0; i < m; i++){
		yout[i] = yout[i]*h6;
		yout[i] = y[i]+h6*dym[i]+yout[i];
	}
}

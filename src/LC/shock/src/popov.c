#include "popov.h"
#include "constant.h"
#include <stdio.h>
#include <math.h>

double L_popov(double t, double Mej, double Eej, double R0, double X_H)
{
	double Tion = 5054., sigma_SB = (P_A)*(P_C)/4.0;
	double vsc, kappa;
	double te, td, ta, ti, Lambda;

	double L;

	kappa = 0.2*(1.+X_H);
	vsc = sqrt(10.*Eej/3./Mej);
	te = R0/vsc;
	td = 9.*kappa*Mej/(4.*pow(M_PI, 3.0)*(P_C)*R0);
	ta = sqrt(2.*td*te);
	Lambda = 2.*sigma_SB*81.*sqrt(10.)/pow(M_PI, 5.0)/(P_C)/(P_C)/sqrt(3.)*kappa*kappa*pow(Mej, 1.5)*pow(Tion, 4.0)/sqrt(Eej)/R0;
	ti = ta/sqrt(1.+Lambda);

	if(t > ti){
		L = 8.*M_PI*sigma_SB*pow(Tion, 4.0)*vsc*vsc*(ti*t*(1+ti*ti/3./ta/ta)-pow(t, 4.0)/3./ta/ta);
	}
	else{
		L = 0.5*Eej/td*exp(-t*t/ta/ta);
	}
	L = fmax(L, 0.0);
	
	return L;
}

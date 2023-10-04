#include "nickel.h"
#include "constant.h"
#include "pars.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


extern pars pdt;

//This subroutine is based on Valenti+2008, MNRAS, 383, 1485

double func_A(double z, double tau_ni, double tau_co, double tau_m)
{
	double y = tau_m/(2.*tau_ni);

	return 2.*z*exp(-2.*z*y+z*z);
}

double func_B(double z, double tau_ni, double tau_co, double tau_m)
{
	double y = tau_m/(2.*tau_ni);
	double s = tau_m*(tau_co-tau_ni)/(2.*tau_co*tau_ni);

	return 2.*z*exp(-2.*z*(y-s)+z*z);
}

double intg_using_simpson(double a, double b, double tau_ni, double tau_co, double tau_m, int n, double (*f)())
{
	double sum = 0., h;
	int i;

	h = (b-a)/(double)n/2.;

	sum = (*f)(a, tau_ni, tau_co, tau_m)+(*f)(b, tau_ni, tau_co, tau_m);
	for(i = 1; i < n; i++){
		sum += 4.*(*f)(a+(2.*(double)i-1.)*h, tau_ni, tau_co, tau_m)+2.*(*f)(a+2.*(double)i*h, tau_ni, tau_co, tau_m);
	}

	sum += 4.*(*f)(a+(2.*(double)n-1.)*h, tau_ni, tau_co, tau_m);
	sum *= h/3.;
	
	return sum;
}

double Lbol_before(double Mni, double Mej, double Eej, double kappa, double t)
{
	const double beta = 13.8;
	double tau_m;
	double x;
	double L;
	int n = 100;
	
	tau_m = sqrt(kappa/beta/(P_C))*pow(6.*Mej*Mej*Mej/5./Eej, 0.25);
	x = t/tau_m;
	L = Mni*exp(-x*x);
	L *= (EPS_NI-EPS_CO)*intg_using_simpson(0., x, TAU_NI, TAU_CO, tau_m, n, func_A)+EPS_CO*intg_using_simpson(0., x, TAU_NI, TAU_CO, tau_m, n, func_B);

	return L;
}

double Lbol_after(double Mni, double Mej, double Eej, double kappa, double t)
{
	const double beta = 13.8;
	double tau_m;
	double x;
	double L;
	double eps;
	double S_Ni_gamma, S_Co_gamma, S_Co_posi_gamma, S_Co_posi_ke;
	double one_minus_exp_F, one_minus_exp_G;
	double F, G;
	double C = (pdt.n-3.)*(pdt.n-3.)/(8.*M_PI*(pdt.n-1.)*(pdt.n-5.));
	double kappa_gam = 0.027;

	F = sqrt(C*kappa_gam*Mej*Mej/Eej);
	G = 16.1*F;
	eps = Mni*(EPS_CO)*(exp(-t/(TAU_CO))-exp(-t/(TAU_NI)));
	S_Ni_gamma = Mni*(EPS_NI)*exp(-t/(TAU_NI));

	tau_m = sqrt(kappa/beta/(P_C))*pow(6.*Mej*Mej*Mej/5./Eej, 0.25);
	x = t/tau_m;

	one_minus_exp_F = 1.-exp(-(F/t)*(F/t));
	one_minus_exp_G = 1.-exp(-(G/t)*(G/t));

	S_Co_gamma = 0.81*eps*one_minus_exp_F;
	S_Co_posi_gamma = 0.164*eps*one_minus_exp_F*one_minus_exp_G;
	S_Co_posi_ke = 0.036*eps*one_minus_exp_G;

	L = S_Ni_gamma+S_Co_gamma+S_Co_posi_gamma+S_Co_posi_ke;

	return L;
}



void search_smooth_point(double Mni, double Mej, double Eej, double kappa, double *tp)
{
	const double beta = 13.8;
	const double day = 86400.;
	double a, b, c;
	double tau_m;
	double x;
	double L;
	double t = day;
	double L1, L2;
	double tol = 1.e-06;
	int count = 0;
	
	tau_m = sqrt(kappa/beta/(P_C))*pow(6.*Mej*Mej*Mej/5./Eej, 0.25);
	x = t/tau_m;

	a = 0.01*tau_m;
	b = 10.*tau_m;
	if(fabs(Mni/(M_SUN)) > 0.00001){
		do{
			count++;
			c = (a+b)/2.;
			L1 = Lbol_before(Mni, Mej, Eej, kappa, a)/Lbol_after(Mni, Mej, Eej, kappa, a);
			L2 = Lbol_before(Mni, Mej, Eej, kappa, b)/Lbol_after(Mni, Mej, Eej, kappa, b);
			L = Lbol_before(Mni, Mej, Eej, kappa, c)/Lbol_after(Mni, Mej, Eej, kappa, c);
			if((1.-L1)*(1.-L) < 0.){
				b = c;
			}
			else if((1.-L2)*(1.-L) < 0.){
				a = c;
			}
			else{
				printf("no solution. count = %d.\n", count);
				exit(EXIT_FAILURE);
			}
		}while((b-a)/day > tol);
	}
	else{
		c = 0.;
	}

	*tp = c;
}


//return luminosity.
//This formula is valid until t = t_diff.
//If t > t_diff, this subroutine returns 0.
//Mni [g], Mej [g]
double rad_from_decay_of_nico(double Mni, double Mej, double Eej, double kappa, double t)
{
	const double beta = 13.8;
	double tau_m, L, point;
	
	tau_m = sqrt(kappa/beta/(P_C))*pow(6.*Mej*Mej*Mej/5./Eej, 0.25);

	search_smooth_point(Mni, Mej, Eej, kappa, &point);
	
	if(point > 0.0001){
		if(t < point){
			L = Lbol_before(Mni, Mej, Eej, kappa, t);
		}
		else{
			L = Lbol_after(Mni, Mej, Eej, kappa, t);
		}
	}
	else{
		L = 0.;
	}

	return L;
}

#include "function.h"

extern pars pdt;

double r_early(double t)
{
	double A, B;
	double n, s, delta, M_ej, E_ej, D;
	n = pdt.n;
	s = pdt.s;
	delta = pdt.delta;
	M_ej = pdt.M_ej;
	E_ej = pdt.E_ej;
	D = pdt.D;
	A = pow(2.*(5.-delta)*(n-5.)*E_ej, (n-3.)/2.);
	B = pow((3.-delta)*(n-3.)*M_ej, (n-5.)/2.);
	return pow((3.-s)*(4.-s)/(4.*M_PI*D*(n-4.)*(n-3.)*(n-delta))*(A/B), 1./(n-s))*pow(t, (n-3.)/(n-s));
}

double v_early(double t)
{
	return r_early(t)/t*(pdt.n-3.)/(pdt.n-pdt.s);
}

double rho_csm(double r)
{
	return pdt.D*pow(r, -pdt.s);
}

double rho_ej(double r, double t)
{
	double A, B;
	double n, delta, M_ej, E_ej, D;
	n = pdt.n;
	delta = pdt.delta;
	M_ej = pdt.M_ej;
	E_ej = pdt.E_ej;
	D = pdt.D;
	if(t > r/pdt.v_t){
		A = pow(2.0*(5.0-delta)*(n-5.0)*E_ej, (delta-3.)/2.);
		B = pow((3.0-delta)*(n-3.0)*M_ej, (delta-5.)/2.);
		return 1.0/(4.0*M_PI*(n-delta))*(A/B)*pow(t, delta-3.0)*pow(r, -delta);
	}
	else{
		A = pow(2.0*(5.0-delta)*(n-5.0)*E_ej, (n-3.)/2.);
		B = pow((3.0-delta)*(n-3.0)*M_ej, (n-5.)/2.);
		return 1.0/(4.0*M_PI*(n-delta))*(A/B)*pow(t, n-3.0)*pow(r, -n);
	}
}

double tau_ej(double r, double t)
{
	double A, B;
	double rhoej_coeff_d, rhoej_coeff_n;
	double n, delta, M_ej, E_ej, D, v_t;
	n = pdt.n;
	delta = pdt.delta;
	M_ej = pdt.M_ej;
	E_ej = pdt.E_ej;
	D = pdt.D;
	v_t =pdt.v_t;
	A = pow(2.0*(5.0-delta)*(n-5.0)*E_ej, (delta-3.)/2.);
	B = pow((3.0-delta)*(n-3.0)*M_ej, (delta-5.)/2.);
	rhoej_coeff_d = A/B/(4.*M_PI*(n-delta));
	A = pow(2.0*(5.0-delta)*(n-5.0)*E_ej, (n-3.)/2.);
	B = pow((3.0-delta)*(n-3.0)*M_ej, (n-5.)/2.);
	rhoej_coeff_n = A/B/(4.*M_PI*(n-delta));

	return rhoej_coeff_d*0.34*pow(t, -3.+delta)*pow(v_t*t, 1.-delta)/(1.-delta)
		+rhoej_coeff_n*0.34*(pow(r, -n+1.)-pow(v_t*t, -n+1.))/(1.-n)*pow(t, -3.+n);
}

double func_M_ej(double r, double t)
{
	double A, B;
	double n, delta, M_ej, E_ej, D;
	n = pdt.n;
	delta = pdt.delta;
	M_ej = pdt.M_ej;
	E_ej = pdt.E_ej;
	D = pdt.D;
	A = 2.0*(5.0-delta)*(n-5.0);
	B = (3.0-delta)*(n-3.0);
	if(r > pdt.v_t*t){
		return pow(t/r,n-3.0)*pow(A*E_ej, (n-3.0)/2.0)/pow(B*M_ej, (n-5.0)/2.0)/((n-delta)*(n-3.0));
	}
	else{
		return pow(1./pdt.v_t, n-3.0)*pow(A*E_ej, (n-3.0)/2.0)/pow(B*M_ej, (n-5.0)/2.0)/((n-delta)*(n-3.0))
		+(pow(A*E_ej, (delta-3.)/2.)/pow(B*M_ej, (delta-5.)/2.))*(pow(pdt.v_t, -delta+3.)-pow(r/t, -delta+3.))/((n-delta)*(n-3.));
	}
}

double T_gas(double y[])
{
	double mu[2], T[2];
	double rho;
	int count = 0, count_max = 100;
	rho = exp(y[1]);
	T[0] = (MU)/(P_K)*y[2]/rho;
	do{
		count++;
		T[1] = T[0];
		mu[0] = mmw(rho, T[1]);
		T[0] = (mu[0]*(MH))/(P_K)*y[2]/rho;
//		T[0] = (T[0]+10.*T[1])/11.;
		mu[1] = mmw(rho, T[0]);
		T[0] = (mu[0]+mu[1])/2.*(MH)/(P_K)*y[2]/rho;
	}while(fabs(1.-T[0]/T[1]) > 1.e-09 && count < count_max);
	if(count == count_max){
//		printf("count = %d\n", count);
	}
	return T[0];
}

double func_chi(double yl[])
{
	double a = fabs(yl[4]/((P_C)*yl[3]));
	if(a <= 1.){
		return 3./(5.-2.*a+2.*sqrt((1.-a)*(4.-a)));
	}
	else{
//		printf("imag number\n");
		return 1.;
	}
}

double func_Erad(double yl[])
{
	double chi = func_chi(yl);
	return yl[3]/chi;
}

double T_rad(double yl[])
{
	double Erad = func_Erad(yl);
	return pow(Erad/(P_A), 0.25);
}

double func_chi_f(double f)
{
	return (3.0+4.0*f*f)/(5.0+2.0*sqrt(4.0-3.0*f*f));
}

double func_tau_csm(double r)
{
	double rho = rho_csm(r), kappa_es = 0.2*(1.+X);

	return 1./(pdt.s-1.)*kappa_es*rho*r;
}

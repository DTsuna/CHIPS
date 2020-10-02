#include <stdio.h>
#include "pars.h"

double t_exp;
pars pdt;

//substitute parameters to struct **pars**.
pars setpars(double n, double delta, double E_ej, double M_ej, double v_w, double t_ini, double R_p)
{
	pars pdt;

	pdt.n = n;
	pdt.delta = delta;
	pdt.E_ej = E_ej;
	pdt.M_ej = M_ej;
	pdt.v_w = v_w;
	pdt.t_ini = t_ini;
	pdt.R_p = R_p;
	pdt.v_t = sqrt(2.0*(5.0-delta)*(n-5.0)/((3.0-delta)*(n-3.0))*E_ej/M_ej);

	return pdt;
}

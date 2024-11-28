#include <stdio.h>
#include <math.h>
#include "pars.h"
#include "constant.h"

double t_exp;
pars pdt;

//substitute parameters to struct **pars**.
pars setpars(double n, double delta, double E_ej, double M_ej, double M_ni, double s, double v_w, double R_p)
{
	pars pdt;

	pdt.n = n;
	pdt.s = s;
	pdt.delta = delta;
	pdt.E_ej = E_ej;
	pdt.M_ej = M_ej*(M_SUN);
	pdt.M_ni = M_ni*(M_SUN);
	pdt.v_w = v_w;
	pdt.R_p = R_p;
	pdt.v_t = sqrt(2.0*(5.0-delta)*(n-5.0)/((3.0-delta)*(n-3.0))*E_ej/pdt.M_ej);

	return pdt;
}

#ifndef INCLUDED_PARS_H
#define INCLUDED_PARS_H

typedef struct{
	double n;
	double delta;
	double s;
	double E_ej;
	double M_ej;
	double M_ni;
	double v_w;
	double R_p;
	double v_t;
}pars;

pars setpars(double, double, double, double, double, double, double, double);

#endif

#ifndef INCLUDED_FUNCTION_H
#define INCLUDED_FUNCTION_H

double r_early(double);
double t_early(double);
double rho_csm(double);
double rho_ej(double, double);
double func_M_ej(double, double);
double set_r_ini(const char*);
double set_r_diff(const char*);
double func_M_csm(double, double);
double p_tot(double[]);
double v_wind(double);

#endif

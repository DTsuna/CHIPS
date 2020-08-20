#ifndef INCLUDED_FUNCTION_H
#define INCLUDED_FUNCTION_H

#include <stdio.h>
#include <math.h>
#include "main.h"
#include "opacity.h"

double r_early(double);
double v_early(double);
double rho_csm(double);
double rho_ej(double, double);
double tau_ej(double, double);
double func_M_ej(double, double);
double T_gas(double[]);
double T_rad(double[]);
double func_chi(double[]);
double func_Erad(double[]);
double func_chi_f(double);
double func_tau_csm(double);

#endif

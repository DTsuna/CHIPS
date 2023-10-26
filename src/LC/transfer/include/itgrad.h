#ifndef INCLUDED_ITGRAD_H
#define INCLUDED_ITGRAD_H

#define CFL_MAX 0.9
#define CFL_MIN 1.0E-4
#define R_MAX 1.E5

void integ_adv_eadtra(double F_ini, double E[], double T_g[], double r[], double dt, const int n, int *flag);
double output_R(double E[], double T_g[], double r[], int i);
double Flux(double F_ini, double E[], double T_g[], double r[], int i, const int n);
double func_radtra(double F_ini, double E[], double T_g[], double r[], double dt, int i, const int n);
void Jacob_func_radtra(double F_ini, double E[], double T_g[], double r[], double a[], double b[], double c[], double func[], double dt, const int n);

#endif

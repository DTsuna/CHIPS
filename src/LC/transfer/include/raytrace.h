#ifndef INCLUDED_RAYTRACE_H
#define INCLUDED_RAYTRACE_H

#include "opacity.h"
#include "pars.h"

int imax(int a, int b);
double ds_path(double b, double r[], int j);
double alpha_ff_nu(double nu, double rho, double T);
double beta_nu(double nu, double rho, double T);
double alpha_nu(double nu, double rho, double T, opacity op);
double Planck_func(double nu, double T);
double integ_ray_tracing(double b, double nu, double *tau_nu, double F_ini, double r[], double rho[], double T[], int n, opacity op);
double Lum_nu(double r_init, double r_out, double nu, double F_ini, double r[], double rho[], double T[], int n, char *filenamenu, int nnu);
void calc_lum(double t, double r_init, double r_out, double F_ini, double r[], double rho[], double T[], int n, char *filename, double abmag[], pars pdt);
int jmin_func(double b, double r[], int n);

#endif

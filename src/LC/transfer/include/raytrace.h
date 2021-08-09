#ifndef INCLUDED_RAYTRACE_H
#define INCLUDED_RAYTRACE_H

#include "opacity.h"

int imax(int a, int b);
double ds_path(double b, double r[], int j);
double alpha_ff_nu(double nu, double rho, double T);
double alpha_p_beta_nu(double nu, double rho, double T, opacity op);
double Planck_func(double nu, double T);
double integ_ray_tracing(double b, double nu, double r[], double rho[], double T[], double r_sh[], double rho_sh[], double T_sh[], int n, int n_sh, opacity op);
double Lum_nu(double r_init, double r_out, double nu, double r[], double rho[], double T[], double r_sh[], double rho_sh[], double T_sh[], int n, int n_sh);
void calc_lum(double r_init, double r_out, double r[], double rho[], double T[], double r_sh[], 
	double rho_sh[], double T_sh[], int n, int n_sh, char *filename, double abmag[]);
int jmin_func(double b, double r[], int n);

#endif

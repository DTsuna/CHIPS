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
void interp_self_similar_values(double *A, double *E_rev, double *E_for);
void set_abundance(void);
void gen_opacity_sc(void);
double sigma_saha(double R, double T);
void sigma_mu_saha(double R, double T, double *sigma, double *mu);

#endif

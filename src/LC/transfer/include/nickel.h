#ifndef INCLUDED_NICKEL_H
#define INCLUDED_NICKEL_H

#define TAU_NI 7.605E+05
#define TAU_CO 9.822E+06
#define EPS_NI 3.9E+10
#define EPS_CO 6.78E+09

double func_A(double z, double tau_ni, double tau_co, double tau_m);
double func_B(double z, double tau_ni, double tau_co, double tau_m);
double integ_using_simpson(double a, double b, double tau_ni, double tau_co, double tau_m, int n, double (*f)());
double Lbol_before(double Mni, double Mej, double Eej, double kappa, double t);
double Lbol_after(double Mni, double Mej, double Eej, double kappa, double t);
void search_smooth_point(double Mni, double Mej, double Eej, double kappa, double *tp);
double rad_from_decay_of_nico(double Mni, double Mej, double Eej, double kappa, double t);

#endif

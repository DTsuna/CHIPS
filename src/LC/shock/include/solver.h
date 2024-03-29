#ifndef INCLUDED_SOLVER_H
#define INCLUDED_SOLVER_H

double func_dM(double, double, double, double);
void solver_rev(double, double[], double[], int*);
void solver_for(double[], double[], double[], int, int*);
void solver(double, double[], double[], int, int*);
void set_phys(double, double, double[], double[], double[], int);
void solver_rev_outp(double x_ini, double int_phys[], double egn[], int *info, FILE *fp);
void solver_for_outp(double int_phys[], double ext_phys[], double egn[], int flag, int *info, FILE *fp);
void solver_outp(double x_ini, double phys[], double egn[], int flag, int *info, FILE *fp);

#endif

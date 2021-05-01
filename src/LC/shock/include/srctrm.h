#ifndef INCLUDED_SRCTRM_H
#define INCLUDED_SRCTRM_H

double func_rhs(double, double, double);
void diff_eq_src(double[], double[], double, double, double[]);
void jacob(double, double, double, double, double[]);
void itg_src(double[], double[], double, double, double);

#endif

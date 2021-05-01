#ifndef INCLUDED_ITGADVE_H
#define INCLUDED_ITGADVE_H

double func_lambda(double);
void matrix_E(double, double[], double[], double[], double[], double, double[], double[], double[], const int);
void init_E(double, double, double[], double[], double, double[], const int);
void itg_adv_E(double, double, double[], double[], double[], double[], double, const int);

#endif

#ifndef INCLUDED_SOLVER_H
#define INCLUDED_SOLVER_H

double func_dM(double, double, double, double);
void solver_rev(double, double[], double[]);
void solver_for(double[], double[], double[], int);
void solver(double, double[], double[], int);
void set_phys(double, double, double[], double[], double[], int);

#endif

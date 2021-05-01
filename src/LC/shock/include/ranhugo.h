#ifndef INCLUDED_RANHUGO_H
#define INCLUDED_RANHUGO_H

double func_mass(double[]);
double func_momentum(double[]);
double func_energy(double[]);
void set_func(double[], double[], double[]);
void jacobian_rh(double[], double[]);

#endif

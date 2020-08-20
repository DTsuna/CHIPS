#ifndef INCLUDED_SOLVER_NEQ_H
#define INCLUDED_SOLVER_NEQ_H

void solver_rev(double, double[], double[], FILE*);
void solver_for(double[], double[], double[], FILE*);
void solver(double, double[], double[], FILE*, FILE*);
void sub_Mass_shocked(double, double, double[], double[], double[], double[], double, void(*)(double, double[], double[]));
void dimless2dim(double[], double[], double[]);
void dim2dimless(double[], double[], double[]);
void set_phys_cd(double, double[], double[], double, double[]);
double set_final_dx(double, double[], double[], double, double, double, double[]);
double func_dMass(double, double, double[], double[], double[]);
void jump_calc(double[], double[], double[]);
double pressure_rad(double[]);

#endif

#ifndef INCLUDED_SOLVER_H
#define INCLUDED_SOLVER_H

#include <stdio.h>
#include "main.h"
#include "irk.h"
#include "diff_eq.h"
#include "boundary.h"

void solver_rev(double, double[], double[]);
void solver_rev_outp(double, double[], double[], FILE*);
void solver_for(double[], double[], double[]);
void solver_f(double[], double[], double[], FILE*);
void solver_for_outp(double[], double[], double[], FILE*);
void err_itr(double*, double, double[], double[], double[], int*, void(*)(double, double[], double[]));
void jump_calc(double[], double[], double[]);
void dimless2dim(double[], double[], double[]);
void dim2dimless(double[], double[], double[]);
double dx_estim(double, double[], double[], double, void(*)(double, double[], double[]));
double ana_sol_kconst_rho(double, double, double[]);
double ana_sol_kconst_pgas(double, double, double[]);
double ana_sol_kconst_v(double, double, double[]);
double ana_sol_kconst_erad(double, double, double[]);
double ana_sol_kconst_flux(double, double, double[]);
double ana_sol_kprho_rho(double, double, double[]);
double ana_sol_kprho_pgas(double, double, double[]);
double ana_sol_kprho_v(double, double, double[]);
double ana_sol_kprho_erad(double, double, double[]);
double ana_sol_kprho_flux(double, double, double[]);

#endif

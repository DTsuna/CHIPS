#ifndef INCLUDED_IRK_H
#define INCLUDED_IRK_H

#define N_24 2
#define N_36 3

#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include "W4.h"

double a_24[N_24*N_24];
double b_24[N_24], c_24[N_24];

double a_36[N_36*N_36];
double b_36[N_36], c_36[N_36];

void func_setup(double[], double[], int, int, double, double, double[], double[], double[], void(*)(double, double[], double[]));
void irk24(double[], int, double, double, double[], double, double*, void(*)(double, double[], double[]));
void irk36(double[], int, double, double, double[], double, double*, void(*)(double, double[], double[]));

#endif

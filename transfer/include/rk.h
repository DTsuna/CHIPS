#ifndef INCLUDED_RK_H
#define INCLUDED_RK_H

#include <math.h>
#include <stdio.h>

void rk4(double[], int, double, double, double[], void(*)(double, double[], double[]));

#endif 

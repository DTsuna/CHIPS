#ifndef INCLUDED_W4_H
#define INCLUDED_W4_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

void get_UL(double*, double*, double*, int);
void get_U_inv(double*, double*, int);
void get_L_inv(double*, double*, int);
void get_x_p(double*, double*, double*, double*, double*, double*, double*, double, int);
void get_dx_dp(double*, double*, double*, double*, double*, double*, double*, double, int);
void get_itr_x(double*, double*, double*, double*, double, int);
void get_itr_x_tol(double*, double*, double*, double*, double, double, int);

#endif

#ifndef INCLUDED_OPACITY_H
#define INCLUDED_OPACITY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dcht.h"
#include "main.h"

#define X 0.737
#define Y 0.249

typedef struct{
	double T[128];
	double R[128];
	double kappa[4096*4];
	int imax;
	int jmax;
}opacity;

void set_opacity(const char*, opacity*);
double kappa_r(double, double);
double sigma_sc(double, double);
double kappa_p(double, double);
double mmw(double, double);
double sigma_sc_high(double, double);
double kappa_ff(double, double);
double j_ff(double, double);

#endif

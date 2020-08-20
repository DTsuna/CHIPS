#ifndef INCLUDED_BNDR_H
#define INCLUDED_BNDR_H

#include <stdio.h>
#include <math.h>
#include "function.h"
#include "main.h"

void boundary_rev(double*, double[], double[]);
void boundary_for(double*, double[], double[]);
void set_yrev(double*, double[], double[]);
void set_yfor(double*, double[], double[]);

#endif

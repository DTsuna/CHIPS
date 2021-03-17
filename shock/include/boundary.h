#ifndef INCLUDED_BOUNDARY_H
#define INCLUDED_BOUNDARY_H

void set_init_rev(double, double[], double[]);
void set_init_for(double, double[], double[]);
void set_init_down_rev(double, double[], double[], double[]);
void set_init_down_for(double, double[], double[], double[]);
void boundary(double, double[], double[], int, int*);

#endif

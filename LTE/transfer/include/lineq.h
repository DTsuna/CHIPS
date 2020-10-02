#ifndef INCLUDED_LINEQ_H
#define INCLUDED_LINEQ_H

#include <float.h>

void mat_gauss(double*, double*, int);
void get_PLU(double*, int*, int);
void mat_PLU(double*, double*, int*, int);
void mat_inv(double*, int);

void copy_mat(double*, double*, int, int);
void copy_mat_t(double*, double*, int, int);
void print_mat(char*, double*, int, int);
void print_mat_z(char*, double*, int, int);
void print_vec(char*, double*, int);
void calc_ax(double*, double*, double*, int);
void product_mat(double*, double*, double*, int, int, int);

#endif

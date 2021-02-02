#ifndef INCLUDED_ITGADVU_H
#define INCLUDED_ITGADVU_H

void init_U(double[], double[], const int);
void matrix_U(double[], double[], double[], double[], double[], double, double[], double[], double[], const int);
void itg_adv_U(double[], double[], double[], double[], double[], double, const int);

#endif

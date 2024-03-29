#ifndef INCLUDED_SAHA_H
#define INCLUDED_SAHA_H

#define CHI_HI 2.16935E-11
#define CHI_HeI 3.92213E-11
#define CHI_HeII 8.67901E-11

void saha(double, double, double*, double*);
void get_num_density(double rho, double T,  double ndens[]);
void set_tables(const char *openfile, double rho[], double U[], double mu[], int *row, int *col);

#endif

#ifndef INCLUDED_SAHA_WR_H
#define INCLUDED_SAHA_WR_H

void get_Jacobi_wr(double rho, double Tem, double n_ion[], double fHeofT, double fHeIIofT, double fCofT, double fCIIofT, double fCIIIofT, double fCIVofT, double fOofT,double fOIIofT, double fOIIIofT, double fOIVofT, double n_e, double Jac[], double X_He, double X_C, double X_O);
void Saha_wr_U(double rho, double U, double *mu, double *Tem);
double set_n_el(double n_ion[]);
void set_f_ion(double rho, double Tem, double f_ion[], double n_ion[], double fHeofT, double fHeIIofT, double fCofT, double fCIIofT, double fCIIIofT, double fCIVofT, double fOofT,double fOIIofT, double fOIIIofT, double fOIVofT, double X_He, double X_C, double X_O);
void Saha_U(double rho, double U, double *mu, double *T);

#endif

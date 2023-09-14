#include "saha_wr.h"
#include "constant.h"
#include "lineq.h"
#include "saha.h"
#include "opacity.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

// input parameters
// double X_He=0.28, X_C=0.5, X_O=0.2; // mass fraction of helium, carbon and oxygen

// constant parameters

void gen_muTem_table(int discriminant)
{
	FILE *fp;
	const int row = 21;
	const int col = 19;
	double rho_array[row];
	double U_array[col];
	double mu, T;
	double mu_array[row*col], T_array[row*col];
	int i, j;

	for(i = 0; i < row; i++){
		rho_array[i] = pow(10., -20.+0.5*(double)i);
	}

	for(i = 0; i < col; i++){
		U_array[i] = pow(10., 10.5+0.25*(double)i);
	}

	for(i = 0; i < row; i++){
		for(j = 0; j < col; j++){
			if(discriminant == 2){
				Saha_wr_U(rho_array[i], U_array[j], &mu, &T);
			}
			else if(discriminant == 0 || discriminant == 1){
				Saha_U(rho_array[i], U_array[j], &mu, &T);
			}
			else{
				printf("Cannot generate tables for mu, T. exit.\n");
				exit(EXIT_FAILURE);
			}
			mu_array[i*col+j] = mu;
			T_array[i*col+j] = T;
//			printf("mu = %f, T = %e K, rho = %e g/cc, U = %e \n", mu, T, rho_array[i], U_array[j]);
		}
	}

	fp = fopen("./LCFiles/mu.txt", "w");
	for(i = 0; i < col; i++){
		fprintf(fp, "%.2f", log10(U_array[i]));
		if(i < col-1){
			fprintf(fp, " ");
		}
	}
	fprintf(fp, "\n");
	for(i = 0; i < row; i++){
		fprintf(fp, "%.2f ", log10(rho_array[i]));
		for(j = 0; j < col; j++){
			fprintf(fp, "%.4f", mu_array[i*col+j]);
			if(j < col-1){
				fprintf(fp, " ");
			}
		}
		fprintf(fp, "\n");
	}

	fclose(fp);

	fp = fopen("./LCFiles/Tem.txt", "w");
	for(i = 0; i < col; i++){
		fprintf(fp, "%.2f", log10(U_array[i]));
		if(i < col-1){
			fprintf(fp, " ");
		}
	}
	fprintf(fp, "\n");
	for(i = 0; i < row; i++){
		fprintf(fp, "%.2f ", log10(rho_array[i]));
		for(j = 0; j < col; j++){
			fprintf(fp, "%.4f", log10(T_array[i*col+j]));
			if(j < col-1){
				fprintf(fp, " ");
			}
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}


void Saha_wr_U(double rho, double U, double *mu_outp, double *T_outp)
{
// determine electron density by solving Saha's eq for He, C, O    !!!!!
	FILE *fp;
	int i_ion, j_ion;
	int pvt[100];
	double n_el;
	double fofT, fCofT, fHeofT, fOofT; //f(T) = (2.*pi*m_e*k_B*T/h**2)^(3/2), fA(T)=f(T)*exp(-xi_A/k_B/T)
	double fCIIofT, fHeIIofT, fOIIofT, fCIIIofT, fOIIIofT, fCIVofT, fOIVofT;
	double Jac[100];
	double n_HeII, n_CII, n_OII, n_HeII_new, n_CII_new, n_OII_new; // 1st ionization state
	double n_HeIII, n_CIII, n_OIII, n_HeIII_new, n_CIII_new, n_OIII_new; // 2nd ionization
	double n_CIV, n_OIV, n_CIV_new, n_OIV_new; // 3rd ionization
	double n_CV, n_OV, n_CV_new, n_OV_new; // 4th ionization
	double n_ion[10], n_ion_new[10], n_ion_old[10], f_ion[10];
	double eps_Saha, max_diff, max_density; // for convergence check
	double eps, tol = 1.e-8;
	double Tem, mu[2];
	double n_He, n_C, n_O;
	double X_H, X_He, X_C, X_O, X_dummy[4];
	int i = 0, count = 0;

	fp = fopen("input/abundance/abundance_for_tablegen.txt", "r");
	while(fscanf(fp, "%lf", X_dummy+i) != EOF){
		i++;
	}

	X_He = X_dummy[1];
	X_C = X_dummy[2];
	X_O = X_dummy[3];

	n_He = X_He*rho/(4.*(MH));
	n_C = X_C*rho/(12.*(MH));
	n_O = X_O*rho/(16.*(MH));
	
// initialize ion densities as vectors
// in order, HeII, HeIII, CII, ..., CV, OII, ... OV
	for(i_ion = 0; i_ion < 10; i_ion++){
		n_ion[i_ion] = rho/(MH)/5.;
	}
	mu[0] = X_He/4.*(1.+n_ion[0]+2.*n_ion[1])/n_He
		+X_C/12.*(1.+n_ion[2]+2.*n_ion[3]+3.*n_ion[4]+4.*n_ion[5])/n_C
		+X_O/16.*(1.+n_ion[6]+2.*n_ion[7]+3.*n_ion[8]+4.*n_ion[9])/n_O;
	mu[0] = 1./mu[0];
	mu[1] = mu[0];

	eps_Saha = 1.e0;

	while(eps_Saha > 1.0e-12){
		count++;
		mu[0] = mu[1];
		mu[1] = X_He/4.*(1.+n_ion[0]+2.*n_ion[1])/n_He
			+X_C/12.*(1.+n_ion[2]+2.*n_ion[3]+3.*n_ion[4]+4.*n_ion[5])/n_C
			+X_O/16.*(1.+n_ion[6]+2.*n_ion[7]+3.*n_ion[8]+4.*n_ion[9])/n_O;
		mu[1] = 1./mu[1];
		mu[1] = (9.*mu[0]+mu[1])/10.;
		Tem = 2./3.*mu[1]*(MH)/(P_K)*U;
//		printf("eps = %e, mu = %f, T = %e K.\n", eps, mu[1], Tem);
// temperature-dependent part
		fofT = (2.*M_PI*(P_E)*(P_K)*Tem/(P_H)/(P_H))*sqrt(2.*M_PI*(P_E)*(P_K)*Tem/(P_H)/(P_H));
		fHeofT = fofT * exp(-(XI_He)*1.602e-12/((P_K)*Tem)); //1.6d-12: eV to erg
		fHeIIofT = fofT * exp(-(XI_HeII)*1.602e-12/((P_K)*Tem));
		fCofT = fofT * exp(-(XI_C)*1.602e-12/((P_K)*Tem));
		fCIIofT = fofT * exp(-(XI_CII)*1.602e-12/((P_K)*Tem));
		fCIIIofT = fofT * exp(-(XI_CIII)*1.602e-12/((P_K)*Tem));
		fCIVofT = fofT * exp(-(XI_CIV)*1.602e-12/((P_K)*Tem));
		fOofT = fofT * exp(-(XI_O)*1.602e-12/((P_K)*Tem));
		fOIIofT = fofT * exp(-(XI_OII)*1.602e-12/((P_K)*Tem));
		fOIIIofT = fofT * exp(-(XI_OIII)*1.602e-12/((P_K)*Tem));
		fOIVofT = fofT * exp(-(XI_OIV)*1.602e-12/((P_K)*Tem));
// electron density
		n_el = set_n_el(n_ion);
// number density dependent part, we set this as a vector
// in order, HeII, HeIII, CII, ..., CV, OII, ... OV
		set_f_ion(rho, Tem, f_ion, n_ion, fHeofT, fHeIIofT, fCofT, fCIIofT, fCIIIofT, fCIVofT, fOofT, fOIIofT, fOIIIofT, fOIVofT, X_He, X_C, X_O);
// calculate Jacobian matrix and its inverse
		get_Jacobi_wr(rho, Tem, n_ion, fHeofT, fHeIIofT, fCofT, fCIIofT, fCIIIofT, fCIVofT, fOofT, fOIIofT, fOIIIofT, fOIVofT, n_el, Jac, X_He, X_C, X_O);
		get_PLU(Jac, pvt, 10);
		mat_PLU(Jac, f_ion, pvt, 10);
		for(i_ion = 0; i_ion < 10; i_ion++){
			n_ion_old[i_ion] = n_ion[i_ion];
			n_ion[i_ion] = n_ion[i_ion]-f_ion[i_ion];
			n_ion_new[i_ion] = n_ion[i_ion];
		}

// multiply the inverse matrix and fionofn to get ionofn_new
	
//error is estimated so that the most dominant value is obtained properly
		for(i_ion = 0; i_ion < 10; i_ion++){
			if(i_ion == 0){
				max_diff = fabs(n_ion_new[0]-n_ion_old[0]);
				max_density = n_ion_new[0];
			}
			else{
				max_diff = fmax(max_diff, fabs(n_ion_new[i_ion]-n_ion_old[i_ion]));
				max_density = fmax(max_density, n_ion_new[i_ion]);
			}
		}
		eps_Saha = max_diff/max_density;
//		printf("eps_Saha = %e\n", eps_Saha);
	}

	for(i_ion = 0; i_ion < 10; i_ion++){
		n_ion[i_ion] = n_ion_new[i_ion];
	}
	n_el = set_n_el(n_ion);
	*T_outp = Tem;
	*mu_outp = mu[1];

	fclose(fp);
}


double set_n_el(double n_ion[])
{
	return (n_ion[0]+2.*n_ion[1])+(n_ion[2]+2.*n_ion[3]+3.*n_ion[4]+4.*n_ion[5])+(n_ion[6]+2.*n_ion[7]+3.*n_ion[8]+4.*n_ion[9]);
}


void set_f_ion(double rho, double Tem, double f_ion[], double n_ion[], double fHeofT, double fHeIIofT, double fCofT, double fCIIofT, double fCIIIofT, double fCIVofT, double fOofT,double fOIIofT, double fOIIIofT, double fOIVofT, double X_He, double X_C, double X_O)
{
	double n_el;

	n_el = set_n_el(n_ion);

	f_ion[0] = n_ion[0]*n_el-(X_He*rho/(4.*(MH))-n_ion[0]-n_ion[1])*(4.*fHeofT);
	f_ion[1] = n_ion[1]*n_el-n_ion[0]*fHeIIofT;
	f_ion[2] = n_ion[2]*n_el-(X_C*rho/(12.*(MH))-n_ion[2]-n_ion[3]-n_ion[4]-n_ion[5])*(4.*fCofT);
	f_ion[3] = n_ion[3]*n_el-n_ion[2]*(1.*fCIIofT);
	f_ion[4] = n_ion[4]*n_el-n_ion[3]*(4.*fCIIIofT);
	f_ion[5] = n_ion[5]*n_el-n_ion[4]*(1.*fCIVofT);
	f_ion[6] = n_ion[6]*n_el-(X_O*rho/(16.*(MH))-n_ion[6]-n_ion[7]-n_ion[8]-n_ion[9])*(1.6*fOofT);
	f_ion[7] = n_ion[7]*n_el-n_ion[6]*(0.5*fOIIofT);
	f_ion[8] = n_ion[8]*n_el-n_ion[7]*(4.*fOIIIofT);
	f_ion[9] = n_ion[9]*n_el-n_ion[8]*(1.*fOIVofT);
}

void get_Jacobi_wr(double rho, double Tem, double n_ion[], double fHeofT, double fHeIIofT, double fCofT, double fCIIofT, double fCIIIofT, double fCIVofT, double fOofT,double fOIIofT, double fOIIIofT, double fOIVofT, double n_e, double Jac[], double X_He, double X_C, double X_O)
{
	double f_ion[10], f_ion_1[10];
	double dn_ion[10];
	int i, j;

	for(i = 0; i < 10; i++){
		dn_ion[i] = 1.e-06*n_ion[i];
	}

	set_f_ion(rho, Tem, f_ion, n_ion, fHeofT, fHeIIofT, fCofT, fCIIofT, fCIIIofT, fCIVofT, fOofT, fOIIofT, fOIIIofT, fOIVofT, X_He, X_C, X_O);
	for(j = 0; j < 10; j++){
		n_ion[j] += dn_ion[j];
		set_f_ion(rho, Tem, f_ion_1, n_ion, fHeofT, fHeIIofT, fCofT, fCIIofT, fCIIIofT, fCIVofT, fOofT, fOIIofT, fOIIIofT, fOIVofT, X_He, X_C, X_O);
		for(i = 0; i < 10; i++){
			Jac[10*i+j] = (f_ion_1[i]-f_ion[i])/dn_ion[j];
		}
		n_ion[j] -= dn_ion[j];
	}
}

void Saha_U(double rho, double U, double *mu, double *T)
{
	FILE *fp;
	double GAMMA = 1.6666666666666666;
        double mu_tmp[2] = {};
        double x;
        double T_tmp = 2000.;
        double n_e, n_H, n_He, n_HI, n_HII, n_HeI, n_HeII, n_HeIII;
	double X_H, X_He, X_dummy[4];
	int i = 0, count = 0;

	fp = fopen("input/abundance/abundance_for_tablegen.txt", "r");
	while(fscanf(fp, "%lf", X_dummy+i) != EOF){
		i++;
	}

	X_H = X_dummy[0];
	X_He = X_dummy[1];

	n_H = X_H*rho/(MH);
	n_He = X_He*rho/(4.*(MH));
        mu_tmp[0] = 0.5;
        mu_tmp[1] = 1.;
        T_tmp = (GAMMA-1.)*mu_tmp[0]*U/P_K;
	if(U > rho*3000./((GAMMA)-1.)*(P_K)/(1.3*(MH))){
        	while(fabs(mu_tmp[0]-mu_tmp[1]) > 1.e-13){
                	mu_tmp[0] = mu_tmp[1];
                	T_tmp = (GAMMA-1.)*mu_tmp[0]*MH*U/P_K;
                	x = 2.*M_PI*P_E*P_K*T_tmp/(P_H*P_H);
                	n_e = rho/(mu_tmp[0]*MH)-(X_H+X_He/4.)*rho/MH;
                	n_HII = pow(x, 1.5)*exp(-CHI_HI/((P_K)*T_tmp))/(n_e+pow(x, 1.5)*exp(-CHI_HI/((P_K)*T_tmp)))*n_H;
                	n_HI = n_H-n_HII;
                	n_HeI = pow(1.+4./n_e*pow(x, 1.5)*exp(-CHI_HeI/((P_K)*T_tmp))+4./(n_e*n_e)*pow(x, 3.)*exp(-(CHI_HeI+CHI_HeII)/((P_K)*T_tmp)), -1.)*n_He;
                	n_HeII = 4.*n_HeI/n_e*pow(x, 1.5)*exp(-CHI_HeI/((P_K)*T_tmp));
                	n_HeIII = n_HeII/n_e*pow(x, 1.5)*exp(-CHI_HeII/((P_K)*T_tmp));
                	mu_tmp[1] = pow(X_H*(1.+n_HII/n_H)+X_He/4.*(1.+n_HeII/n_He+2.*n_HeIII/n_He), -1.);
                	mu_tmp[1] = (mu_tmp[1]+mu_tmp[0])/2.;
                	mu_tmp[1] = (mu_tmp[1]+mu_tmp[0])/2.;
                	mu_tmp[1] = (mu_tmp[1]+mu_tmp[0])/2.;
                	mu_tmp[1] = (mu_tmp[1]+mu_tmp[0])/2.;
        	}
		*mu = mu_tmp[0];
		*T = T_tmp;
	}
	else{
		*mu = pow(X_H*1.+(X_He)/4., -1.);
		*T = (GAMMA-1.)*(*mu)*(MH)*U/(P_K);
	}
	fclose(fp);
}

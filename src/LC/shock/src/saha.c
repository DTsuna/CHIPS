#include "saha.h"
#include "opacity.h"
#include "constant.h"

extern double X, Y;

//n_e=nd[0], n_HII = nd[1], n_HeII = nd[2], n_HeIII = nd[3];
void get_num_density(double rho, double T,  double ndens[])
{
	double mu_tmp[2] = {};
	double x, x_to_3ov2;
	double n_e, n_H, n_He, n_HI, n_HII, n_HeI, n_HeII, n_HeIII;
	double Boltz_HI, Boltz_HeI, Boltz_HeII;

	mu_tmp[0] = 0.5;
	mu_tmp[1] = 1.;
	// precompute Boltzmann factors for speed
	Boltz_HI = exp(-CHI_HI/(P_K*T));
	Boltz_HeI = exp(-CHI_HeI/(P_K*T));
	Boltz_HeII = exp(-CHI_HeII/(P_K*T));
	x = 2.*M_PI*P_E*P_K*T/(P_H*P_H);
	x_to_3ov2 = x * sqrt(x);

	while(fabs(mu_tmp[1]-mu_tmp[0]) > 1.e-15){
		mu_tmp[0] = mu_tmp[1];
		n_H = X*rho/MH;
		n_He = Y/4.*rho/MH;
		n_e = rho/(mu_tmp[0]*MH)-(X+Y/4.)*rho/MH;
		n_HII = x_to_3ov2 * Boltz_HI / (n_e + x_to_3ov2 * Boltz_HI) * n_H;
		n_HI = n_H-n_HII;
		n_HeI = n_He / (1.+4./n_e*x_to_3ov2*Boltz_HeI+4./(n_e*n_e)*x_to_3ov2*x_to_3ov2*Boltz_HeI*Boltz_HeII);
		n_HeII = 4.*n_HeI/n_e* x_to_3ov2 * Boltz_HeI;
		n_HeIII = n_He - n_HeI - n_HeII;
		mu_tmp[1] = 1./ (X*(1.+n_HII/n_H)+Y/4.*(1.+n_HeII/n_He+2.*n_HeIII/n_He));
		mu_tmp[1] = (mu_tmp[1]+mu_tmp[0])/2.;
		mu_tmp[1] = (mu_tmp[1]+mu_tmp[0])/2.;
	}

	ndens[0] = n_e;
	ndens[1] = n_HII;
	ndens[2] = n_HeII;
	ndens[3] = n_HeIII;
}

void saha(double rho, double U, double *mu_outp, double *T_outp)
{
	static double rho_array[30], Udivrho_array[30], mu_array[900], T_array[900];
	static int count = 0;
	static int row, col;

	int i, j;
	double a, b, c;

	if(count == 0){
		set_tables("./LCFiles/mu.txt", rho_array, Udivrho_array, mu_array, &row, &col);
		set_tables("./LCFiles/Tem.txt", rho_array, Udivrho_array, T_array, &row, &col);
		count++;
	}

	if(U/rho < pow(10., Udivrho_array[0])){
		*mu_outp = mu_array[0];
		*T_outp = 2./3.*(*mu_outp)*(MH)/(P_K)*U/rho;
	}
	else if(pow(10., Udivrho_array[0]) <= U/rho && U/rho <= pow(10., Udivrho_array[col-1])){
		for(j = 0; j < col-1; j++){
			if(pow(10., Udivrho_array[j]) <= U/rho && U/rho <= pow(10., Udivrho_array[j+1])){
				break;
			}
		}
		if(rho < pow(10., rho_array[0])){
			*mu_outp = (mu_array[j+1]-mu_array[j])/(Udivrho_array[j+1]-Udivrho_array[j])*(log10(U/rho)-Udivrho_array[j])+mu_array[j];
			*T_outp = 2./3.*(*mu_outp)*(MH)/(P_K)*U/rho;
		}
		else if(pow(10., rho_array[0]) <= rho && rho < pow(10., rho_array[row-1])){
			for(i = 0; i < row-1; i++){
				if(pow(10., rho_array[i]) <= rho && rho <= pow(10., rho_array[i+1])){
					break;
				}
			}
			a = ((mu_array[col*(i+1)+j])-(mu_array[col*i+j]))/(rho_array[i+1]-rho_array[i]);
			b = ((mu_array[col*i+j+1])-(mu_array[col*i+j]))/(Udivrho_array[j+1]-Udivrho_array[j]);
			c = ((mu_array[col*(i+1)+j+1])-(mu_array[col*(i+1)+j])-(mu_array[col*i+j+1])
				+(mu_array[col*i+j]))/((rho_array[i+1]-rho_array[i])*(Udivrho_array[j+1]-Udivrho_array[j]));
			*mu_outp = mu_array[col*i+j]+a*(log10(rho)-rho_array[i])
				+b*(log10(U/rho)-Udivrho_array[j])
				+c*(log10(rho)-rho_array[i])*(log10(U/rho)-Udivrho_array[j]);
			*T_outp = 2./3.*(*mu_outp)*(MH)/(P_K)*U/rho;
		}
		else{
			*mu_outp = (mu_array[(row-1)*col+j+1]-mu_array[(row-1)*col+j])/
					(Udivrho_array[j+1]-Udivrho_array[j])*(log10(U/rho)-Udivrho_array[j])+mu_array[(row-1)*col+j];
			*T_outp = 2./3.*(*mu_outp)*(MH)/(P_K)*U/rho;
		}
	}
	else{
		*mu_outp = mu_array[col-1];
		*T_outp = 2./3.*(*mu_outp)*(MH)/(P_K)*U/rho;
	}
}

	

void set_tables(const char *openfile, double rho[], double U[], double mu[], int *row, int *col)
{
	FILE *fp;
	char filename[256];
	char buf[256], *ascii;
	int i = 0, k = 0, l = 0;

	snprintf(filename, 256, "%s", openfile);
	if((fp = fopen(filename, "r")) == NULL){
		printf("ERROR: Can't open opacity file \"%s\". Check whether it exists.\n", filename);
		exit(EXIT_FAILURE);
	}
	else{
//		printf("Opacity file \"%s\" was set.\n", filename);
	}
	fgets(buf, 256, fp);
	U[0] = atof(strtok(buf, " "));
	i = 1;
	while(1){
		ascii = strtok(NULL, " ");
		if(ascii == NULL){
			break;
		}
		else{
			U[i] = atof(ascii);
			i++;
		}
	}
	*col = i;
	while(fgets(buf, 256, fp) != NULL){
		k = 0;
		rho[l] = atof(strtok(buf, " "));
		while(1){
			ascii = strtok(NULL, " ");
			if(ascii == NULL){
				break;
			}
			else{
				mu[l*(*col)+k] = atof(ascii);
				k++;
			}
		}
		l++;
	}
	*row = l;
	fclose(fp);
}

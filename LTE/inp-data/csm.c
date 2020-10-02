#include <stdio.h>
#include <math.h>

#define M_SUN (1.989E+33)

double rho_csm(double, double);

int main(void)
{
	double Mdot = 1.00e-02*(M_SUN)/365.25/24.0/86400.0000;
	double v_w = 1.00e+07;
	double s = 1.5;
	double r = 1.00e+14, dr = 1.00e+14;
	FILE *fp;

	fp = fopen("csm_profile.txt", "w");
	while(r-dr < 1.00e+16){
		fprintf(fp, "%.3e %.3e\n", r, rho_csm(r, s));
		r += dr;
	}

	fclose(fp);
	return 0;
}

double rho_csm(double r, double s)
{
	return 1.00e-13*pow(r/1.00e+14, -s);
}

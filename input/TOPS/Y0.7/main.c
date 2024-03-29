#include <stdio.h>
#include <string.h>
#include <math.h>

#define eV 1.16045E+04

int main(void)
{
	FILE *fp, *fw, *fw1;
	char filename[][64] = {"5E-04.txt", "1E-03.txt", "3E-03.txt", 
				"1E-02.txt", "3E-02.txt", "1E-01.txt", 
				"3E-01.txt", "1E0.txt", "3E0.txt", "1E+01.txt", "1E+02.txt"
	};
	int i = 0, j = 0, k;
	double dummy[5];
	double temp[] = {5.0e-01*(eV), 1.0*(eV), 3.0*(eV), 
			1.0e+01*(eV), 3.0e+01*(eV), 1.0e+02*(eV),
			3.0e+02*(eV), 1.e+03*(eV), 3.e+03*(eV), 1.e+04*(eV), 1.e+05*(eV)};

	double rho[100], kappa[100], kabs[100];

	fw = fopen("Y0.7.txt", "w");
	fw1 = fopen("Y0.7.abs.txt", "w");

	for(i = 0; i < 11; i++){
		j = 0;
		fp = fopen(filename[i], "r");
		while(fscanf(fp, "%lf %lf %lf %lf %lf", &rho[j], &kappa[j], &kabs[j], &dummy[3], &dummy[4]) != EOF){
			j++;
		}
		if(i == 0){
			for(k = 0; k < j; k++){
				fprintf(fw, "%1.4f ", log10(rho[k]/pow(temp[i], 1.5)));
				fprintf(fw1, "%1.4f ", log10(rho[k]/pow(temp[i], 1.5)));
			}
			fprintf(fw, "\n");
			fprintf(fw1, "\n");
		}
		fprintf(fw, "%1.4f %1.4f", log10(temp[i]), log10(kappa[0]));
		fprintf(fw1, "%1.4f %1.4f", log10(temp[i]), log10(kabs[0]/rho[0]*pow(temp[i], 3.5)));
		for(k = 1; k < j; k++){
			fprintf(fw, " %1.4f", log10(kappa[k]));
			fprintf(fw1, " %1.4f", log10(kabs[k]/rho[k]*pow(temp[i], 3.5)));
		}
		fprintf(fw, "\n");
		fprintf(fw1, "\n");
		fclose(fp);
	}

	fclose(fw);
	fclose(fw1);

	return 0;
}

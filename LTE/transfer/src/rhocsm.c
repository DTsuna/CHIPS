#include <stdio.h>
#include <math.h>
#include "opacity.h"
#include "rhocsm.h"
#include "constant.h"

char csm[256];

double rho_csm(double r)
{
	static double r_c[1000], rho_c[1000];
	static int flag = 0, nsize;

	int i = 0;
	double rho;
	double dammy[5];

	if(flag == 0){
		FILE *fp;
		char filename[256];
		fp = fopen(csm, "r");
		fgets(filename, 512, fp);
		while(fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &dammy[0], &dammy[1], &r_c[i], &dammy[2], &rho_c[i], &dammy[3], &dammy[4]) != EOF){
//			rho_c[i]/=3.15e+07;
			i++;
		}
		nsize = i;
		flag = 1;
		fclose(fp);
	}
	for(i = 0; i < nsize-1; i++){
		if(r < r_c[i+1] && r >= r_c[i]){
			break;
		}
	}
	rho = (log(rho_c[i+1])-log(rho_c[i]))/(log(r_c[i+1])-log(r_c[i]))*(log(r)-log(r_c[i]))+log(rho_c[i]);
	rho = exp(rho);

	return rho;
}

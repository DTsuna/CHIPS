#include <stdio.h>
#include <math.h>
#include "matrix.h"

void inv_3(double *a, double *inv_a, int NSIZE)
{
	double buf;
	int i, j, k;
	
	for(i = 0; i < NSIZE; i++){
		for(j = 0; j < NSIZE; j++){
			inv_a[NSIZE*i+j] =(i == j)?1.:0.;
		}
	}
	for(i = 0; i < NSIZE; i++){
		buf = 1./a[NSIZE*i+i];
		for(j = 0; j < NSIZE; j++){
			a[NSIZE*i+j] *= buf;
			inv_a[NSIZE*i+j] *= buf;
		}
		for(j = 0; j < NSIZE; j++){
			if(i != j){
				buf = a[NSIZE*j+i];
				for(k = 0; k < NSIZE; k++){
					a[NSIZE*j+k] -= a[NSIZE*i+k]*buf;
					inv_a[NSIZE*j+k] -= a[NSIZE*i+k]*buf;
				}
			}
		}
	}
}

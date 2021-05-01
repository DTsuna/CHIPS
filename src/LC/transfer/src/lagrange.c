#include <stdio.h>
#include <math.h>
#include "lagrange.h"

void lagrange_polynomial(double t, int n, double x[], double y[], double *result) //The array size of x, y should be n.
{
	double sum, p;
	int i, j;
	sum = 0.0;

	// 総和
	for(i = 0; i < n; i++){
		p = y[i];
		// 総積
		for(j = 0; j < n; j++){
			if(i != j){
				p *= (t-x[j])/(x[i]-x[j]);
			}
		}
		sum += p;
	}
//	printf("sum = %e\n", sum);
	*result = sum;
}

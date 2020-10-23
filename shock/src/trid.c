#include "trid.h"
#include <math.h>

void trid_matrix_algorithm(double a[], double b[], double c[], double d[], double x[], int n) //a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]
{
	int i;
	double w;

	for(i = 1; i < n; i++){
		w = a[i]/b[i-1];
		b[i] = b[i]-w*c[i-1];
		d[i] = d[i]-w*d[i-1];
	}

	x[n-1] = d[n-1]/b[n-1];
	for(i = n-2; i >= 0; --i){
		x[i] = (d[i]-c[i]*x[i+1])/b[i];
	}
}

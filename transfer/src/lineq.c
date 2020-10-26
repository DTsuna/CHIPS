#include "lineq.h"
#include <math.h>
#include <stdio.h>

/******************************************************************
For more details, please refer to:
http://www-it.sci.waseda.ac.jp/teachers/w405201/CPR2/cprogram18.pdf
These subroutines are extracted from the above URL.
******************************************************************/


void copy_mat(double *src, double *dest, int m, int n)
{
	int i, j;

	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			dest[n*i+j] = src[n*i+j];
		}
	}
}

void copy_mat_t(double *src, double *dest, int m, int n)
{
	int i, j;

	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			dest[i+m*j] = src[n*i+j];
		}
	}
}

void print_mat(char *name, double *src, int m, int n)
{
	int i, j;

	printf("\n%s\n", name);
	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			printf(" %-9.3g", src[n*i+j]);
		}
		printf("\n");
	}
}

void print_vec(char *name, double *src, int n)
{
	int i;

	printf("\n%s\n", name);
	for(i = 0; i < n; i++){
		printf(" %-9.3g", src[i]);
		printf("\n");
	}
}

void calc_ax(double *A, double *x, double *Ax, int n)
{
	int i, j;
	for(i = 0; i < n; i++){
		Ax[i] = 0.;
		for(j = 0; j < n; j++){
			Ax[i] += A[n*i+j]*x[j];
		}
	}
}

void product_mat(double *A, double *B, double *C, int M, int L, int N)
{
	/* return C(M, N) = A(M, L)B(L, N) */
	int i, j, k;

	for(i = 0; i < M; i++){
		for(j = 0; j < N; j++){
			C[N*i+j] = 0;
			for(k = 0; k < L; k++){
				C[N*i+j] += A[L*i+k]*B[N*k+j];
			}
		}
	}
}

double _z(double x)
{
	return fabs(x) < DBL_EPSILON ? 0:x;
}

void print_mat_z(char *name, double *src, int m, int n)
{
	int i, j;
	printf("\n%s\n", name);
	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			printf(" %-9.3g", _z(src[n*i+j]));
		}
		printf("\n");
	}
}

void get_PLU(double *a, int *pvt, int n)
{
	int i, j, k, p;
	double c, r;

	for(k = 0; k < n-1; k++){
		p = k;
		for(i = k; i < n; i++){
			if(fabs(a[p*n+k]) < fabs(a[i*n+k])){
				p = i;
			}
		}
		if(fabs(a[p*n+k]) < DBL_EPSILON){
			printf("failure\n");
			exit(0);
		}
		pvt[k] = p;
		if(p != k){
			for(j = k; j < n; j++){
				c = a[k*n+j];
				a[k*n+j] = a[p*n+j];
				a[p*n+j] = c;
			}
		}
		for(i = k+1; i < n; i++){
			r = a[i*n+k]/a[k*n+k];
			a[i*n+k] = r;
			for(j = k+1; j < n; j++){
				a[i*n+j] -= r*a[k*n+j];
			}
		}
	}
}

void mat_PLU(double *a, double *b, int *pvt, int n)
{
	int i, j, k;
	double c, sum;

	for(k = 0; k < n-1; k++){
		c = b[k];
		b[k] = b[pvt[k]];
		b[pvt[k]] = c;
		for(i = k+1; i < n; i++){
			b[i] -= a[i*n+k]*b[k];
		}
	}

	for(i = n-1; i >= 0; i--){
		sum = 0.;
		for(j = i+1; j < n; j++){
			sum += a[i*n+j]*b[j];
		}
		b[i] = (b[i]-sum)/a[i*n+i];
	}
}

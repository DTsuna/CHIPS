#include "calc_array.h"

void add_array(double A[], double B[], double C[], int n)
{
	int i;
	for(i = 0; i < n; i++){
		C[i] = A[i]+B[i];
	}
}

void subtract_array(double A[], double B[], double C[], int n)
{
	int i;
	for(i = 0; i < n; i++){
		C[i] = A[i]-B[i];
	}
}

void product_array(double A[], double B[], double C[], int n)
{
	int i;
	for(i = 0; i < n; i++){
		C[i] = A[i]*B[i];
	}
}

void devide_array(double A[], double B[], double C[], int n)
{
	int i;
	for(i = 0; i < n; i++){
		C[i] = A[i]/B[i];
	}
}

double min_array(double A[], int n)
{
	double Amin = A[0];
	int i;
	for(i = 1; i < n; i++){
		if(A[i] < Amin){
			Amin = A[i];
		}
	}
	return Amin;
}

double max_array(double A[], int n)
{
	double Amax = A[0];
	int i;
	for(i = 1; i < n; i++){
		if(A[i] > Amax){
			Amax = A[i];
		}
	}
	return Amax;
}


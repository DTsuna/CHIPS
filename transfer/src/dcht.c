#include "dcht.h"

int dcht(double x_in, double *x, int n)
{
	int i;
	int imin_l, imax_l;
	imin_l = 0;
	imax_l = n-1;
	i = (imin_l+imax_l)/2;
	if(x_in < x[0] || x_in > x[n-1]){
	//	printf("error!\n");
		return n;
	}
	else if(x_in >= x[0] && x_in < x[1]){
		return 0;
	}
	else if(x_in >= x[n-2] && x_in <= x[n-1]){
		return n-2;
	}
	else{
		while(1){
			if(x_in < x[i]){
				imax_l = i;
				i = (i+imin_l)/2;
			}
			else if(x_in >= x[i+1]){
				imin_l = i+1;
				i = (i+imax_l+1)/2;
			}
			else{
				break;
			}
		}
		return i;
	}
}

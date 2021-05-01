#include "nr.h"
#include "lineq.h"

void nr_itr(double *x, double *dx, double *J, double *func, int n)
{
	int i, pvt[n];
	for(i = 0; i < n; i++){
		func[i] *= -1.;
	}
	get_PLU(J, pvt, n);
	mat_PLU(J, func, pvt, n);
	for(i = 0; i < n; i++){
		dx[i] = func[i];
	}
}

#include <math.h>
#include "opacity.h"
#include "itgadvu.h"

/*
Variable 'nsize' is the size of array.
*/
void init_U(double U[], double func[], const int nsize)
{
	int i;
	for(i = 0; i < nsize; i++){
		func[i] = U[2*i];
	}
}

void matrix_U(double r[], double U[], double rho[], double rho_ed[], double v_w, double dt, 
	double a[], double b[], double c[], const int nsize)
{
	int i;
	double dr[nsize], dlnrhodr;
	double gamma = 5.0/3.0;

	for(i = 0; i < nsize; i++){
		dr[i] = r[i+1]-r[i];
	}

	dlnrhodr = (log(rho[1])-log(rho_ed[0]))/(2.*dr[0]);
	a[0] = 0.0;
	b[0] = 1.-gamma/(gamma-1.)*v_w*dlnrhodr*dt-v_w/2.*dt/dr[0];
	c[0] = v_w/2.*dt/dr[0];

	for(i = 1; i < nsize-1; i++){
		dlnrhodr = (log(rho[i+1])-log(rho[i-1]))/(2.*dr[i]);
		a[i] = -v_w/2.*dt/dr[i];
		b[i] = 1.-gamma/(gamma-1.)*v_w*dlnrhodr*dt;
		c[i] = v_w/2.*dt/dr[i];
	}

	/*matrix N-1行目*/
	dlnrhodr = (log(rho_ed[1])-log(rho[nsize-2]))/(2.*dr[nsize-1]);
	a[nsize-1] = -v_w/2.*dt/dr[nsize-1];
	b[nsize-1] = 1.-gamma/(gamma-1.)*v_w*dlnrhodr*dt+v_w/2.*dt/dr[nsize-1];
	c[nsize-1] = 0.;
}

void itg_adv_U(double r[], double U[], double rho[], double rho_ed[], double v_w, double dt, const int nsize)
{
	double func[nsize];
	double a[nsize], b[nsize], c[nsize], x[nsize];
	int i;

	init_U(U, func, nsize);
	matrix_U(r, U, rho, rho_ed, v_w, dt, a, b, c, nsize);

	trid_matrix_algorithm(a, b, c, func, x, nsize);
	for(i = 0; i < nsize; i++){
		U[2*i+1] = x[i];
	}
}

#include "W4.h"
#include "opacity.h"
#include "srctrm.h"
#include "saha.h"
#include "constant.h"
#include "nr.h"

#define W4

double func_rhs(double E, double U, double rho)
{
	double mu, T, kappa;
	saha(rho, U, &mu, &T);
	kappa = kappa_p(rho, T);
	return kappa*rho*(P_C)*((P_A)*pow(T, 4.)-E);
}

void diff_eq_src(double E[], double U[], double rho, double dt, double func[])
{
	double f[2];
	f[0] = func_rhs(E[0], U[0], rho);
	f[1] = func_rhs(E[1], U[1], rho);
	func[0] = E[1]-E[0]-dt*f[1];
	func[1] = U[1]-U[0]+dt*f[1];
}

void jacob(double E, double U, double rho, double dt, double J[])
{
	double mu, T, kappa;
	double dfdE, dfdU;
	
	saha(rho, U, &mu, &T);
	kappa = kappa_p(rho, T);
	dfdE = -kappa*rho*(P_C);
	dfdU = 4.*kappa*rho*(P_A)*(P_C)*pow(T, 4.)/U+kappa*rho*(P_A)*(P_C)*pow(T, 4.)/(2.*U); //あってる？？

	J[0] = 1.-dt*dfdE;
	J[1] = -dt*dfdU;
	J[2] = dt*dfdE;
	J[3] = 1.+dt*dfdU;
}

void itg_src(double E[], double U[], double rho, double dt, double tol)
{
	double err;
	double x[2], p[2] = {}, J[4], func[2];
	double dtau = 0.5;
	double mu, T, kappa;
	int count = 0, count_max = 300;

	saha(rho, U[0], &mu, &T);
	kappa = kappa_p(rho, T);
	do{
		diff_eq_src(E, U, rho, dt, func);
		jacob(E[1], U[1], rho, dt, J);
		x[0] = E[1]/E[0]; x[1] = U[1]/U[0];

		J[0] *= E[0]/E[0]; J[1] *= U[0]/E[0];
		J[2] *= E[0]/E[0]; J[3] *= U[0]/E[0];

		func[0] /= E[0];
		func[1] /= E[0];

#ifdef W4
		get_itr_x(x, p, J, func, dtau, 2);
		E[1] = x[0]*E[0]; U[1] = x[1]*U[0];
#else
		nr_itr(x, p, J, func, 2);
		E[1] = (x[0]+p[0])*E[0]; U[1] = (x[1]+p[1])*U[0];
#endif


		err = func[0]*func[0]/E[1]/E[1]+func[1]*func[1]/U[1]/U[1];
		err = sqrt(err/2.);
		count++;
		if(count == count_max){
			break;
		}
		if(isnan(E[1])){
			break;
		}
	}while(err > tol);
}

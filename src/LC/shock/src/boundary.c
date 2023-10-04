#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "boundary.h"
#include "constant.h"
#include "pars.h"
#include "ranhugo.h"
#include "function.h"
#include "W4.h"
#include "nickel.h"

extern double t_exp;
extern pars pdt;
extern double X;

//set physical quantities at upstream of reverse shock.
void set_init_rev(double x, double y_up[], double egn[])
{
	double kappa = 0.2*(1.+X);
	y_up[0] = x/t_exp-egn[0];
	y_up[1] = rho_ej(x, t_exp);
	y_up[2] = 0.0;
	y_up[3] = rad_from_decay_of_nico(pdt.M_ni, pdt.M_ej, pdt.E_ej, kappa, t_exp)/(4.*M_PI*x*x);
}

//set physical quantities at upstream of forward shock.
void set_init_for(double x, double y_up[], double egn[])
{
	y_up[0] = pdt.v_w-egn[1];
	y_up[1] = rho_csm(x);
	y_up[2] = 0.0;
	y_up[3] = egn[2];
}

//guess initial value at downstream of reverse shock.
void set_init_down_rev(double x, double y_up[], double y_down[], double egn[])
{
	double kappa = 0.2*(1.+X);
	y_down[0] = y_up[0]/7.0;
	y_down[1] = y_up[1]*7.0;
	y_down[2] = 6.0/7.0*y_up[1]*y_up[0]*y_up[0];
	y_down[2] = pow(3.0*y_down[2]/(P_A), 0.25);
	y_down[3] = rad_from_decay_of_nico(pdt.M_ni, pdt.M_ej, pdt.E_ej, kappa, t_exp)/(4.*M_PI*x*x);
}

//guess initial value at downstream of forward shock.
void set_init_down_for(double x, double y_up[], double y_down[], double egn[])
{
	y_down[0] = y_up[0]/7.0;
	y_down[1] = y_up[1]*7.0;
	y_down[2] = 6.0/7.0*y_up[1]*y_up[0]*y_up[0];
	y_down[2] = pow(3.0*y_down[2]/(P_A), 0.25);
	y_down[3] = egn[2];
}


/*
The following subroutine is used to calculate boundary conditions at both shocks.
Note that boundary conditions at reverse shock are calculated when the variable *flag* is 0.
If flag = 1, then at forward shock are calculated.
If flag != 0 or 1, then warning is output and ABEND... be careful.
*/
void boundary(double x, double y[], double egn[], int flag, int *info)
{
	int i, j;
	int count = 0, count_max = 100;
	double y_up[4], y_fd[4];
	double func[3], func_fd[3];
	double J[9];
	double a[3], b[3] = {};
	double dtau = 0.5;
	double err, tol = 1.0e-10;

	if(flag == 0){
		set_init_rev(x, y_up, egn);
		set_init_down_rev(x, y_up, y, egn);
		set_init_down_rev(x, y_up, y_fd, egn);
	}
	else if(flag == 1){
		set_init_for(x, y_up, egn);
		set_init_down_for(x, y_up, y, egn);
		set_init_down_for(x, y_up, y_fd, egn);
	}
	else{
		printf("flag error: wrong boundary conditions.\n");
		exit(EXIT_FAILURE);
	}

	func_fd[0] = func_mass(y);
	func_fd[1] = func_momentum(y);
	func_fd[2] = func_energy(y);

	do{
		err = 0.0;
		set_func(y_up, y, func);
		jacobian_rh(y, J);
		for(i = 0; i < 3; i++){
			a[i] = y[i]/y_fd[i];
			func[i] /= func_fd[i];
			for(j = 0;  j < 3; j++){
				J[3*i+j] *= y_fd[j]/func_fd[i];
			}
		}
		get_itr_x(a, b, J, func, dtau, 3);
		for(i = 0; i < 3; i++){
			y[i] = a[i]*y_fd[i];
			err += b[i]*b[i];
		}
		count++;
		if(count == count_max){
			printf("iteration failure. exit. err = %e (boundary conditions)\n", err);
			*info = 1;
			break;
		}
	}while(err > tol);
}

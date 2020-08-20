#ifndef INCLUDED_MAIN_H
#define INCLUDED_MAIN_H

#define P_K (1.38064852E-16)
#define P_C (2.99792458E+10)
#define P_A (7.56591670E-15)
#define P_G (6.67259000E-08)
#define P_E (9.10938356E-28)
#define P_H (6.6260755E-27)
#define MH (1.660539040E-24)
#define SIGMA_TH (0.66524587158E-24)
#define M_SUN (1.989E+33)
#define MU 0.62112*(MH)

#define N 5

typedef struct{
	double n;
	double s;
	double delta;
	double E_ej;
	double M_ej;
	double M_Mloss;
	double v_w;
	double t_ini;
	double D;
	double R_p;
	double alpha;
	double v_t;
	double t_t;
}pars;

void egn_init(double[]);

double t_exp;

#endif

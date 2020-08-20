#include <stdio.h>
#include <omp.h>
#include <time.h>
#include "W4.h"
#include "irk.h"
#include "rk.h"
#include "main.h"
#include "file_rw.h"
#include "opacity.h"
#include "function.h"
#include "solver_neq.h"

#define NSIZE 3

pars pdt; //Parameter set.
double xfd, yrev[N], yfor[N];
//These values are used only to normalize equations and need not be modified while calculating.
//They are set in subroutine boundary_***.

void predict_egn(double, double[], double[], double[], double[], double[], double[], int, int);
void itr_egn(double, double[], FILE*, FILE*, double);

int main(void)
{
	double cpu_time, time1;
	double int_phys[4], ext_phys[4], delta_phys[4];
	double egn[4];
	double y[2*N+1];
	double x_ini;
	double dt;
	double t_elp[NSIZE] = {}, urs_elp[NSIZE] = {}, ufs_elp[NSIZE] = {}, Ffs_elp[NSIZE] = {}, rfs_elp[NSIZE] = {};
	FILE *frev, *ffor, *fegn;
	int i, j, count = 0;

	time1 = omp_get_wtime();
	readin_pars(&pdt); //ファイルを読み込んでパラメータを設定する
	mk_dir_nsd(pdt);
	file_open_pars(&frev, pdt, "rev", "w");
	file_open_pars(&ffor, pdt, "for", "w");
	file_open_pars(&fegn, pdt, "egn", "w");
	t_exp = pdt.t_ini;
	
	egn_init(egn);
	x_ini = r_early(t_exp);
	xfd = x_ini;
	set_yrev(&x_ini, yrev, egn);
	set_yfor(&x_ini, yfor, egn);
	solver_rev(x_ini, int_phys, egn, frev);
	solver_for(int_phys, ext_phys, egn, ffor);
//	printf("%e %e %e\n", ext_phys[0], ext_phys[1], ext_phys[2]);
	while(t_exp < 86400.*100.){
		printf("t = %f d.\n", t_exp/86400.);
		xfd = x_ini;
		itr_egn(x_ini, egn, frev, ffor, time1);
		solver(x_ini, egn, delta_phys, frev, ffor);
		fprintf(fegn, "%f %.15e %.15e %.15e %.15e %.15e\n", t_exp/86400., x_ini, egn[3], egn[0], egn[1], egn[2]);
		printf("%f %.15e %.15e %.15e %.15e %.15e\n", t_exp/86400., x_ini, egn[3], egn[0], egn[1], egn[2]);
		for(i = 0; i < NSIZE-1; i++){
			t_elp[i] = t_elp[i+1];
			urs_elp[i] = urs_elp[i+1];
			ufs_elp[i] = ufs_elp[i+1];
			Ffs_elp[i] = Ffs_elp[i+1];
			rfs_elp[i] = rfs_elp[i+1];
			printf("t_elp = %f, urs_elp = %e, ufs_elp = %e, Ffs_elp = %e, rfs_elp = %e\n", 
				t_elp[i], urs_elp[i], ufs_elp[i], Ffs_elp[i], rfs_elp[i]);
		}
		t_elp[NSIZE-1] = log10(t_exp/86400.);
		urs_elp[NSIZE-1] = log10(egn[0]);
		ufs_elp[NSIZE-1] = log10(egn[1]);
		Ffs_elp[NSIZE-1] = log10(egn[2]);
		rfs_elp[NSIZE-1] = log10(egn[3]);
		printf("t_elp = %f, urs_elp = %e, ufs_elp = %e, Ffs_elp = %e, rfs_elp = %e\n", 
			t_elp[NSIZE-1], urs_elp[NSIZE-1], ufs_elp[NSIZE-1], Ffs_elp[NSIZE-1], rfs_elp[NSIZE-1]);
		dt = 0.2*t_exp;
		t_exp += dt;
		x_ini += egn[0]*dt;
		egn[3] += egn[1]*dt;
		fclose(fegn);
		file_open_pars(&fegn, pdt, "egn", "a");
		count++;
		if(count > NSIZE-1){
			predict_egn(t_exp/86400., egn, t_elp, urs_elp, ufs_elp, Ffs_elp, rfs_elp, NSIZE, 0);
			printf("urs_elp = %e, ufs_elp = %e, Ffs_elp = %e, rfs_elp = %e\n", egn[0], egn[1], egn[2], egn[3]);
		}
	}


	fclose(frev);
	fclose(ffor);
	fclose(fegn);

	cpu_time = omp_get_wtime()-time1;
	printf("CPU time = %f s\n", cpu_time);
	return 0;
}

void predict_egn(double t, double egn[], double t_elp[], double urs_elp[], double ufs_elp[], double Ffs_elp[], double rfs_elp[], int n, int m)
{
	int i;
	lagrange_polynomial(log10(t), n-m, t_elp+m, urs_elp+m, &egn[0]);
	lagrange_polynomial(log10(t), n-m, t_elp+m, ufs_elp+m, &egn[1]);
	lagrange_polynomial(log10(t), n-m, t_elp+m, Ffs_elp+m, &egn[2]);
	lagrange_polynomial(log10(t), n-m, t_elp+m, rfs_elp+m, &egn[3]);
	for(i = 0; i < 4; i++){
		egn[i] = pow(10., egn[i]);
	}
}

void egn_init(double egn[])
{
	egn[0] = v_early(t_exp);
	egn[1] = 1.01*egn[0];
	egn[2] = 2.e+14;
	egn[3] = 1.1*r_early(t_exp);
//	egn[0] = 7.142638813761748e+08;
//	egn[1] = 7.146234110234046e+08;
//	egn[2] = 4.732713276491505e+13;
//	egn[3] = 1.340088805774638e+14;
//	2.000000 1.331438235914642e+14 1.340088805774638e+14 7.142638813761748e+08 7.146234110234046e+08 4.732713276491505e+13
}

void itr_egn(double x_ini, double egn[], FILE *fprev, FILE *fpfor, double time1)
{
	int i, j;
	double err;
	double x = x_ini;
	double delta_phys[4], dphysfd[4], d_dphys[16];
	double ddpdegn[16], inv[16];
	double egnfd[4], degn[4];
	double dtau = 0.2;
	double p[4] = {0., 0., 0., 0.};
	double legn[16];
	double tol = 0.03;
	double dammy[4];

	set_yrev(&x, yrev, egn);
	set_yfor(&x, yfor, egn);
	solver(x_ini, egn, delta_phys, fprev, fpfor);
	xfd = x_ini;
	for(i = 0; i < 4; i++){
		dphysfd[i] = fabs(delta_phys[i]);
		egnfd[i] = fabs(egn[i]);
		degn[i] = 1.e-05*egnfd[i];
	}
	degn[0] *= 0.01;
	degn[1] *= 0.01;
	degn[3] *= 0.01;
	do{
		err = 0.;
		printf("x_ini/t-u_rs = %e\n", x_ini/t_exp-egn[0]);
		#pragma omp parallel for
		for(i = 0; i < 5; i++){
			if(i != 4){
				memcpy(legn+4*i, egn, sizeof(double)*4);
//				printf("%e %e %e %e\n", legn[4*i], legn[4*i+1], legn[4*i+2], legn[4*i+3]);
				legn[4*i+i] += degn[i];
				solver(x_ini, legn+4*i, d_dphys+4*i, NULL, NULL);
//				printf("%e %e %e %e\n", d_dphys[4*i], d_dphys[4*i+1], d_dphys[4*i+2], d_dphys[4*i+3]);
				legn[4*i+i] -= degn[i];
				printf("i = %d, dv = %e,  dp = %e, dF = %e, dr = %e\n", i, d_dphys[4*i], d_dphys[4*i+1], d_dphys[4*i+2], d_dphys[4*i+3]);
			}
			else{
				solver(x_ini, egn, delta_phys, NULL, NULL);
				printf("i = %d, dv = %e,  dp = %e, dF = %e, dr = %e\n", 
				i, delta_phys[0], delta_phys[1], delta_phys[2], delta_phys[3]);
			}
		}
				solver(x_ini, egn, dammy, fprev, fpfor);
		

		for(i = 0; i < 4; i++){
			for(j = 0; j < 4; j++){
				ddpdegn[4*j+i] = (d_dphys[4*i+j]-delta_phys[j])/degn[i];
				ddpdegn[4*j+i] *= egnfd[i]/dphysfd[j];
			}
		}
		for(i = 0; i < 4; i++){
			for(j = 0; j < 4; j++){
				printf("%e ", ddpdegn[4*i+j]);
			}
			printf("\n");
		}
		for(i = 0; i < 4; i++){
			egn[i] /= egnfd[i];
			delta_phys[i] /= dphysfd[i];
		}
		get_itr_x_tol(egn, p, ddpdegn, delta_phys, dtau, tol, 4);
//		get_itr_x(egn, p, ddpdegn, delta_phys, dtau, 4);
//		inv_3(ddpdegn, inv, 4);
//		for(i = 0; i < 4; i++){
//			for(j = 0; j < 4; j++){
//				printf("%e ", inv[4*i+j]);
//			}
//			printf("\n");
//		}
		for(i = 0; i < 4; i++){
			egn[i] *= egnfd[i];
		}
		delta_phys[0] *= dphysfd[0]/egn[0];
		delta_phys[1] *= dphysfd[1]/(rho_csm(egn[3])*egn[1]*egn[1]);
		delta_phys[2] *= dphysfd[2]/egn[2];
		delta_phys[3] *= dphysfd[3]/(egn[3]-x_ini);
		for(i = 0; i < 4; i++){
			err += delta_phys[i]*delta_phys[i];
//			delta_phys[i] *= dphysfd[i];
		}
		delta_phys[0] *= egn[0];
		delta_phys[1] *= (rho_csm(egn[3])*egn[1]*egn[1]);
		delta_phys[2] *= egn[2];
		delta_phys[3] *= (egn[3]-x_ini);
		err = sqrt(err/4.0);
		printf("u_rs = %.15e, u_fs = %.15e, F_fs = %.15e, r_fs = %.15e, dv = %e, dp = %e, dF = %e, dr = %e\n", 
			egn[0], egn[1], egn[2], egn[3], delta_phys[0], delta_phys[1], delta_phys[2], delta_phys[3]);
		printf("%e, %e, %e, %e\n", p[0], p[1], p[2], p[3]);
		printf("err = %e, elapsed time = %f min.\n", err, (omp_get_wtime()-time1)/60.);
	}while(err > 1.e-03);
}

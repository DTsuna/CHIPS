#include "file_rw.h"

void readin_pars(pars *pdt)
{
	FILE *fp;
	fp = fopen("./pars.txt", "r");
	if(fp == NULL){
		printf("FATAL ERROR: Can't read file.\n");
		exit(EXIT_FAILURE);
	}
	if(fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		&pdt->t_ini, &pdt->E_ej, &pdt->M_ej, &pdt->M_Mloss, &pdt->n, &pdt->s, &pdt->delta, &pdt->v_w, &pdt->alpha) == EOF){
		printf("no pars? Look into your parameter file and check parameters exist.\n");
		exit(EXIT_FAILURE);
	}
	else{
		pdt->M_ej *= M_SUN;
		pdt->D = (3.-pdt->s)/1.99*pow(10., 16.*pdt->s-14.)*pdt->M_Mloss*(1.e+07/pdt->v_w);
		pdt->t_ini *= 86400.;
		pdt->R_p = 2.*(P_G)*pdt->M_ej/(pdt->v_w*pdt->v_w);
		pdt->v_t = sqrt(2.*(5.-pdt->delta)*(pdt->n-5.)/((3.-pdt->delta)*(pdt->n-3.))*pdt->E_ej/pdt->M_ej);
	}
	printf("Set parameters:\n");
	printf("E_ej: %e erg/s, M_ej: %f M_sun, Mdot: %f M_sun/yr, n: %f, s: %f, delta: %f, v_w: %e cm/s, alpha: %f\n", 
		pdt->E_ej, pdt->M_ej/M_SUN, pdt->M_Mloss, pdt->n, pdt->s, pdt->delta, pdt->v_w, pdt->alpha);
	fclose(fp);
}

void mk_dir_nsd(pars pdt) //特定のパラメータに対してディレクトリを作成する
{
	char mk_dir[256];
	snprintf(mk_dir, 256, "./outp-data/n%ds%dd%d", (int)(round(pdt.n*100.)), (int)(round(pdt.s*100.)), (int)(round(pdt.delta*100.)));
	if(mkdir(mk_dir, 0700) == 0){
		printf("succeeded to make directory \"%s\"\n", mk_dir);
	}
	else{
		printf("directory \"%s\" already exists.\n", mk_dir);
	}
}

void file_open_pars(FILE **fp, pars pdt, const char *region, const char *mode) //特定のパラメータに対してアウトプットファイルを作成する
{
	char filename[256], mk_dir[256];
	
	snprintf(mk_dir, 256, "./outp-data/n%ds%dd%d", (int)(round(pdt.n*100.)), (int)(round(pdt.s*100.)), (int)(round(pdt.delta*100.)));
	snprintf(filename, 256, "%s/%s_%d_%d_%d_%d.txt", mk_dir, region, 
		(int)(round(pdt.E_ej*1.e-48)), (int)(round(pdt.M_ej/M_SUN*100.)), (int)(round(pdt.M_Mloss*1.e+07)), (int)(round(pdt.v_w*1.e-06)));

	if((*fp = fopen(filename, mode)) == NULL){
		printf("\x1b[33m");
		printf("ERROR: Can't open output file.[m\n");
		exit(EXIT_FAILURE);
	}
	else{
		printf("Output file \"%s\" was created.\n", filename);
	}
}

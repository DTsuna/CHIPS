#include <stdio.h>
#include <math.h>

#define CELL_NUM 10000
#define ELEM_NUM 19
#define CUT 1         	

#define MSUN 1.9884e33
#define PI 3.1415926535
#define G 6.6743e-8


int main(){
	printf("Program Start\n");
	FILE *file_in;
	file_in = fopen("hoge2.txt", "r");

	int ret ;
	int i = 0;
	int zone;
	double m_msun[CELL_NUM + 5], m_g[CELL_NUM + 5], dmass_g[CELL_NUM + 5], radius[CELL_NUM + 5], density[CELL_NUM + 5], pressure[CELL_NUM + 5], temperature[CELL_NUM + 5];
	double x[ELEM_NUM][CELL_NUM + 5];
	double temp_m_msun, temp_m_g, temp_dmass_g, temp_radius, temp_density, temp_pressure, temp_temperature;
	double temp_x[ELEM_NUM];

	while( (ret = fscanf( file_in, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &zone, &temp_m_msun, &temp_m_g, &temp_dmass_g, &temp_radius, &temp_density, &temp_pressure, &temp_temperature, &temp_x[0], &temp_x[1], &temp_x[2], &temp_x[3], &temp_x[4], &temp_x[5], &temp_x[6], &temp_x[7], &temp_x[8], &temp_x[9], &temp_x[10], &temp_x[11], &temp_x[12], &temp_x[13], &temp_x[14], &temp_x[15], &temp_x[16], &temp_x[17], &temp_x[18])) != EOF && i < 40000){
		m_msun[i] = temp_m_msun;
		m_g[i] = temp_m_g;
		dmass_g[i] = temp_dmass_g;
		radius[i] = temp_radius;
		density[i] = temp_density;
		pressure[i] = temp_pressure;
		temperature[i] = temp_temperature;
		for (int j = 0; j < ELEM_NUM; ++j)
		{
			x[j][i] = temp_x[j];
		}
		i++;
	}
	int data_size = i ;
	fclose(file_in);



	for (int i = 1; i < CELL_NUM +1; ++i)
	{
		pressure[i] = sqrt(pressure[i-1]*pressure[i]);
	}

	double grav;
	double grav_sum_pressure[CELL_NUM +1];

	grav_sum_pressure[CELL_NUM] = pressure[CELL_NUM];
	grav_sum_pressure[0] = pressure[0];


	for (int i = CELL_NUM-1; i > 0; --i)
	{
		grav = (G*(m_g[i] + dmass_g[i]*0.5)*density[i+1])/(pow((0.5*radius[i] + 0.5*radius[i+1]),2));
		grav_sum_pressure[i] = grav_sum_pressure[i+1] + grav*(radius[i+1] - radius[i]);
	}

/*
	FILE*file_out2;
	file_out2 = fopen("grav_integrate_mod.txt", "w");
	for (int i = 0; i < CELL_NUM +1; ++i)
	{
		fprintf(file_out2, "%d %e %e %e\n", i , pressure[i], grav_sum_pressure[i], (pressure[i] - grav_sum_pressure[i])/pressure[i]);
	}
*/


	FILE *file_out;
	file_out = fopen("hoge3.txt", "w");
	fprintf(file_out, "%d\n", CELL_NUM +1);
	fprintf(file_out, "grid  m_r(Msun)     m_r(g)     dmass(g)     radius(cm)  density(g/cm^3) pressure(erg/cm^3) temperature(K) h1          he3           he4           c12           n14           o16           ne20           mg24          si28          s32          ar36          ca40          ti44          cr48          cr56          fe52         fe54         fe56         ni56\n" );

	int k;
	for(k = 0; k < CELL_NUM +1; k++){
		fprintf(file_out,"%d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", k+1, m_g[k]/MSUN, m_g[k], dmass_g[k], radius[k], density[k], grav_sum_pressure[k], temperature[k], x[0][k], x[1][k], x[2][k], x[3][k], x[4][k], x[5][k], x[6][k], x[7][k], x[8][k], x[9][k], x[10][k], x[11][k], x[12][k], x[13][k], x[14][k], x[15][k], x[16][k], x[17][k], x[18][k]);
	}
	fclose(file_out);

	return 0;
}

#include <stdio.h>
#include <math.h>

#define CELL_NUM 10000 //hydroで用いるセル数指定
#define ELEM_NUM 19
#define INPUT_MAX 20000

#define MASS_CUT 4.0 //in solar mass unit

#define MSUN 1.9884e33
#define PI 3.1415926535
#define G 6.6743e-8

double max(double n1,double n2){
	if (n1 > n2)
	{
		return n1;
	}
	else
	{
		return n2;
	}
}

int main(){
	printf("Program Start\n");
	FILE *file_in;
	file_in = fopen("hoge1.txt", "r");

	int ret ;
	int i = 0;
	int zone;
	double m_msun[CELL_NUM + 5], m_g[CELL_NUM + 5], dmass_g[CELL_NUM + 5], radius[CELL_NUM + 5], density[CELL_NUM + 5], pressure[CELL_NUM + 5], temperature[CELL_NUM + 5];
	double x[ELEM_NUM][CELL_NUM + 5];

	double temp_m_msun, temp_m_g, temp_dmass_g, temp_radius, temp_density, temp_pressure, temp_temperature;
	double temp_x[ELEM_NUM];


	while( (ret = fscanf( file_in, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &zone, &temp_m_msun, &temp_m_g, &temp_dmass_g, &temp_radius, &temp_density, &temp_pressure, &temp_temperature, &temp_x[0], &temp_x[1], &temp_x[2], &temp_x[3], &temp_x[4], &temp_x[5], &temp_x[6], &temp_x[7], &temp_x[8], &temp_x[9], &temp_x[10], &temp_x[11], &temp_x[12], &temp_x[13], &temp_x[14], &temp_x[15], &temp_x[16], &temp_x[17], &temp_x[18])) != EOF && i < INPUT_MAX){
		if(temp_m_msun > MASS_CUT)
		{
			m_msun[i] = temp_m_msun;
			m_g[i] = temp_m_g;
			dmass_g[i] = temp_dmass_g;
			radius[i] = temp_radius;
			density[i] = temp_density;
			pressure[i] = temp_pressure;
			temperature[i] = temp_temperature;

			double cut = 1e-40;
			for (int j = 0; j < ELEM_NUM; ++j)
			{
				x[j][i] = max(temp_x[j],cut);
			}
			i++;
		}
	}
	int data_size = i ;
	printf("mesa data i = %d\n",i);
	fclose(file_in);


	double dmass = (m_g[data_size - 1] - m_g[0])/(CELL_NUM);
	double new_m_g[CELL_NUM + 1];

	/*dmassを一定にする場合
	for (int i = 0; i < CELL_NUM + 1; ++i)
	{
		new_m_g[i] = m_g[0] + i *dmass;
	}
	*/


	//dmassに傾斜をつける場合
	double new_dmass[CELL_NUM + 1];
	new_dmass[0] = m_g[0];
	for (int m = 1; m < CELL_NUM + 1; ++m)
	{
		new_dmass[m] = ( (-1.6/(CELL_NUM-1) )* (m-1) + 1.8)*dmass;
	}	
	for (int i = 0; i < CELL_NUM + 1; ++i)
	{
		new_m_g[i] = 0;
		for (int ii = 0; ii < i +1; ++ii)
		{
			new_m_g[i] = new_m_g[i] + new_dmass[ii];
		}
	}
	//ここまで


	printf("old M = %e \n", m_g[data_size -1]);
	printf("new M = %e \n", new_m_g[CELL_NUM ]);	

	
	double new_radius[CELL_NUM + 1];

	double new_pressure[CELL_NUM + 1];
	double new_temperature[CELL_NUM + 1];

	double new_x[ELEM_NUM][CELL_NUM+1];

	new_radius[0] = radius[0];
	new_radius[CELL_NUM] = radius[data_size - 1];

	new_pressure[0] = pressure[0];
	new_pressure[CELL_NUM] = pressure[data_size - 1];
	new_temperature[0] = temperature[0];
	new_temperature[CELL_NUM] = temperature[data_size - 1];

	for (int j = 0; j < ELEM_NUM; ++j)
	{
		new_x[j][0] = x[j][0];
		new_x[j][CELL_NUM] = x[j][data_size - 1];
	}

	double log_new_rad, log_new_pre, log_new_temp;
	double log_new_x[ELEM_NUM];

	for (int i = 1; i < CELL_NUM; ++i)
	{
		int j = 0;
		while (new_m_g[i] > m_g[j+1]){
			j = j+1;
		}


		log_new_rad =  log10(radius[j]) + ((log10(radius[j+1]) - log10(radius[j]))*(new_m_g[i] - m_g[j]))/(m_g[j+1] - m_g[j]);

		log_new_pre =  log10(pressure[j]) + ((log10(pressure[j+1]) - log10(pressure[j]))*(new_m_g[i] - m_g[j]))/(m_g[j+1] - m_g[j]);
		log_new_temp =  log10(temperature[j]) + ((log10(temperature[j+1]) - log10(temperature[j]))*(new_m_g[i] - m_g[j]))/(m_g[j+1] - m_g[j]);

		for (int k = 0; k < ELEM_NUM; ++k)
		{
			log_new_x[k] = log10(x[k][j]) + ((log10(x[k][j+1]) - log10(x[k][j]))*(new_m_g[i] - m_g[j]))/(m_g[j+1] - m_g[j]);
		}

		new_radius[i] = pow(10,log_new_rad);
		new_pressure[i] = pow(10,log_new_pre);
		new_temperature[i] = pow(10,log_new_temp);

		for (int l = 0; l < ELEM_NUM; ++l)
		{
			new_x[l][i] = pow(10,log_new_x[l]);
		}

	}

	double new_density[CELL_NUM + 1];
	new_density[0] = density[0];
	for (int i = 1; i <CELL_NUM +1; ++i)
	{
		new_density[i] = (3*(new_m_g[i] - new_m_g[i-1]))/(4*PI*(pow(new_radius[i],3) - pow(new_radius[i-1],3)));
	}


	FILE *file_out;
	file_out = fopen("hoge2.txt", "w");

	int k;
	for(k = 0; k < CELL_NUM +1; k++){
		fprintf(file_out,"%d %e %e %e %e %e %e %e ", k+1, new_m_g[k]/MSUN, new_m_g[k], new_dmass[k], new_radius[k], new_density[k], new_pressure[k], new_temperature[k]);
		fprintf(file_out,"%e %e %e %e %e %e %e %e %e %e ", new_x[0][k], new_x[1][k], new_x[2][k], new_x[3][k], new_x[4][k], new_x[5][k], new_x[6][k], new_x[7][k], new_x[8][k], new_x[9][k]);
		fprintf(file_out,"%e %e %e %e %e %e %e %e %e\n", new_x[10][k], new_x[11][k], new_x[12][k], new_x[13][k], new_x[14][k], new_x[15][k], new_x[16][k], new_x[17][k], new_x[18][k]);
	}
	fclose(file_out);

	return 0;
}
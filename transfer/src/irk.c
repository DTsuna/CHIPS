#include "irk.h"
#include "calc_array.h"

void func_setup(double g[], double y[], int rank, int dim, double h, double x, double k[], double a[], double c[], void (*dervs)(double, double[], double[]))
{
	int p, q, r;
	double y0[rank*dim];
	for(q = 0; q < rank; q++){
		memcpy(y0+q*dim, y, sizeof(double)*dim);
	}

//	printf("1\n");
	for(p = 0; p < dim; p++){
		for(q = 0; q < rank; q++){
			for(r = 0; r < rank; r++){
				y0[q*dim+p] += h*a[q*rank+r]*k[r*dim+p];
			}
		}
	}
	
//	printf("2\n");
	for(q = 0; q < rank; q++){
		(*dervs)(x+c[q]*h, y0+q*dim, g+q*dim);
		for(p = 0; p < dim; p++){
			g[q*dim+p] = k[q*dim+p]-g[q*dim+p];
		}
	}
//	printf("3\n");
}

void irk24(double y[], int m, double h, double x, double yout[], double tol, double *errl, void (*dervs)(double, double[], double[]))
{
	a_24[0] =   0.2500000000000000;
	a_24[1] = -0.03867513459481286;
	a_24[2] =   0.5386751345948129;
	a_24[3] =   0.2500000000000000;
	b_24[0] =   0.5000000000000000;
	b_24[1] =   0.5000000000000000;
	c_24[0] =   0.2113248654051871;
	c_24[1] =   0.7886751345948129;

	double k[N_24*m], k1[N_24*m], dk[N_24*m], delta_k[N_24*m];
	double g[N_24*m], g0[N_24*m], g1[N_24*m], dgdk[(N_24)*m*(N_24)*m], inv_dgdk[(N_24)*m*(N_24)*m];
	double dydx[m];
	double err = 0.;
	double mink;
	double kk[N_24*m];
	int count_max = 300;
	int p, q, r, p0, q0, i, j, count = 0;
	(*dervs)(x, y, dydx);
	for(p = 0; p < N_24; p++){
		memcpy(k+p*m, dydx, sizeof(double)*m);
		memcpy(dk+p*m, dydx, sizeof(double)*m);
	}
	for(p = 0; p < m; p++){
		for(q = 0; q < N_24; q++){
			dk[q*m+p] = fabs(k[q*m+p])*1.e-04;
			k[q*m+p] = 0.;
			//dk[q*m+p] = 1.e-04;
			k1[q*m+p] = 0.;
		}
	}

	do{
		err = 0.;
		func_setup(g, y, N_24, m, h, x, k, a_24, c_24, dervs);
		for(p = 0; p < m; p++){
			for(q = 0; q < N_24; q++){
				k[q*m+p] += dk[q*m+p];
				func_setup(g0, y, N_24, m, h, x, k, a_24, c_24, dervs);
				k[q*m+p] -= 2.*dk[q*m+p];
				func_setup(g1, y, N_24, m, h, x, k, a_24, c_24, dervs);
				k[q*m+p] += dk[q*m+p];
				for(p0 = 0; p0 < m; p0++){
					for(q0 = 0; q0 < N_24; q0++){
						dgdk[(q0*m+p0)*((N_24)*m)+q*m+p] = (g0[q0*m+p0]-g1[q0*m+p0])/(2.*dk[q*m+p]);
					}
				}
			}
		}
		get_itr_x(k, k1, dgdk, g, 0.15, N_24*m);
	//	get_x_LH(k, k1, dgdk, g, 0.2, N_24*m);
		for(i = 0; i < N_24*m; i++){
			err += g[i]*g[i];
		}
		count++;
		err = sqrt(err/(double)(m*N_24));

		for(i = 0; i < m*N_24; i++){
			kk[i] = fabs(k[i]);
		}
		mink = max_array(k, m*N_24);
	}while(err > tol && count < count_max);

	for(i = 0; i < m; i++){
		yout[i] = y[i];
		for(j = 0; j < N_24; j++){
			yout[i] += h*b_24[j]*k[j*m+i];
		}
	}
	if(count == count_max){
	//	printf("iteration failure. err = %e\n", err);
	}
	*errl = err;
//	printf("count = %d\n", count);
}


void irk36(double y[], int m, double h, double x, double yout[], double tol, double *err, void (*dervs)(double, double[], double[]))
{
	a_36[0]     =     0.1388888888888889;
	a_36[1]     =   -0.03597666752493889;
	a_36[2]     =   0.009789444015308318;
	a_36[3]     =     0.3002631949808646;
	a_36[4]     =     0.2222222222222222;
	a_36[5]     =   -0.02248541720308681;
	a_36[6]     =     0.2679883337624694;
	a_36[7]     =     0.4804211119693833;
	a_36[8]     =     0.1388888888888889;
	b_36[0]     =     0.2777777777777778;
	b_36[1]     =     0.4444444444444444;
	b_36[2]     =     0.2777777777777778;
	c_36[0]     =     0.1127016653792583;
	c_36[1]     =     0.5000000000000000;
	c_36[2]     =     0.8872983346207417;


/*	a_36[0] =  0.16666666666666666666666666;
	a_36[1] = -0.33333333333333333333333333;
	a_36[2] =  0.16666666666666666666666666;
	a_36[3] =  0.16666666666666666666666666;
	a_36[4] =  0.41666666666666666666666666;
	a_36[5] =  0.08333333333333333333333333;
	a_36[6] =  0.16666666666666666666666666;
	a_36[7] =  0.66666666666666666666666666;
	a_36[8] =  0.16666666666666666666666666;
	b_36[0] =  0.16666666666666666666666666;
	b_36[1] =  0.66666666666666666666666666;
	b_36[2] =  0.16666666666666666666666666;
	c_36[0] =  0.0;
	c_36[1] =  0.50000000000000000000000000;
	c_36[2] =  1.00000000000000000000000000;
*/
	int count_max = 200;
	double k[N_36*m], k1[N_36*m], dk[N_36*m], delta_k[N_36*m];
	double g[N_36*m], g0[N_36*m], g1[N_36*m], dgdk[(N_36)*m*(N_36)*m], inv_dgdk[(N_36)*m*(N_36)*m];
	int pvt[N_36*m];
	double dydx[m];
	int p, q, r, p0, q0, i, j, count = 0;
	(*dervs)(x, y, dydx);
	for(p = 0; p < N_36; p++){
		memcpy(k+p*m, dydx, sizeof(double)*m);
		memcpy(dk+p*m, dydx, sizeof(double)*m);
	}
	for(p = 0; p < m; p++){
		for(q = 0; q < N_36; q++){
				dk[q*m+p] = fabs(k[q*m+p])*1.e-04;
		}
	}

	do{
		*err = 0.;
		func_setup(g, y, N_36, m, h, x, k, a_36, c_36, dervs);
		for(p = 0; p < m; p++){
			for(q = 0; q < N_36; q++){
				k[q*m+p] += dk[q*m+p];
				func_setup(g0, y, N_36, m, h, x, k, a_36, c_36, dervs);
				k[q*m+p] -= 2.*dk[q*m+p];
				func_setup(g1, y, N_36, m, h, x, k, a_36, c_36, dervs);
				k[q*m+p] += dk[q*m+p];
				for(p0 = 0; p0 < m; p0++){
					for(q0 = 0; q0 < N_36; q0++){
						dgdk[(q0*m+p0)*((N_36)*m)+q*m+p] = (g0[q0*m+p0]-g1[q0*m+p0])/(2.*dk[q*m+p]);
					}
				}
			}
		}
		get_itr_x(k, k1, dgdk, g, 0.5, N_36*m);
		for(i = 0; i < N_36*m; i++){
			*err += g[i]*g[i];
		}
		/*get_PLU(dgdk, pvt, N_36*m);
		mat_PLU(dgdk, g, pvt, N_36*m);
		for(i = 0; i < N_36*m; i++){
			k[i] += g[i];
		}*/
		count++;
		*err = sqrt(*err/(double)(m*N_36));
	}while(*err > tol && count < count_max);
	for(i = 0; i < m; i++){
		yout[i] = y[i];
		for(j = 0; j < N_36; j++){
			yout[i] += h*b_36[j]*k[j*m+i];
		}
	}
	if(count == count_max){
		printf("iteration failure.\n");
		printf("err = %e\n", *err);
	/*	for(i = 0; i < m; i++){
			yout[i] = 1.e+300;
		}*/
	}
}

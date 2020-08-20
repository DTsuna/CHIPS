#include "W4.h"

void get_UL(double *J, double *U, double*L, int n)
{
	int i, k, m;
	double sum0, sum1;
	for(i = 0; i < n*n; i++){
		U[i] = 0.;
		L[i] = 0.;
	}
	for(i = n-1; i >= 0; i--){
		U[n*i+n-1] = J[n*i+n-1];
		L[n*(n-1)+i] = J[n*(n-1)+i]/U[n*(n-1)+n-1];
		if(i == n-1){
			L[n*(n-1)+i] = 1.;
		}
	}
	for(k = n-2; k > 0; k--){
		for(i = k; i >= 0; i--){
			sum0 = 0.;
			for(m = k+1; m < n; m++){
				sum0 += U[n*i+m]*L[n*m+k];
			}
			U[n*i+k] = J[n*i+k]-sum0;
			sum1 = 0.;
			for(m = k+1; m < n; m++){
				sum1 += U[n*k+m]*L[n*m+i];
			}
			L[n*k+i] = (J[n*k+i]-sum1)/U[n*k+k];
			if(i == k){
				L[n*k+i] = 1.;
			}
		}
	}
	U[0] = 0.;
	for(i = 1; i < n; i++){
		U[0] += U[i]*L[n*i];
	}
	U[0] = J[0]-U[0];
	L[0] = 1.;
}

void get_U_inv(double *U, double *U_inv, int n)
{
	int i, j, k;
	double sum = 0.;
	
	for(i = 0; i < n*n; i++){
		U_inv[i] = 0.;
	}

	for(i = 0; i < n; i++){
		U_inv[n*i+i] = 1./U[n*i+i];
	}
	
	for(j = n-1;  j > 0; j--){
		for(i = j-1; i >= 0; i--){
			sum = 0.;
			for(k = i+1; k <= j; k++){
				sum += U[n*i+k]*U_inv[n*k+j];
			}
			U_inv[n*i+j] = -sum/U[n*i+i];
		}
	}
}

void get_L_inv(double *L, double *L_inv, int n)
{
	int i, j, k;
	double sum = 0.;

	for(i = 0; i < n*n; i++){
		L_inv[i] = 0.;
	}

	for(i = 0; i < n; i++){
		L_inv[n*i+i] = 1.;
	}

	for(j = 0; j < n-1; j++){
		for(i = j+1; i < n; i++){
			sum = 0.;
			for(k = j; k < i; k++){
				sum += L[n*i+k]*L_inv[n*k+j];
			}
			L_inv[n*i+j] = -sum;
		}
	}
}

void get_x_p(double *x_old, double *p_old, double *x, double *p, double *U_inv, double *L_inv, double *func, double dtau, int n)
{
	int i, j;
	double dx[n], dp[n];
	for(i = 0; i < n; i++){
		dx[i] = 0.;
		dp[i] = 0.;
		for(j = 0; j < n; j++){
			//dx[i] += 0.5*L_inv[n*i+j]*p_old[j];
			dx[i] += dtau*L_inv[n*i+j]*p_old[j];
			//dp[i] += -0.5*U_inv[n*i+j]*func[j];
			dp[i] += -dtau*U_inv[n*i+j]*func[j];
		}
		x[i] = x_old[i]+dx[i];
		p[i] = p_old[i]+dp[i]-2.*dtau*p_old[i];
		//p[i] = dp[i];
	}
	memcpy(x_old, x, sizeof(double)*n);
	memcpy(p_old, p, sizeof(double)*n);
}

void get_dx_dp(double *x_old, double *p_old, double *dx, double *dp, double *U_inv, double *L_inv, double *func, double dtau, int n)
{
	int i, j;
	for(i = 0; i < n; i++){
		dx[i] = 0.;
		dp[i] = 0.;
		for(j = 0; j < n; j++){
			dx[i] += dtau*L_inv[n*i+j]*p_old[j];
			dp[i] += -dtau*U_inv[n*i+j]*func[j];
		}
	}
}


void get_itr_x(double *x, double *p, double *J, double *func, double dtau, int n)
{
	double U[n*n], L[n*n], U_inv[n*n], L_inv[n*n];
	double x_new[n], p_new[n];
	get_UL(J, U, L, n);
	get_U_inv(U, U_inv, n);
	get_L_inv(L, L_inv, n);
	get_x_p(x, p, x_new, p_new, U_inv, L_inv, func, dtau, n);
	memcpy(x, x_new, sizeof(double)*n);
	memcpy(p, p_new, sizeof(double)*n);
}

void get_itr_x_tol(double *x, double *p, double *J, double *func, double dtau, double tol, int n)
{
	double U[n*n], L[n*n], U_inv[n*n], L_inv[n*n];
	double x_new[n], p_new[n];
	double dx[n], dp[n];
	int i;

	get_UL(J, U, L, n);
	get_U_inv(U, U_inv, n);
	get_L_inv(L, L_inv, n);
	get_dx_dp(x, p, dx, dp, U_inv, L_inv, func, dtau, n);
	
	for(i = 0; i < n; i++){
		if(fabs(dx[i]/x[i]) > tol){
			dx[i] = tol*dx[i]/fabs(dx[i])/x[i];
		}
		x[i] += dx[i];
		p[i] = p[i]+dp[i]-2.*dtau*p[i];
	}
}

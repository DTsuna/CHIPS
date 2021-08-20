#include "opacity.h"
#include "constant.h"

#define SWITCH_OP_FREE_FREE 0
#define WARN 0

double X, Y;

void set_opacity(const char *openfile, opacity *op)
{
	FILE *fp;
	char filename[256];
	char buf[256], *ascii;
	int i = 0, k = 0, l = 0;

	snprintf(filename, 256, "%s", openfile);
	if((fp = fopen(filename, "r")) == NULL){
		printf("ERROR: Can't open opacity file \"%s\". Check whether it exists.\n", filename);
		exit(EXIT_FAILURE);
	}
	else{
		printf("Opacity file \"%s\" was set.\n", filename);
	}
	fgets(buf, 256, fp);
	op->R[0] = atof(strtok(buf, " "));
	i = 1;
	while(1){
		ascii = strtok(NULL, " ");
		if(ascii == NULL){
			break;
		}
		else{
			op->R[i] = atof(ascii);
			i++;
		}
	}
	op->jmax = i;
	while(fgets(buf, 256, fp) != NULL){
		k = 0;
		op->T[l] = atof(strtok(buf, " "));
		while(1){
			ascii = strtok(NULL, " ");
			if(ascii == NULL){
				break;
			}
			else{
				op->kappa[l*(op->jmax)+k] = atof(ascii);
				k++;
			}
		}
		l++;
	}
	op->imax = l;
	fclose(fp);
}

double kappa_r(double rho, double T)
{
	double a, b, c;
	//double R = rho/pow(T*1.e-06, 3.0);
	double R = rho/pow(T, 1.5);
	double kappa = 0., sigma, sigma_at_R0;
	int i, j;
	static opacity op = {{0.}, {0.}, {0.}, 0, 0};

	if(op.imax == 0){
		set_opacity("./LCFiles/opacity.txt", &op);
	}
	
	i = op.imax/2;
	j = op.jmax/2;

	if(log10(T) < op.T[0]){
		if(log10(R) < op.R[0]){
			return pow(10., op.kappa[0]);
		}
		else if(log10(R) > op.R[0] && log10(R) < op.R[op.jmax-1]){
			j = dcht(log10(R), op.R, op.jmax);
			kappa = op.kappa[j]+(op.kappa[j+1]-op.kappa[j])/(op.R[j+1]-op.R[j])*(log10(R)-op.R[j]);
			return pow(10., kappa);
		}
		else{
			return pow(10., op.kappa[op.jmax-1]);
		}
	}
	else if(log10(R) > op.R[op.jmax-1]){
		i = dcht(log10(T), op.T, op.imax);
		kappa = op.kappa[op.jmax*i+op.jmax-1]
		+(op.kappa[op.jmax*(i+1)+op.jmax-1]-op.kappa[op.jmax*i+op.jmax-1])/(op.T[i+1]-op.T[i])*(log10(T)-op.T[i]);
#if WARN == 1
		printf("WARN (rosseland): Too high density to extrapolate or interpolate opacity...\n");
		printf("density: %e log10(R): %f\n", rho, log10(R));
#endif
		kappa = pow(10., kappa);
		kappa *= R/pow(10., op.R[op.jmax-1]);
		return kappa;
	}
	else if(log10(T) > op.T[op.imax-1]){
		return pow(10., op.kappa[op.jmax*(op.imax-1)]);
	}
	else if(log10(T) > op.T[0] && log10(T) < op.T[op.imax-1] && log10(R) < op.R[0]){
		i = dcht(log10(T), op.T, op.imax);
		kappa = op.kappa[op.jmax*i]+(op.kappa[op.jmax*(i+1)]-op.kappa[op.jmax*i])/(op.T[i+1]-op.T[i])*(log10(T)-op.T[i]);
		sigma = sigma_sc(rho, T);
		sigma_at_R0 = sigma_sc(pow(10,op.R[0])*pow(T, 1.5),T);
		kappa = pow(10., kappa)-sigma_at_R0; //absorption opacity at R=R0.
		kappa = kappa*R*pow(10., -op.R[0]);
		return fmax(kappa+sigma, sigma);
	}
	else{
		i = dcht(log10(T), op.T, op.imax);
		j = dcht(log10(R), op.R, op.jmax);
		a = ((op.kappa[op.jmax*(i+1)+j])-(op.kappa[op.jmax*i+j]))/(op.T[i+1]-op.T[i]);
		b = ((op.kappa[op.jmax*i+j+1])-(op.kappa[op.jmax*i+j]))/(op.R[j+1]-op.R[j]);
		c = ((op.kappa[op.jmax*(i+1)+j+1])-(op.kappa[op.jmax*(i+1)+j])-(op.kappa[op.jmax*i+j+1])
			+(op.kappa[op.jmax*i+j]))/((op.T[i+1]-op.T[i])*(op.R[j+1]-op.R[j]));
		kappa = op.kappa[op.jmax*i+j]+a*(log10(T)-op.T[i])+b*(log10(R)-op.R[j])+c*(log10(T)-op.T[i])*(log10(R)-op.R[j]);
		return pow(10., kappa);
	}
}

double sigma_sc(double rho, double T)
{
	double a, b, c;
	double R = rho/pow(T*1.e-06, 3.);
	double kappa = 0.;
	static opacity op = {{0.}, {0.}, {0.}, 0, 0};
	int i, j;

	if(op.imax == 0){
		set_opacity("./LCFiles/opacity_sc.txt", &op);
	}

	i = op.imax/2;
	j = op.jmax/2;
	if(log10(T) < op.T[0]){
		return 1.e-08;
	}
	else if(log10(R) > op.R[op.jmax-1]){
#if WARN == 1
		printf("WARN (scattering): Too high density to extrapolate or interpolate opacity...\n");
		printf("density: %e\n", rho);
#endif
		i = dcht(log10(T), op.T, op.imax);
		kappa = op.kappa[op.jmax*i+op.jmax-1]
		+(op.kappa[op.jmax*(i+1)+op.jmax-1]-op.kappa[op.jmax*i+op.jmax-1])/(op.T[i+1]-op.T[i])*(log10(T)-op.T[i]);
		return pow(10., kappa);
	}
	else if(log10(T) > op.T[op.imax-1]){
		return sigma_sc_high(rho, T);
	}
	else if(log10(T) > op.T[0] && log10(T) <= op.T[op.imax-1] && log10(R) < (op.R[0]+1.e-05)){
		if(log10(T) < op.T[op.imax-1]){
			i = dcht(log10(T), op.T, op.imax);
		}
		else{
			i = op.imax-2;
		}
		kappa = op.kappa[op.jmax*i]+(op.kappa[op.jmax*(i+1)]-op.kappa[op.jmax*i])/(op.T[i+1]-op.T[i])*(log10(T)-op.T[i]);
		return pow(10., kappa);
	}
	else{
		i = dcht(log10(T), op.T, op.imax);
		j = dcht(log10(R), op.R, op.jmax);
		a = ((op.kappa[op.jmax*(i+1)+j])-(op.kappa[op.jmax*i+j]))/(op.T[i+1]-op.T[i]);
		b = ((op.kappa[op.jmax*i+j+1])-(op.kappa[op.jmax*i+j]))/(op.R[j+1]-op.R[j]);
		c = ((op.kappa[op.jmax*(i+1)+j+1])-(op.kappa[op.jmax*(i+1)+j])-(op.kappa[op.jmax*i+j+1])
			+(op.kappa[op.jmax*i+j]))/((op.T[i+1]-op.T[i])*(op.R[j+1]-op.R[j]));
		return pow(10., op.kappa[op.jmax*i+j]+a*(log10(T)-op.T[i])+b*(log10(R)-op.R[j])+c*(log10(T)-op.T[i])*(log10(R)-op.R[j]));
	}
}

double kappa_p(double rho, double T)
{
	double a, b, c;
	double R = rho/pow(T, 1.5);
	double kappa = 0.;
	static opacity op = {{0.}, {0.}, {0.}, 0, 0};
	int i, j;
	
	if(op.imax == 0){
		set_opacity("./LCFiles/kappa_p.txt", &op);
	}

#if SWITCH_OP_FREE_FREE == 0
	if(log10(T) < op.T[0]){
		if(log10(R) < op.R[0]){
			return pow(10., op.kappa[0]-2.*op.T[0]+op.R[0]);
		}
		else if(log10(R) > op.R[0] && log10(R) < op.R[op.jmax-1]){
			j = dcht(log10(R), op.R, op.jmax);
			kappa = op.kappa[j]+(op.kappa[j+1]-op.kappa[j])/(op.R[j+1]-op.R[j])*(log10(R)-op.R[j]);
			kappa = kappa-3.5*op.T[0]+log10(rho);
			return pow(10., kappa);
		}
		else{
			return pow(10., op.kappa[op.jmax-1]-2.*op.T[0]+op.R[op.jmax-1]);
		}
	}
	else if(log10(R) > op.R[op.jmax-1]){
		i = dcht(log10(T), op.T, op.imax);
		kappa = op.kappa[op.jmax*i+op.jmax-1]
		+(op.kappa[op.jmax*(i+1)+op.jmax-1]-op.kappa[op.jmax*i+op.jmax-1])/(op.T[i+1]-op.T[i])*(log10(T)-op.T[i]);
#if WARN == 1
		printf("WARN (planck_mean): Too high density to extrapolate or interpolate opacity...\n");
		printf("density: %e log10(R): %f\n", rho, log10(R));
#endif
		kappa = kappa-2.*log10(T)+op.R[op.jmax-1];
		kappa = pow(10., kappa);
		kappa *= R/pow(10., op.R[op.jmax-1]);
		return kappa;
	}
	else if(log10(T) > op.T[op.imax-1]){
		return pow(10., op.kappa[op.jmax*(op.imax-1)]-2.*log10(T)+log10(R));
	}
	else if(log10(T) > op.T[0] && log10(T) < op.T[op.imax-1] && log10(R) < op.R[0]){
		i = dcht(log10(T), op.T, op.imax);
		kappa = op.kappa[op.jmax*i]
		+(op.kappa[op.jmax*(i+1)]-op.kappa[op.jmax*i])/(op.T[i+1]-op.T[i])*(log10(T)-op.T[i]);
		kappa = kappa-2.*log10(T)+log10(R);
		return pow(10., kappa);
	}
	else{
		i = dcht(log10(T), op.T, op.imax);
		j = dcht(log10(R), op.R, op.jmax);
		a = ((op.kappa[op.jmax*(i+1)+j])-(op.kappa[op.jmax*i+j]))/(op.T[i+1]-op.T[i]);
		b = ((op.kappa[op.jmax*i+j+1])-(op.kappa[op.jmax*i+j]))/(op.R[j+1]-op.R[j]);
		c = ((op.kappa[op.jmax*(i+1)+j+1])-(op.kappa[op.jmax*(i+1)+j])-(op.kappa[op.jmax*i+j+1])
			+(op.kappa[op.jmax*i+j]))/((op.T[i+1]-op.T[i])*(op.R[j+1]-op.R[j]));
		kappa = op.kappa[op.jmax*i+j]
			+a*(log10(T)-op.T[i])+b*(log10(R)-op.R[j])+c*(log10(T)-op.T[i])*(log10(R)-op.R[j]);
		kappa = kappa-2.*log10(T)+log10(R);
		return pow(10., kappa);
	}

#elif SWITCH_OP_FREE_FREE == 1
	return j_ff(rho, T)/((P_A)*(P_C)*pow(T, 4.));

#else
	return 0.;

#endif


}

double mmw(double rho, double T)
{
	double a, b, c;
	double R = rho/pow(T*1.e-06, 3.);
	double kappa = 0.;
	static opacity op = {{0.}, {0.}, {0.}, 0, 0};
	int i, j;
	
	if(op.imax == 0){
		set_opacity("./LCFiles/mean_ml_wght.txt", &op);
	}
	i = op.imax/2;
	j = op.jmax/2;
	if(log10(T) < op.T[0]){
		return 1.29870;
	}
	else if(log10(R) > op.R[op.jmax-1]){
		i = dcht(log10(T), op.T, op.imax);
		return op.kappa[op.jmax*(i+1)-1];
	}
	else if(log10(T) > op.T[op.imax-1]){
		return op.kappa[op.jmax*op.imax-1];
	}
	else if(log10(T) > op.T[0] && log10(T) < op.T[op.imax-1] && log10(R) < op.R[0]){
		i = dcht(log10(T), op.T, op.imax);
		kappa = (op.kappa[op.jmax*(i+1)]-op.kappa[op.jmax*i])/(op.T[op.jmax*(i+1)]-op.T[op.jmax*i]);
		kappa = kappa*(log10(T)-op.T[op.jmax*i])+op.kappa[op.jmax*i];
		return kappa;
	}
	else{
		i = dcht(log10(T), op.T, op.imax);
		j = dcht(log10(R), op.R, op.jmax);
		a = ((op.kappa[op.jmax*(i+1)+j])-(op.kappa[op.jmax*i+j]))/(op.T[i+1]-op.T[i]);
		b = ((op.kappa[op.jmax*i+j+1])-(op.kappa[op.jmax*i+j]))/(op.R[j+1]-op.R[j]);
		c = ((op.kappa[op.jmax*(i+1)+j+1])-(op.kappa[op.jmax*(i+1)+j])-(op.kappa[op.jmax*i+j+1])
			+(op.kappa[op.jmax*i+j]))/((op.T[i+1]-op.T[i])*(op.R[j+1]-op.R[j]));
		kappa = op.kappa[op.jmax*i+j]+a*(log10(T)-op.T[i])+b*(log10(R)-op.R[j])+c*(log10(T)-op.T[i])*(log10(R)-op.R[j]);
		return op.kappa[op.jmax*i+j]+a*(log10(T)-op.T[i])+b*(log10(R)-op.R[j])+c*(log10(T)-op.T[i])*(log10(R)-op.R[j]);
	}
}

double sigma_sc_high(double rho, double T)
{
	return 0.2*(1.+X)/(1.+pow(T/4.5e+08, 0.86));
}

double kappa_ff(double rho, double T)
{
	return 5.08e+22*(X+Y)*(1.+X)*rho*pow(T, -7./2.);
}

double j_ff(double rho, double T)
{
	return 6.5e+20*(X+0.5*Y)*(X+Y)*rho*sqrt(T);
}

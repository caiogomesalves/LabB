
#include "runMCMC.h"


double logfY_Binorm(mat S, mat Y, mat L)
{
    double res = as_scalar( trans(Y)*S - accu( L%log(exp(S)+1) ) );
    return(res);
}

mat DlogfY_Binorm(mat S, mat Y, mat L)
{
    return( Y - L%exp(S)/(exp(S)+1) );
}

double logfY_Pois(mat S, mat Y, mat L)
{
    double res = as_scalar( trans(Y)*S - accu(L%exp(S)) );
    return(res);
}

mat DlogfY_Pois(mat S, mat Y, mat L)
{
    return( Y-L%exp(S) );
}

double logft_Pois_mean(mat S, mat Y, mat L)
{
	double res = as_scalar( -accu(L%exp(S)) + accu(Y)*log(accu(L%exp(S))) );
	return(res);
}

double logft_Pois_max(mat S, mat Y, mat L)
{
	double t = as_scalar( max(Y, 0) );
	mat lam = L%exp(S);
    double Pt = 1, Pt1 = 1;
    for(int i=0; i<lam.n_rows; i++)
    {
        Pt *= gsl_cdf_poisson_P(t, lam[i,0]);
        Pt1 *= gsl_cdf_poisson_P(t-1, lam[i,0]);
    }
    double res = log(Pt - Pt1);
    if(res < log(0.0001)) res = log(0.0001);
    return(res);
}

double logft_Pois_min(mat S, mat Y, mat L)
{
	double t = as_scalar( min(Y, 0) );
	mat lam = L%exp(S);
    double Pt = 1, Pt1 = 1;
    for(int i=0; i<lam.n_rows; i++)
    {
        Pt *= 1-gsl_cdf_poisson_P(t, lam[i,0]);
        Pt1 *= 1-gsl_cdf_poisson_P(t-1, lam[i,0]);
    }
    double res = log( (1-Pt) - (1-Pt1) );
    if(res < log(0.0001)) res = log(0.0001);
    return(res);
}

double loglSS(mat S, mat Zinv, mat D, mat m)
{
    double res = as_scalar( - 0.5*trans(S-D*m)*Zinv*(S-D*m) );
    return(res);
}

double logll(mat S, mat Zinv, mat D, mat m, mat Z, mat C, mat Omg, double s, double a)
{
    double res = as_scalar( -0.5*log(det(Z)) - 0.5*trans(S-D*m)*Zinv*(S-D*m) + log(det(C)) + log(det(chol(Omg))) + log(a*(2+s)) );
    return(res);
}


double logPiSig_Halft(double s, double A, double df)
{
    return(-(df+1)/2)*log(1 + s*s/(A*A*df));
}

double logPiSig_InvGamma(double s, double shape, double scale)
{
    return(-(shape+1)*log(s) - scale/s);
}

double logPiSig_Recip(double s, double a, double b)
{
    return(-log(s));
}

mat rhoPowExp(mat T, double a, double k)
{
	return( exp(-pow(T/a,k)) );
}

mat rhoSph(mat T, double a, double k)
{
	mat R = ones<mat>(T.n_rows,T.n_cols);
	double u;
	for(int i=0; i < T.n_rows; i++)
	{
		for(int j=0; j < T.n_cols; j++)
		{
			u = T(i,j);
			if(u > 0)
			{
				if(u < a) R(i,j) = 1 - 1.5 * (u/a) + 0.5 * pow(u/a, 3);
				else R(i,j) = 0 ;
			}
		}
	}
	return( R );
}

mat rhoMatern(mat T, double a, double k)
{
	mat R = ones<mat>(T.n_rows,T.n_cols);
	double u;
	for(int i=0; i < T.n_rows; i++)
	{
		for(int j=0; j < T.n_cols; j++)
		{
			u = T(i,j);
			if(u > 0)
			{
				if(u < 600*a) R(i,j) = (pow(2, 1-k)/gsl_sf_gamma(k)) * pow(u/a, k)* gsl_sf_bessel_Knu(k, u/a);
				else R(i,j) = 0 ;
			}
		}
	}
	return( R );
}



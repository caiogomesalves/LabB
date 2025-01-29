
#include "runMCMC.h"


// Define function as extern with RcppExport
RcppExport SEXP runMCMCBPcpp( SEXP Y_, SEXP L_, SEXP T_, SEXP D_, SEXP run_, SEXP nmLan_,
                              SEXP fam_, SEXP famY_, SEXP famSig_, SEXP par1_, SEXP par2_, 
                              SEXP ifkappa_, SEXP scale_, SEXP mscale_,
                              SEXP sscale_, SEXP ascale_, SEXP kscale_, SEXP alow_, SEXP aup_,
                              SEXP mini_, SEXP sini_, SEXP aini_, SEXP kini_) {
// Input format: Y=<mat>n*1, L=<mat>n*1, T=<mat>n*n, run=int;

try {
// Convert R inputs into C++ types
mat Y = as<mat>(Y_);
mat L = as<mat>(L_);
mat T = as<mat>(T_);
mat D = as<mat>(D_);	// new~~~~~~~~~~~~~~~~~// matrix(n*p)
int run = as<int>(run_);
int nmLan = as<int>(nmLan_);
int fam = as<int>(fam_);
int famY = as<int>(famY_);
int famSig = as<int>(famSig_);
double par1 = as<double>(par1_);
double par2 = as<double>(par2_);
int ifkappa = as<int>(ifkappa_);
double scale = as<double>(scale_);
double mscale = as<double>(mscale_);
double sscale = as<double>(sscale_);
double ascale = as<double>(ascale_);
double kscale = as<double>(kscale_);
double alow = as<double>(alow_);
double aup = as<double>(aup_);
mat mini = as<mat>(mini_); 	// new~~~~~~~~~~~~~~~~~// matrix(p*1)
double sini = as<double>(sini_);
double aini = as<double>(aini_);
double kini = as<double>(kini_);

double klower = 0.05, kupper = 1.95;
    mat (*rhofunc) (mat rhoU, double rhoa, double rhok );
    if(fam == 1 )
    {
        rhofunc = &rhoPowExp;
        kupper = 1.95;
    }
	else
	{
	    if(fam == 2)
        {
            rhofunc = &rhoMatern;
            kupper = 4.95;
        }
        else
        {
            rhofunc = &rhoSph;
        }
	}

    double (*logfY) (mat S, mat Y, mat L);
    mat (*DlogfY) (mat S, mat Y, mat L);
    mat Shat = log(Y/L), A = diagmat(Y);
    if(famY == 1 )
    {
        logfY = &logfY_Pois;
        DlogfY = &DlogfY_Pois;
        Shat = log(Y/L);
        A = diagmat(Y);
    }
    else
    {
        logfY = &logfY_Binorm;
        DlogfY = &DlogfY_Binorm;
        Shat = log(Y/L/(1-Y/L));
        A = diagmat(Y%(1-Y/L));
    }

    double (*logPiSig) (double s, double par1, double par2);
    if(famSig == 1 )
    {
        logPiSig = &logPiSig_Halft;
    }
	else
	{
	    if(famSig == 2)
        {
            logPiSig = &logPiSig_InvGamma;
        }
        else
        {
            logPiSig = &logPiSig_Recip;
        }
	}
	
    /* set up tunning and other running parameters*/
    int n = Y.n_rows;
    //mat D = ones<mat>(n,1);

    /* starting values */
    mat m = mini;
    double s = sini, a = aini, k = kini;
    mat S1 = randn<mat>(Y.n_rows,Y.n_cols);
/*  RNGScope scope;
    NumericVector tmp1 = rnorm(n);
    mat S1 = zeros<mat>(1,n);
    for(int i=0;i<n;i++) S1(0,i) = tmp1(i);
*/
    /* Variables initialised by script */
    mat Z = pow(s,2)* (*rhofunc)(T, a, k);
    mat Zinv = inv(Z);
    mat Zt = inv(Zinv + A);
    mat C = chol(Zt);
    mat B = Zt*Zinv*D;	// n*p
    mat E = Zt*A*Shat;
    mat Ot = inv(Z + inv(A));
    mat Omg = inv(trans(D)*Ot*D );
    mat tildem = chol(trans(D)*Ot*D) * ( m - Omg*trans(D)*Ot*Shat );
    //mat tildem = chol(trans(D)*Ot*D) * ( m - inv(trans(D)*Ot*D) * trans(D)*Ot*Shat );
    mat S = C*S1 + E + B*m;

    /* Clear space for saving posterior samples */
    double accS1 = 0, accs = 0, acca = 0, accm = 0, acck = 0;
    mat ssave = zeros<mat>(1,run), asave = ssave, ksave = ssave;
    mat msave = zeros<mat>(m.n_rows,run), Ssave = zeros<mat>(n,run); // tildemsave = msave; Spsave=Ssave;

double snew, anew, knew, logp, logq;
mat mnew, tildemnew, gradm, gradmnew;
mat gradS, gradSnew;
mat S1new = S1, Snew = S;
mat Znew = Z, Zinvnew = Zinv, Ztnew = Zt, Cnew = C, Bnew = B, Enew = E, Otnew = Ot, Omgnew = Omg;


double ruf_S, ruf_m, ruf_s, ruf_a, ruf_k, rnorm_s, rnorm_a, rnorm_k;



/*  Start the main MCMC loop */
for(int i=0; i<run; i++)
{


    /* Gibbs I: update S1|s,a,m,Y with a Langevin step  */
    /* Updates S1 nmLan times within each Gibbs step */
    for(int w=0; w<nmLan; w++)
    {
        gradS = 0.5*( trans(C)* (*DlogfY)(S, Y, L)  - trans(C)*Zinv*C*S1 - trans(C)*Zinv*(B*m + E - D*m) );
        S1new = S1 + scale*gradS + randn<mat>(Y.n_rows,Y.n_cols) *sqrt(scale);
        Snew = C*S1new + B*m + E;
        gradSnew = 0.5*( trans(C)* (*DlogfY)(Snew, Y, L) - trans(C)*Zinv*S1new - trans(C)*Zinv*(B*m + E - D*m) );
        logp = (*logfY)(S, Y, L) + loglSS(S, Zinv, D, m);
        logq = (*logfY)(Snew, Y, L) + loglSS(Snew, Zinv, D, m);
        ruf_S = as_scalar(randu<vec>(1));
        if ( log( ruf_S ) < logq-logp- as_scalar( 0.5*trans(2*(S1new-S1)+scale*(gradSnew-gradS))*(gradSnew+gradS) )  )
            {
                S1 = S1new;
                S = Snew;
                accS1 = accS1+1;
            }
    }




    /* Gibbs II: update m|S1,a,Y,s with an additive Langevin step  */
    gradm = 0.5*( inv(chol(trans(D)*Ot*D)) * ( trans(B)* (*DlogfY)(S, Y, L) - trans(B-D)*Zinv*(S-D*m) ) );
    tildemnew =  tildem + mscale*gradm+ randn<mat>(m.n_rows,m.n_cols)*sqrt(mscale);
    mnew = m + chol(Omg)*(tildemnew-tildem) ;
    Snew = C*S1+B*mnew+E;
    gradmnew = 0.5*( inv(chol(trans(D)*Ot*D))*( trans(B)* (*DlogfY)(Snew, Y, L) - trans(B-D)*Zinv*(Snew-D*mnew) ) );
    logp = (*logfY)(S, Y, L) + loglSS(S, Zinv, D, m);
    logq = (*logfY)(Snew, Y, L) + loglSS(Snew, Zinv, D, mnew);
    ruf_m = as_scalar(randu<vec>(1));
    if ( log( ruf_m ) < logq-logp- as_scalar( 0.5* trans( 2*(tildemnew-tildem)+mscale*(gradmnew-gradm) ) *(gradmnew+gradm) )   )
    {
        tildem = tildemnew;
        m = mnew;
        S = Snew;
        accm = accm+1;
    }
    msave.col(i) = m;



    /* Gibbs III: update s|S1,a,Y,m with a multiplicative RW step */
    rnorm_s = as_scalar(randn<vec>(1));
    snew = s*exp( rnorm_s*sscale );
    anew = a*pow(snew,2)/pow(s,2);

    Znew = pow(snew,2)* (*rhofunc)(T, anew, k);
    Zinvnew = inv(Znew);
    Ztnew = inv(Zinvnew + A);
    Cnew = chol(Ztnew);
    Bnew = Ztnew*Zinvnew*D;
    Enew = Ztnew*A*Shat;
    Otnew = inv(Znew + inv(A));
    Omgnew = inv(trans(D)*Otnew*D );

    mnew =  chol(Omgnew)*tildem + Omgnew*trans(D)*Otnew*Shat ;
    Snew = Cnew*S1+Bnew*mnew+Enew;

    logp = as_scalar( (*logfY)(S, Y, L) + logll(S, Zinv, D, m, Z, C, Omg, s, a) + (*logPiSig)(s, par1, par2) );
    logq = as_scalar( (*logfY)(Snew, Y, L) + logll(Snew, Zinvnew, D, mnew, Znew, Cnew, Omgnew, snew, anew) + (*logPiSig)(snew, par1, par2)  );

    if(anew < aup && anew > alow)
    {
    	ruf_s = as_scalar(randu<vec>(1));
        if ( log( ruf_s ) < logq+log(snew)+log(anew)-logp-log(s)-log(a) )
        {
            s=snew;
            a=anew;
            Z = Znew;
            Zinv = Zinvnew;
            Zt = Ztnew;
            C = Cnew;
            B = Bnew;
            E = Enew;
            Ot = Otnew;
            Omg = Omgnew;
            m = mnew;
            S = Snew;
            accs=accs+1;
        }
    }
    ssave(0,i)=s;


    /* Gibbs IV: update a|S1,s,Y,m with a multiplicative RW step */
    rnorm_a = as_scalar(randn<vec>(1));
    anew = a*exp(-rnorm_a*ascale);
    snew = s;

    Znew = pow(snew,2)* (*rhofunc)(T, anew, k);
    Zinvnew = inv(Znew);
    Ztnew = inv(Zinvnew + A);
    Cnew = chol(Ztnew);
    Bnew = Ztnew*Zinvnew*D;
    Enew = Ztnew*A*Shat;
    Otnew = inv(Znew + inv(A));
    Omgnew = inv(trans(D)*Otnew*D );

    mnew =  chol(Omgnew)*tildem + Omgnew*trans(D)*Otnew*Shat ;
    Snew = Cnew*S1+Bnew*mnew+Enew;

    logp = as_scalar( (*logfY)(S, Y, L) + logll(S, Zinv, D, m, Z, C, Omg, s, a)   );
    logq = as_scalar( (*logfY)(Snew, Y, L) + logll(Snew, Zinvnew, D, mnew, Znew, Cnew,  Omgnew, snew, anew) );

    if(anew < aup && anew > alow)
    {
    	ruf_a = as_scalar(randu<vec>(1));
        if ( log( ruf_a ) < logq+log(anew)-logp-log(a) )
        {
            s=snew;
            a=anew;
            Z = Znew;
            Zinv = Zinvnew;
            Zt = Ztnew;
            C = Cnew;
            B = Bnew;
            E = Enew;
            Ot = Otnew;
            Omg = Omgnew;
            m = mnew;
            S = Snew;
            acca=acca+1;
        }
    }
    asave(0,i)=a;



    /* Gibbs V: update k|S1,s,Y,m,a with a multiplicative RW step */
if(ifkappa!=0)
{
    rnorm_k = as_scalar(randn<vec>(1));
    knew = k + rnorm_k*kscale;
        while(knew < klower || knew > kupper)
        {
            rnorm_k = as_scalar(randn<vec>(1));
            knew = k + rnorm_k*kscale;
        }

    Znew = pow(s,2)* (*rhofunc)(T, a, knew);
    Zinvnew = inv(Znew);
    Ztnew = inv(Zinvnew + A);
    Cnew = chol(Ztnew);
    Bnew = Ztnew*Zinvnew*D;
    Enew = Ztnew*A*Shat;
    Otnew = inv(Znew + inv(A));
    Omgnew = inv(trans(D)*Otnew*D );

    mnew =  chol(Omgnew)*tildem + Omgnew*trans(D)*Otnew*Shat ;
    Snew = Cnew*S1+Bnew*mnew+Enew;

    logp = as_scalar( (*logfY)(S, Y, L) + logll(S, Zinv, D, m, Z, C, Omg, s, a) );
    logq = as_scalar( (*logfY)(Snew, Y, L) + logll(Snew, Zinvnew, D, mnew, Znew, Cnew, Omgnew, snew, anew) );

    ruf_k = as_scalar(randu<vec>(1));
        if ( log( ruf_k ) < logq - logp )
        {
            k = knew;
            Z = Znew;
            Zinv = Zinvnew;
            Zt = Ztnew;
            C = Cnew;
            B = Bnew;
            E = Enew;
            Ot = Otnew;
            Omg = Omgnew;
            m = mnew;
            S = Snew;
            acck=acck+1;
        }
    ksave(0,i)=k;
}

    Ssave.col(i) = S;
}


if(ifkappa!=0)
{
NumericVector acc = NumericVector::create( _["accS1"] = accS1/run/nmLan,
	_["accm"] = accm/run, _["accs"] = accs/run, _["acca"] = acca/run, _["acck"] = acck/run );

return List::create( _["S.posterior"] = Ssave,
	_["m.posterior"] = msave,
	_["s.posterior"] = ssave,
	_["a.posterior"] = asave,
	_["k.posterior"] = ksave,
	_["AccRate"] = acc );
}
else
{
NumericVector acc = NumericVector::create( _["accS1"] = accS1/run/nmLan,
	_["accm"] = accm/run, _["accs"] = accs/run, _["acca"] = acca/run );

return List::create( _["S.posterior"] = Ssave,
	_["m.posterior"] = msave,
	_["s.posterior"] = ssave,
	_["a.posterior"] = asave,
	_["AccRate"] = acc );
}


} catch( std::exception &ex ){
	forward_exception_to_r( ex );
} catch(...){
	::Rf_error("c++ exception (unknown reason)");
}
return R_NilValue;
}




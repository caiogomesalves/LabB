#ifndef runMCMC_h
#define runMCMC_h

// Add C++ code to pkg src/ directory.
// #include <Rcpp.h>
#include <RcppArmadillo.h>
//#include <RcppGSL.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>

using namespace Rcpp;
using namespace arma;

RcppExport SEXP myPackage_hello_world() ;

double logfY_Binorm(mat S, mat Y, mat L);
mat DlogfY_Binorm(mat S, mat Y, mat L);

double logfY_Pois(mat S, mat Y, mat L);
mat DlogfY_Pois(mat S, mat Y, mat L);

double logft_Pois_mean(mat S, mat Y, mat L);
double logft_Pois_max(mat S, mat Y, mat L);
double logft_Pois_min(mat S, mat Y, mat L);

double loglSS(mat S, mat Zinv, mat D, mat m);
double logll(mat S, mat Zinv, mat D, mat m, mat Z, mat C, mat Omg, double s, double a);

double logPiSig_Halft(double s, double A, double df);
double logPiSig_InvGamma(double s, double shape, double scale);
double logPiSig_Recip(double s, double a, double b);

mat rhoPowExp(mat T, double a, double k);
mat rhoSph(mat T, double a, double k);
mat rhoMatern(mat T, double a, double k);

RcppExport SEXP runMCMCBPcpp( SEXP Y_, SEXP L_, SEXP T_, SEXP D_, SEXP run_, SEXP nmLan_,
                              SEXP fam_, SEXP famY_, SEXP famSig_, SEXP par1_, SEXP par2_,
                              SEXP ifkappa_, SEXP scale_, SEXP mscale_,
                              SEXP sscale_, SEXP ascale_, SEXP kscale_, SEXP alow_, SEXP aup_,
                              SEXP mini_, SEXP sini_, SEXP aini_, SEXP kini_);

RcppExport SEXP runMCMCpartialPoiscpp( SEXP Y_, SEXP L_, SEXP T_, SEXP D_, SEXP run_, SEXP nmLan_,
                              SEXP fam_, SEXP famY_, SEXP famT_, SEXP ifkappa_, SEXP scale_, SEXP mscale_,
                              SEXP sscale_, SEXP ascale_, SEXP kscale_, SEXP alow_, SEXP aup_,
                              SEXP mini_, SEXP sini_, SEXP aini_, SEXP kini_);

#endif

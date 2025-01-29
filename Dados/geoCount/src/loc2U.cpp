
#include "runMCMC.h"


// Define function as extern with RcppExport
RcppExport SEXP loc2Ucpp( SEXP loc_) {
try {
// Convert R inputs into C++ types
mat loc = as<mat>(loc_);

int n = loc.n_rows;
mat U = zeros<mat>(n,n);

for(int i=0; i<n-1; i++)
{
	for(int j=i+1; j<n; j++)
	{
		U(i,j) = sqrt( pow(loc(i,0)-loc(j,0),2) + pow(loc(i,1)-loc(j,1),2) );
		U(j,i) = U(i,j);
	}
}

return( wrap(U) );

} catch( std::exception &ex ){
	forward_exception_to_r( ex );
} catch(...){
	::Rf_error("c++ exception (unknown reason)");
}
return R_NilValue;
}



// Define function as extern with RcppExport
RcppExport SEXP locUloccpp( SEXP loc_, SEXP locp_ ) {
try {
// Convert R inputs into C++ types
mat loc = as<mat>(loc_);
mat locp = as<mat>(locp_);

int n = loc.n_rows, np = locp.n_rows;
mat U = zeros<mat>(np,n);

for(int i=0; i<np; i++)
{
	for(int j=0; j<n; j++)
	{
		U(i,j) = sqrt( pow(locp(i,0)-loc(j,0),2) + pow(locp(i,1)-loc(j,1),2) );
	}
}

return( wrap(U) );

} catch( std::exception &ex ){
	forward_exception_to_r( ex );
} catch(...){
	::Rf_error("c++ exception (unknown reason)");
}
return R_NilValue;
}






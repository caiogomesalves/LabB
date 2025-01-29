
#include "runMCMC.h"

// Define function as extern with RcppExport
RcppExport SEXP helloWorldcpp() {
try {
    CharacterVector x = CharacterVector::create( "Hello, world!", "Loading is successful!" )  ;
    NumericVector y   = NumericVector::create( 0, 1 ) ;
    List z            = List::create( x, y ) ;

    return z ;

} catch( std::exception &ex ){
	forward_exception_to_r( ex );
} catch(...){
	::Rf_error("c++ exception (unknown reason)");
}
return R_NilValue;
}




#include <Rcpp.h>
using namespace Rcpp;

#include <sosFilter.c>


// [[Rcpp::export]]
NumericVector sosFilter_cpp(NumericMatrix sos, double gain, NumericVector x, NumericMatrix zi) 
{
    int nx = x.size();
    int nsec = sos.nrow();
    
    sos = transpose(sos);
    zi = transpose(zi);
    
    NumericVector y = NumericVector(nx);
    
    for(int i = 0; i < nx; i++)
        y[i] = sosFilter_next(x[i], &zi[0], &sos[0], nsec, gain);
    
    return y;
}




// [[Rcpp::export]]
NumericVector UpSosDn_cpp(NumericVector x, NumericMatrix sos, double gain, int p, int q)
{
    int nx = x.size();
    int nsec = sos.nrow();
    sos = transpose(sos);
    int ny = nx * p / q;
    NumericVector y = NumericVector(ny);
    upSosDn(&x[0], nx, &sos[0], nsec, gain, p, q, &y[0]); 
    return y;

}
  
  
  
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/

#include <Rcpp.h>
using namespace Rcpp;



double sectionFormIItrans(double *state, double section[6], double value)
{
  double *b = section;
  double *a = section + 3;
  double res = b[0] * value + state[0];
  state[0] = b[1] * value - a[1] * res + state[1];
  state[1] = b[2] * value - a[2] * res;
  return res;
}


double sosFilter_next(double x, double *state, double *sos, int nsec, double gain) 
{
    for (int r = 0; r < nsec; r++)
      x = sectionFormIItrans(state + r * 2, sos + r * 6, x);
    return x * gain;
}



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


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/

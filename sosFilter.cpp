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


void sosFilter_c(double *x, int nx, double *state, double *sos, int nsec, double gain, double *y) 
{
  for(int i = 0; i < nx; i++)
  {          
    double value = x[i];
    for (int r = 0; r < nsec; r++)
       value = sectionFormIItrans(state + r * 2, sos + r * 6, value);
    y[i] = value * gain;
  }
  
}




// [[Rcpp::export]]
NumericVector sosFilter_cpp(NumericMatrix sos, double gain, NumericVector x, NumericMatrix zi) 
{
  int nx = x.size();
  int nsec = sos.nrow();
  
  sos = transpose(sos);
  zi = transpose(zi);

  NumericVector y = NumericVector(nx);
  
  sosFilter_c(&x[0], nx, &zi[0], &sos[0], nsec, gain, &y[0]);

  return y;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/

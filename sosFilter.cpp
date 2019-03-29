#include <Rcpp.h>
using namespace Rcpp;



double sectionFormIItrans(double *state, double section[6], double value)
{
  double *b = section;
  double *a = section + 3;
  double res = b[0] * value + state[1];
  state[1] = b[1] * value - a[1] * res + state[2];
  state[2] = b[2] * value - a[2] * res;
  
  return res;
}


void sosFilter_c(double *sos, int nsec, double gain, double *x, int nx, double *zi, double *y) 
{
  
  double stateVec[nsec * 3];
  double* state = &stateVec[0];
  
  for (int r = 0; r < nsec; r++)
    for(int i = 0; i < 2; i++)
      state[r * 3 + 1 + i] = zi == 0 ? 0 : zi[r * 2 + i];
  
  for(int i = 0; i < nx; i++)
  {          
    double value = x[i];
    for (int r = 0; r < nsec; r++)
       value = sectionFormIItrans(state + r * 3, sos + r * 6, value);
    y[i] = value * gain;
  }
  
}




// [[Rcpp::export]]
NumericVector sosFilter_cpp(NumericMatrix sos, double gain, NumericVector x, NumericMatrix zi) 
{
  int nx = x.size();
  
  NumericVector y = NumericVector(nx);
  
  int nsec = sos.nrow();

  sos = transpose(sos);
  zi = transpose(zi);
  
  NumericVector stateVec = NumericVector(nsec * 3);

  sosFilter_c(&sos[0], nsec, gain, &x[0], nx, &zi[0], &y[0]);

  return y;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/

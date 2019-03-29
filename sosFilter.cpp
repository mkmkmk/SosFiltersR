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


// [[Rcpp::export]]
NumericVector UpSosDn_cpp(NumericVector x, NumericMatrix sos, double gain, int p, int q)
{
    int nx = x.size();
    int nsec = sos.nrow();
    
    sos = transpose(sos);
    
    double state[nsec * 2];
    for(int i = 0; i < nsec * 2; i++)
        state[i] = 0;
    
    int ny = nx * p / q;
    NumericVector y = NumericVector(ny);
    
    int iy = 0;
    int iq = 0;
    
    for(int i = 0; i < nx; i++)
    {
        double samp = x[i];
        for(int j = 0; j < p; j++)
        {
           double res = sosFilter_next(samp, state, &sos[0], nsec, gain);
           iq++;
           if(iq == q)
           {
              iq = 0;
              y[iy] = res;
              iy++;
              if(iy == ny)
                 break;
           }
        }
        if(iy == ny)
           break;
    }
    
    return y;

}
  
  
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/

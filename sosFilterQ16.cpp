#include <Rcpp.h>
using namespace Rcpp;

#include <sosFilterQ16.c>
 

// [[Rcpp::export]]
IntegerVector sosFilterQ16_cpp(NumericMatrix sos, NumericVector gains, IntegerVector x, NumericMatrix zi) 
{
    int nx = x.size();
    int nsec = sos.nrow();
     
    sos = transpose(sos);
    zi = transpose(zi);
    
    IntegerVector y = IntegerVector(nx);

    int32_t sosQ16[sos.size()];
    for(int i = 0; i < sos.size(); i++)
        sosQ16[i] = sos[i] * (1ll << Q_BITS);
    
    int32_t ziQ16[zi.size()];
    for(int i = 0; i < zi.size(); i++)
        ziQ16[i] = zi[i] * (1ll << Q_BITS);
    
    int32_t gainsQ16[gains.size()];
    for(int i = 0; i < gains.size(); i++)
        gainsQ16[i] = gains[i] * (1ll << Q_BITS);
    
    if(0) 
    {
        printf("sos.size() = %d \n", (int)sos.size());
        for(int i = 0; i < sos.size(); i++)
            printf("sosQ16[i] = %d \n", (int)sosQ16[i]);
        for(int i = 0; i < zi.size(); i++)
            printf("ziQ16[i] = %d \n", (int)ziQ16[i]);
    }
    
    for(int i = 0; i < nx; i++)
        y[i] = sosFilterQ16_next(x[i], ziQ16, sosQ16, nsec, gainsQ16);
    
   
    return y;
}
   
 

// [[Rcpp::export]]
IntegerVector UpSosDnQ16_cpp(IntegerVector x, NumericMatrix sos, int p, int q)
{
    int nx = x.size();
    int nsec = sos.nrow();
    sos = transpose(sos);
    int ny = nx * p / q;
    IntegerVector y = IntegerVector(ny);
    
    int32_t state[nsec * 2];
    for(int i = 0; i < nsec * 2; i++)
        state[i] = 0;
    
    int32_t sosQ16[sos.size()];
    for(int i = 0; i < sos.size(); i++)
        sosQ16[i] = sos[i] * (1ll << Q_BITS);
  
    upSosDnQ16(&x[0], nx, sosQ16, nsec, 0, state, p, q, &y[0]); 
    return y;

}

// [[Rcpp::export]]
IntegerVector UpSosDnQ16_g_cpp(IntegerVector x, NumericMatrix sos, NumericVector gains, int p, int q)
{
    int nx = x.size();
    int nsec = sos.nrow();
    sos = transpose(sos);
    int ny = nx * p / q;
    IntegerVector y = IntegerVector(ny);
    
    int32_t state[nsec * 2];
    for(int i = 0; i < nsec * 2; i++)
        state[i] = 0;
    
    int32_t sosQ16[sos.size()];
    for(int i = 0; i < sos.size(); i++)
        sosQ16[i] = sos[i] * (1ll << Q_BITS);
    
    int32_t gainsQ16[gains.size()];
    for(int i = 0; i < gains.size(); i++)
        gainsQ16[i] = gains[i] * (1ll << Q_BITS);
    
    
    upSosDnQ16(&x[0], nx, sosQ16, nsec, gainsQ16, state, p, q, &y[0]); 
    return y;
    
}
  
  
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

library(signal)
library(Rcpp)
source("sos.R")

{
    fs = 1000
    
    x = (100 + 1:200 %% 30) * 10
    
    
    tf = butter(2, c(10, 80) / (fs / 2), type = "pass")
    sos = tf2sos(tf)
    
    ft_ref = filter(tf, x, init.x = rep(x[1], length(tf$b) - 1))
    plot(ft_ref, type='l', col = 'blue4')
    
    #ft_flo = sosFilter_cpp(sos$sos, sos$g, x, zi = sos_zi(sos) * x[1])
    #lines(ft_flo, type='l', col = 'blue')
    #stopifnot(all(log2(abs(as.vector(sos_zi(sos) * x[1]) * 2^16)) < 31))
    
    
    #gains = rep(sos$g^(1/nrow(sos$sos)), nrow(sos$sos))
    gains = rep(1, nrow(sos$sos))
    gains[nrow(sos$sos)] = sos$g
    ft_ver = sosFilterQ16_cpp(sos$sos, gains, x, zi = sos_zi(sos) * x[1])
    
    lines(ft_ver, type='l', col = 'red4')
}

# -----------------------

{
    fsamp = 200
    fsig = 3.7
    
    t = 1:(fsamp * 1.1) / fsamp
    x = 200 * ( sin(t * 2 * pi * fsig) + 0.4 * sin(t * 2 * pi * fsamp * .27))
    tf = butter(6, 0.8 / 4)
    sos = tf2sos(tf)
    plot(x, type='l', col = 'blue4')
    
    ft_ref = filter(tf, x)
    plot(ft_ref, type='l', col = 'blue4')
    
    gains = rep(sos$g ^ (1 / nrow(sos$sos)), nrow(sos$sos))
    stopifnot(all(2^16 * gains >= 1))
        
    ft_ver = sosFilterQ16_cpp(sos$sos, gains, x, zi = 0 * sos_zi(sos))
    lines(ft_ver, type = 'l', col = 'red4')
    
    
    
}

*/

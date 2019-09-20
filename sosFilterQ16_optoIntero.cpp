#include <Rcpp.h>
using namespace Rcpp;

#include <sosFilterQ16.c>
 

// [[Rcpp::export]]
IntegerVector sosFilterQ16_direct_cpp(IntegerVector sos, int nsec, IntegerVector x) 
{
    int nx = x.size();

    IntegerVector y = IntegerVector(nx);
    
    int stateSize = nsec * 2;
    int32_t state[stateSize];
    for(int i = 0; i < stateSize; i++)
        state[i] = 0;
    
    for(int i = 0; i < nx; i++)
        y[i] = sosFilterQ16_next(x[i], state, &sos[0], nsec, 0);

    return y;
}
 
  
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R


if(T)
{

    if(T)
    {
        # // Q.16 resamp filter 1/10, order 4, gain x1
        # #define RESAMP_P (1)
        # #define RESAMP_Q (10)
        # #define N_SEC (2)
        # int32_t sos[] = { 887, 1774, 887, 65536, -103234, 41047, 887, 1774, 887, 65536, -115921, 54145 };
        nsec = 2
        sos = c( 887, 1774, 887, 65536, -103234, 41047, 887, 1774, 887, 65536, -115921, 54145 )
    }
    

    if(F)
    {
        # // Q.16 resamp filter 1/10, order 6, gain x1
        # #define RESAMP_P (1)
        # #define RESAMP_Q (10)
        # #define SOS_N_SEC (3)
        # int32_t SOS[] = { 889, 1778, 889, 65536, -102365, 40149, 889, 1778, 889, 65536, -107968, 45934, 889, 1778, 889, 65536, -119277, 57610 };
        nsec = 3
        sos = c( 889, 1778, 889, 65536, -102365, 40149, 889, 1778, 889, 65536, -107968, 45934, 889, 1778, 889, 65536, -119277, 57610 )
    }
    
    if(F)
    {
        # // Q.16 resamp filter 1/10, order 8, gain x1
        # #define RESAMP_P (1)
        # #define RESAMP_Q (10)
        # #define SOS_N_SEC (4)
        # int32_t SOS[] = { 890, 1780, 890, 65536, -102060, 39835, 890, 1780, 890, 65536, -105201, 43077, 890, 1780, 890, 65536, -111543, 49625, 890, 1780, 890, 65536, -121080, 59471 };        
        nsec = 4
        sos = c(890, 1780, 890, 65536, -102060, 39835, 890, 1780, 890, 65536, -105201, 43077, 890, 1780, 890, 65536, -111543, 49625, 890, 1780, 890, 65536, -121080, 59471)
    }
    
    sig = c( 46349, 45978, 46566, 45964, 45743, 44980, 46402, 48867, 45838, 51432, 46220, 45183, 46344, 47163, 47808, 44410, 42837, 45190, 48560, 48712, 49516, 46797, 47468, 45661, 46804, 43527, 45990, 37271, 44083, 40125, 53000, 42629, 47534)
    sig = rep(sig, 10)
    sig = sample(sig)
    
    #sig = 47200 - 10000 * runif(500)
    
    code0 = 47200
    fsig = sosFilterQ16_direct_cpp(sos, nsec, sig - code0) + code0
    #fsig = fsig + code0
    if(F)
        debugonce(filter)
    fsigMA = filter(rep(1/10, 10), a = 1, x=sig)
    
    
    
    plot(sig, type = 'l', col='blue2') # , ylim = range(sig, fsig, fsigMA))
    lines(fsigMA, type = 'l', col='darkgreen', lwd=2)
    lines(fsig, type = 'l', col='red3', lwd=2)
    
    code2microAmp = function(code)
    {
        Vref = 4.096
        Rf = 2200
        Vcom = 2.95
        (Vcom - Vref * code / 2^16) / Rf * 1e6
    }
    plot(code2microAmp(sig), type = 'l', col='blue2') #, ylim = range(sig, fsig))
    lines(code2microAmp(fsigMA), type = 'l', col='darkgreen', lwd=2)
    lines(code2microAmp(fsig), type = 'l', col='red4')
  
    
}

*/

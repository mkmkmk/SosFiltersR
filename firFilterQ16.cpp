#include <Rcpp.h>
using namespace Rcpp;


#define Q_BITS 16

int firFilterQ16_next(int x, int *state, int *s_idx, int *firQ16, int fir_siz) 
{
    int res = 0;
    int *istate = state + *s_idx;
    
    firQ16 += fir_siz;
    for (int i = 0; i < fir_siz; i++)
    {
        res += ((int64_t)*--firQ16 * *istate) >> Q_BITS;
        if(++istate == state + fir_siz)
            istate = state;
    }
    *istate = x << Q_BITS;
    
    if (++*s_idx == fir_siz) 
        *s_idx = 0;

    return res >> Q_BITS;
}

// [[Rcpp::export]]
IntegerVector firFilterQ16_cpp(IntegerVector firQ16, IntegerVector x) 
{
    int nx = x.size();
    int fir_siz = firQ16.size();
    
    IntegerVector y = IntegerVector(nx);
    
    int s_idx = 0;
    int state[fir_siz];
    for(int i = 0; i < fir_siz; i++)
        state[i] = 0;
    
    for(int i = 0; i < nx; i++)
        y[i] = firFilterQ16_next(x[i], state, &s_idx, &firQ16[0], fir_siz);

    return y;
}
 

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

library(signal)
library(Rcpp)

# -----------------------
if(F)
{
    fsamp = 200
    fsig = 3.7
   
    t = 1:(fsamp * 1.1) / fsamp
    x = 200 * ( sin(t * 2 * pi * fsig) + 0.4 * sin(t * 2 * pi * fsamp * .27))
    x = round(x)
    # Py:
    # down = 4
    # h = signal.firwin(2 * 10 * down + 1, .8 / down, window=('kaiser', 5.0))
    # h_Q32 = np.round(h * 2**16).astype(int)
    fir = c( 0,   -15,   -33,   -42,   -32,     0,    48,    94,   111,
             81,     0,  -109,  -204,  -234,  -165,     0,   212,   387,
             436,   303,     0,  -380,  -687,  -767,  -529,     0,   662,
             1199,  1348,   939,     0, -1216, -2268, -2648, -1944,     0,
             3000,  6536,  9872, 12254, 13117, 12254,  9872,  6536,  3000,
             0, -1944, -2648, -2268, -1216,     0,   939,  1348,  1199,
             662,     0,  -529,  -767,  -687,  -380,     0,   303,   436,
             387,   212,     0,  -165,  -234,  -204,  -109,     0,    81,
             111,    94,    48,     0,   -32,   -42,   -33,   -15,     0)
    
    plot(x, type='l', col = 'blue4')

    ft_ref = filter(fir / 2^16, 1, x)
    ft_ref = as.vector(ft_ref)
    ft_ref = ts(ft_ref, start = 2)
    plot(ft_ref, type = 'l', col = 'blue4')
    
        
    ft_ver = firFilterQ16_cpp(fir, x)
    lines(ft_ver, type = 'l', col = 'red4')
    
}

if(F)
{
    
    sig = 1e3 + 5 * runif(100)
    sig = round(sig)
    # Py:
    # down = 4
    # h = signal.firwin(2 * 10 * down + 1, .8 / down, window=('kaiser', 5.0))
    # h_Q32 = np.round(h * 2**16).astype(int)
    fir = c( 0,   -15,   -33,   -42,   -32,     0,    48,    94,   111,
             81,     0,  -109,  -204,  -234,  -165,     0,   212,   387,
             436,   303,     0,  -380,  -687,  -767,  -529,     0,   662,
             1199,  1348,   939,     0, -1216, -2268, -2648, -1944,     0,
             3000,  6536,  9872, 12254, 13117, 12254,  9872,  6536,  3000,
             0, -1944, -2648, -2268, -1216,     0,   939,  1348,  1199,
             662,     0,  -529,  -767,  -687,  -380,     0,   303,   436,
             387,   212,     0,  -165,  -234,  -204,  -109,     0,    81,
             111,    94,    48,     0,   -32,   -42,   -33,   -15,     0)
    
    ft_ref = filter(fir / 2^16, 1, sig)
    ft_ref = as.vector(ft_ref)
    ft_ref = ts(ft_ref, start = 2)
    plot(ft_ref, type = 'l', col = 'blue4')
    
    ft_ver = firFilterQ16_cpp(fir, sig)
    lines(ft_ver, type = 'l', col = 'red4')
    
    if(F)
        print(system.time(firFilterQ16_cpp(fir, 1e3 + 5 * runif(1e7))))
    
}

if(F)
{
    fsamp = 100
    fsig = 3.7
    t = 1:(fsamp * 1.0) / fsamp
    x = 200 * ( sin(t * 2 * pi * fsig) + 1 * sin(t * 2 * pi * fsamp * .27))
    sig = 400 + runif(100) + x
    sig = round(sig)
    h = rep(1/10, 10) + 1:10 / 50
    h = round(h * 2^16)
    
    ft_ref = filter(h / 2^16, a = 1, x = sig)
    ft_ref = as.vector(ft_ref)
    ft_ref = round(ft_ref)
    ft_ref = ts(ft_ref, start = 2)
    plot(ft_ref, type='l', col = 'blue4')

    ft_ver = firFilterQ16_cpp((h), sig)
    lines(ft_ver, type = 'l', col = 'red4')
    
    
}

*/

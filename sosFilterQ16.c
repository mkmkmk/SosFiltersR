
#include <stdio.h>
#include <stdint.h>


#define Q_BITS 16


int32_t sectionFormIItrans(int32_t *state, int32_t section[6], int64_t value)
{
    int32_t *b = section;
    int32_t *a = section + 3;
    int64_t res = ((b[0] * value) >> Q_BITS) + state[0];
    state[0] = ((b[1] * value) >> Q_BITS) - ((a[1] * res) >> Q_BITS) + state[1];
    state[1] = ((b[2] * value) >> Q_BITS) - ((a[2] * res) >> Q_BITS);
    return res;
}


int32_t sosFilterQ16_next(int32_t x, int32_t *state, int32_t *sos, int nsec, int32_t* gainsQ16) 
{
    int64_t ix = (int64_t)x << Q_BITS;
    for (int r = 0; r < nsec; r++)
    {
        ix = sectionFormIItrans(state + r * 2, sos + r * 6, ix);
        ix = (ix * gainsQ16[r]) >> Q_BITS;
    }
    //return (ix * gainQ16) >> (Q_BITS * 2);
    return ix >> Q_BITS;
}


void upSosDnQ16(int32_t *x, int nx, int32_t *sos, int nsec, int32_t *gainsQ16, int32_t *state, int p, int q, int32_t *y) 
{
    int ny = nx * p / q;
    
    int iy = 0;
    int iq = 0;
    
    for(int i = 0; i < nx; i++)
    {
        int32_t samp = x[i];
        for(int j = 0; j < p; j++)
        {
            int32_t res = sosFilterQ16_next(samp, state, sos, nsec, gainsQ16);
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
    
}


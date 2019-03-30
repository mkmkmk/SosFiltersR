


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




void upSosDn(double *x, int nx, double *sos, int nsec, double gain, double *state, int p, int q, double *y) 
{
      int ny = nx * p / q;
    
    int iy = 0;
    int iq = 0;
    
    
    for(int i = 0; i < nx; i++)
    {
        double samp = x[i];
        for(int j = 0; j < p; j++)
        {
            double res = sosFilter_next(samp, state, sos, nsec, gain);
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


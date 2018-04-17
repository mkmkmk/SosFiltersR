
cplxpair = function(inp)
{
    thresh = 100 * .Machine$double.eps
    inpCp = inp[abs(Im(inp)) > thresh]
    inpRe = inp[abs(Im(inp)) < thresh]
    c(inpCp[order(Re(inpCp), Im(inpCp))], inpRe)
}


# matrix(cplxpair(c(exp (2i*pi*0:4/5), 1:6)))
# cplxpair (exp (2i*pi*0:4/5)) == exp (2i*pi*c(3, 2, 4, 1, 0)/5)


cplxreal = function (z, thresh = 100*.Machine$double.eps)
{
    ## interesting for testing: if nargin<2, thresh=1E-3; }
    if (length(z) == 0)
    {
        return(list(zc = c(), zr = c()))
    }else
    {
        zcp = cplxpair(z) # sort complex pairs, real roots at end
        nz = length(z)
        nzrsec = 0
        i = nz
        while (i > 0) # determine no. of real values
        {
            if(abs(Im(zcp[i])) > thresh) break
            zcp[i] = Re(zcp[i])
            nzrsec = nzrsec + 1
            i = i - 1
        }
        nzsect2 = nz - nzrsec
        if ((nzsect2 %% 2) != 0)
        {
            error('cplxreal: Odd number of complex values!');
        }
        nzsec = nzsect2 / 2
        
        if(nzsect2 > 2)
            zc = zcp[seq(2,nzsect2, 2)]
        else
            zc = NULL
        
        if(nzsect2 < nz)
            zr = zcp[(nzsect2 + 1):nz]
        else
            zr = NULL
    }
    
    list(zc = zc, zr = zr)
    
}


# cplxreal(roots(c(1, 0, 0, 1)))
# 0.5+0.8660254i 
# -1+0i



zp2sos = function(z,p,k)
{
    #if nargin<3, k=1; endif
    #if nargin<2, p=[]; endif
    
    zcpl = cplxreal(z)
    zc = zcpl$zc
    zr = zcpl$zr
    
    pcpl = cplxreal(p)
    pc = pcpl$zc
    pr = pcpl$zr
    
    ## zc,zr,pc,pr
    
    nzc = length(zc)
    npc = length(pc)
    
    nzr = length(zr)
    npr = length(pr)
    
    ## Pair up real zeros:
    if (nzr)
    {
        if ((nzr %% 2)==1)
        {
            zr = c(zr, 0)
            nzr = nzr + 1
        }      
        
        nzrsec = nzr / 2
        zrms = -zr[seq(1, nzr - 1, 2)] - zr[seq(2, nzr, 2)]
        zrp = zr[seq(1, nzr - 1, 2)] * zr[seq(2, nzr, 2)]
    }else{
        nzrsec = 0
    }
    
    ## Pair up real poles:
    if (npr)
    {
        if ((npr %% 2)==1){
            pr = c(pr, 0)
            npr = npr + 1 
        }
        nprsec = npr / 2
        prms = -pr[seq(1,npr-1,2)] - pr[seq(2,npr,2)]
        prp = pr[seq(1,npr-1,2)] * pr[seq(2,npr,2)]
    }else{
        nprsec = 0;
    }
    
    nsecs = max(nzc + nzrsec, npc + nprsec)
    
    ## Convert complex zeros and poles to real 2nd-order section form:
    
    if(nzc > 0)
    {
        zcm2r = -2 * Re(zc)
        zca2 = abs(zc) ^ 2
    }else
    {
        zcm2r = NULL
        zca2 = NULL
    }
    
    if(npc > 0)
    {
        pcm2r = -2 * Re(pc)
        pca2 = abs(pc) ^ 2
    }else
    {
        pcm2r = NULL
        pca2 = NULL
    }
    
    
    sos = matrix(rep(0, nsecs*6),ncol=6)
    sos[,1] = rep(1,nsecs) # all 2nd-order polynomials are monic
    sos[,4] = rep(1,nsecs)
    
    nzrl = nzc + nzrsec # index of last real zero section
    nprl = npc + nprsec # index of last real pole section
    
    for (i in 1:nsecs)
    {
        
        if (i<=nzc) # lay down a complex zero pair:
        {
            
            sos[i,2:3] = c(zcm2r[i], zca2[i])
            
        }else if (i<=nzrl) # lay down a pair of real zeros:
        {
            sos[i,2:3] = c(zrms[i-nzc], zrp[i-nzc])
        }
        
        if (i<=npc) # lay down a complex pole pair:
        {
            sos[i,5:6] = c(pcm2r[i], pca2[i])
        }
        else if (i<=nprl) # lay down a pair of real poles:
        {
            sos[i,5:6] = c(prms[i-npc], prp[i-npc])
        }
    }
    
    list(sos = sos, g = k) 
    
}


tf2sos = function(tf)
{
    z = roots(tf$b)
    p = roots(tf$a)
    k = tf$b[1] / tf$a[1]
    zp2sos(z, p, k)
}


#ret val.
1




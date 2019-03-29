rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
print(getwd())
library(signal)
library(Rcpp)
source("sos.R")

sourceCpp('sosFilter.cpp')

#--------------- up_sos_dn test

fsamp = 200
fsig = 3.7


if (T)
{
    p = 2
    q = 1
} else
{
    p = 3
    q = 20
}


t = 1:(fsamp * 1.5) / fsamp
sig = sin(t * 2 * pi * fsig) + 0.4 * sin(t * 2 * pi * fsamp * .27)

plot(t, sig, type = 'l', col = "blue")

#dec_tf = cheby1(8, 0.01, 0.8 / q)
dec_tf = butter(8, 0.8 * min(1, p / q))
dec_sos = tf2sos(dec_tf)

if(F)
{
    #ref_ftSig = hackFilter(dec_tf$b, dec_tf$a, sig, rep(0,100))
    ref_ftSig = sosFilter(dec_sos, sig)
    lines(t[1:length(ref_ftSig)], ref_ftSig, type = 'l', col = "red4")
}

if(F)
    paste(as.vector(t(dec_sos$sos)), collapse = ", ")

resSig = UpSosDn_cpp(sig, dec_sos$sos, dec_sos$g, p, q)

t2 = 0:(length(resSig) - 1) / fsamp * length(sig) / length(resSig)
lines(t2, resSig, type = 'l', col = "red4")


if(F)
{
    setwd("../Resample/")
    source("resampOct.R")
    refSig = resampOct(x = sig, h = resampFilter(p, q), p = p, q = q)
    plot(t2, refSig, type = 'l', col = "blue4")
    lines(t2, resSig, type = 'l', col = "red4")
}



rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
print(getwd())
library(signal)
library(Rcpp)
source("sos.R")

sourceCpp('sosFilter.cpp')
sourceCpp('sosFilterQ16.cpp')


#--------------- up_sos_dn test

fsamp = 200
fsig = 3.7


p = 2
q = 3

p = 1
q = 5

p = 3
q = 5

t = 1:(fsamp * 1.1) / fsamp
sig = 100 * (sin(t * 2 * pi * fsig) + 0.4 * sin(t * 2 * pi * fsamp * .27))

plot(t, sig, type = 'l', col = "blue")

#dec_tf = cheby1(8, 0.01, 0.8 / q)

# https://en.wikipedia.org/wiki/Upsampling#Interpolation_filter_design

dec_tf = butter(8, 0.8 / max (p, q))
dec_sos = tf2sos(dec_tf)

# wykres tłumienia
if(F)
{ 
    char = freqz(filt = dec_tf, Fs = 16, n = 2 ^ 16)
    char$db = 10 * log10(Mod(char$h))
    #plot(char$f, char$db, type = 'l', ylim = c(-13, 0))
    plot(char$f, 1 / Mod(char$h), type = 'l', ylim = c(0, 100), xlim = c(0, 2))
    
    10^(13/10)
    10 * log10(20)
}

if(F)
{
    gsos = dec_sos$sos
    gsos[1, 1:3] = gsos[1,1:3] * dec_sos$g
    gsos[, 1:3]
    gsos[, 4:6]
}


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
lines(t2, resSig, type = 'l', col = "darkgreen")

gains = rep(dec_sos$g ^ (1 / nrow(dec_sos$sos)), nrow(dec_sos$sos))

# tak jak w pythonie - całe wzmocnienie na pierwszej sekcji
if (F)
{
    gains = c(dec_sos$g, rep(1, nrow(dec_sos$sos) - 1))
    gains[gains < 2^-16] = 2^-16
}

stopifnot(all(2^16 * gains >= 1))

# włączenie gains do sos
for(r in 1:nrow(dec_sos$sos))
{
    dec_sos$sos[r, 1:3] = dec_sos$sos[r, 1:3] * gains[r]
    gains[r] = 1
}
dec_sos$g = 1

resSig16 = UpSosDnQ16_cpp(sig, dec_sos$sos, p, q)
lines(t2, resSig16, type = 'l', col = "red2")



if(F)
{
    setwd("../Resample/")
    source("resampOct.R")
    refSig = resampOct(x = sig, h = resampFilter(p, q), p = p, q = q)
    plot(t2, refSig, type = 'l', col = "blue4")
    lines(t2, resSig, type = 'l', col = "red4")
}


if(F)
{
    p = 2
    q = 3
    dec_tf = butter(8, 0.8 / max (p, q))
    dec_sos = tf2sos(dec_tf)
    sosCoef = as.vector(t(dec_sos$sos))
    cat(sprintf("// resamp filter %s/%s\n", p, q))
    cat(sprintf("#define RESAMP_P (%g)\n", p))
    cat(sprintf("#define RESAMP_Q (%g)\n", q))
    cat(sprintf("#define N_SEC (%g)\n", nrow(dec_sos$sos)))
    cat(sprintf("#define GAIN (%.6g)\n", dec_sos$g))
    cat(sprintf("double sos[] = { %s };\n", paste(sprintf("%.10g", sosCoef), collapse = ", ")))
    
    
}



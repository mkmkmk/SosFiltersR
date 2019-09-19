rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
print(getwd())
library(signal)
library(Rcpp)
source("sos.R")

sourceCpp('sosFilter.cpp')


p = 1

q = 16/5
q = 4
q = 2
order = 8
gainMul = 1

dec_tf = butter(order, 0.8 / max (p, q))
dec_sos = tf2sos(dec_tf)

# korekcja gain-a
dec_sos$g = dec_sos$g * gainMul


if(F)
{
    sosCoef = as.vector(t(dec_sos$sos))
    cat("\n\n//-------\n")
    cat(sprintf("// float resamp filter %d/%d, order %d, gain x%g\n", p, q, order, gainMul))
    cat(sprintf("#define RESAMP_P (%g)\n", p))
    cat(sprintf("#define RESAMP_Q (%g)\n", q))
    cat(sprintf("#define N_SEC (%g)\n", nrow(dec_sos$sos)))
    cat(sprintf("#define GAIN (%.6g)\n", dec_sos$g))
    cat(sprintf("double sos[] = { %s };\n", paste(sprintf("%.10g", sosCoef), collapse = ", ")))
    cat("\n")
    
} else # Q16.16 
{
    # rozdzielenie gains
    gains = rep(dec_sos$g ^ (1 / nrow(dec_sos$sos)), nrow(dec_sos$sos))
    
    # włączenie gains do sos
    for(r in 1:nrow(dec_sos$sos))
    {
        dec_sos$sos[r, 1:3] = dec_sos$sos[r, 1:3] * gains[r]
        gains[r] = 1
    }
    dec_sos$g = 1

    sosCoef = as.vector(t(dec_sos$sos))
    sosCoef = round(sosCoef * 2^16)
    
    br = range(log2(abs(sosCoef)))
    cat(sprintf("bit range: [%.4g, %.4g]\n", br[1], br[2]))
    
    cat("\n\n//-------\n")
    cat(sprintf("// Q.16 resamp filter %d/%d, order %d, gain x%g\n", p, q, order, gainMul))
    cat(sprintf("#define RESAMP_P (%g)\n", p))
    cat(sprintf("#define RESAMP_Q (%g)\n", q))
    cat(sprintf("#define SOS_N_SEC (%g)\n", nrow(dec_sos$sos)))
    cat(sprintf("int32_t SOS[] = { %s };\n", paste(sprintf("%.10g", sosCoef), collapse = ", ")))
    cat("\n")
    
}




rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
print(getwd())
library(signal)
library(Rcpp)
source("sos.R")

sourceCpp('sosFilter.cpp')

#--------------- decym test
fsamp = 1000
fsig = 3.7

q = 20
t = 1:(fsamp * 1) / fsamp
sig = sin(t * 2 * pi * fsig) + 0.4 * sin(t * 2 * pi * fsamp * .27)

plot(t, sig, type = 'l', col = "blue")

#dec_tf = cheby1(8, 0.01, 0.8 / q)
dec_tf = butter(8, 0.8 / q)
dec_sos = tf2sos(dec_tf)

if(F)
    paste(as.vector(t(dec_sos$sos)), collapse = ", ")

dec_zi = sos_zi(dec_sos) * mean(sig)

decSig = sig
#decSig = sosFilter(dec_sos, sig, zi = dec_zi)
decSig = sosFilter_cpp(dec_sos$sos, dec_sos$g, sig, zi = dec_zi)
decSig = decSig[seq(1, length(sig), q)]
t2 = t[seq(1, length(sig), q)]
lines(t2, decSig, type = 'l', col = "red")


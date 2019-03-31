rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
print(getwd())
library(signal)
library(Rcpp)
source("sos.R")

sourceCpp('sosFilter.cpp')


p = 2
q = 3
dec_tf = butter(8, 0.8 / max (p, q))
dec_sos = tf2sos(dec_tf)
sosCoef = as.vector(t(dec_sos$sos))
cat("\n\n//-------\n")
cat(sprintf("// resamp filter %s/%s\n", p, q))
cat(sprintf("#define RESAMP_P (%g)\n", p))
cat(sprintf("#define RESAMP_Q (%g)\n", q))
cat(sprintf("#define N_SEC (%g)\n", nrow(dec_sos$sos)))
cat(sprintf("#define GAIN (%.6g)\n", dec_sos$g))
cat(sprintf("double sos[] = { %s };\n", paste(sprintf("%.10g", sosCoef), collapse = ", ")))
cat("\n")


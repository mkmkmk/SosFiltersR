rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
print(getwd())
library(signal)
library(Rcpp)
source("sos.R")

sourceCpp('sosFilter.cpp')


tf = butter(4 / 2, c(1.5, 23) / (1000 / 2))
z = roots(tf$b)
p = roots(tf$a)
k = tf$b[1] / tf$a[1]

zp2sos(z,p,k)
tf2sos(tf)


tf = butter(4, c(0.2, 0.6))
tf2sos(tf)
# [b, a] = butter(4, [0.2, 0.6]);
# [z,p,k] = tf2zp(b, a)
# [sos_var,g] = zp2sos(z, p, k)
# sos_var =
#     1.0000000   2.0001266   1.0001266   1.0000000   0.4575326   0.6245384
#     1.0000000   1.9998734   0.9998734   1.0000000  -0.0073608   0.1956041
#     1.0000000  -2.0000000   1.0000000   1.0000000  -0.8822287   0.3301957
#     1.0000000  -2.0000000   1.0000000   1.0000000  -1.3945409   0.7466724
# 
# g =  0.046583
#
# [z,p,g]=butter(4, [0.2, 0.6]);
# [sos_var,g] = zp2sos(z, p, k)
# sos_var =
#     1.0000000   2.0000000   1.0000000   1.0000000   0.4575326   0.6245384
#     1.0000000   2.0000000   1.0000000   1.0000000  -0.0073608   0.1956041
#     1.0000000  -2.0000000   1.0000000   1.0000000  -0.8822287   0.3301957
#     1.0000000  -2.0000000   1.0000000   1.0000000  -1.3945409   0.7466724
# g =  0.0041592




fs = 1000
tf2sos(butter(4 / 2, c(1.8, 30) / (fs/2)))
# [z,p,g]=butter(4/2, [1.8, 30]/(fs/2));
# [sos_var,g] = zp2sos(z, p, k)
# sos_var =
#     1.00000   2.00000   1.00000   1.00000  -1.76226   0.79070
#     1.00000  -2.00000   1.00000   1.00000  -1.98425   0.98439
# g =  0.0041592


tf = butter(4, 0.6, type = "low")
tf2sos(tf)
# [z,p,g]=butter(4, 0.6);
# [sos_var,g] = zp2sos(z, p, k)
# sos_var =
#     1.000000   2.000000   1.000000   1.000000   0.453120   0.466326
#     1.000000   2.000000   1.000000   1.000000   0.328976   0.064588
# g =  0.046583


tf = butter(4, 10 / (fs/2), type = "low")
tf2sos(tf)
# [z,p,g]=butter(4,  10 / (fs/2), type = "low");
# [sos_var,g] = zp2sos(z, p, k)

# { {1,2,1}, {1,-1.8866095826215084000,0.89033973628402607000}, 8.9848614639706479000e-007,  1  },
# { {1,2,1}, {1,-1.9492159580258410000,0.95306989532789022000}, 1, 1 }  }




#---------------------
#analiza 1 sekcji
tf = butter(2, 1 / (fs/2), type = "high")
x = 1000 + 1:5000 %% 100
plot(filter(tf, x, init.x = rep(0, 2)), type = 'l', col = 'red')
lines(filter(tf, x, init.x = rep(x[1], 2)), type = 'l', col = 'darkgray')


#flt2 = filter(tf, c(rep(x[1], 2), x))

flt2 = hackFilter(tf$b, tf$a, x, rep(x[1], 2))
lines(flt2, type = 'l', col = 'green')

#debugonce("sosFilter")
prepSos = list(g = 1, sos = t(c(tf$b, tf$a)))
# flt3 = sosFilter(prepSos, x, states = t(c(26704772, 26704772, 26703997)))

flt3 = sosFilter(prepSos, x)
lines(flt3, type = 'l', col = 'black')

flt4 = sosFilter(prepSos, x, zi = x[1] *biquad_zi(tf$b, tf$a))
lines(flt4, type = 'l', col = 'black')




# ---------------
# 2 sekcje

x = 1000 + 1:5000 %% 100

tf = butter(2, c(1,8) / (fs/2), type = "pass")
sos = tf2sos(tf)


plot(filter(tf, x, init.x = rep(0, 4)), type = 'l', col = 'red')
lines(filter(tf, x, init.x = rep(x[1], 4)), type = 'l', col = 'darkgray')

#plot(filter(tf, x, init.x = rep(x[1], 4)), type = 'l', col = 'red')


ft_ref = filter(tf, x, init.x = rep(x[1], 4))
ft_ver = sosFilter(sos, x)

#plot(ft_ref, type='l', col = 'blue', ylim = range(ft_ref))
lines(ft_ver, type = 'l', col = 'red')

# debugonce("sosFilter")


ft_ver2 = sosFilter(sos, x, zi = sos_zi(sos) * x[1])
lines(ft_ver2, type='l', col = 'blue')

plot(ft_ref, type='l', col = 'blue')
lines(ft_ver2, type='l', col = 'red')

ft_ver3 = sosFilter_cpp(sos$sos, sos$g, x, zi = sos_zi(sos) * x[1])
plot(ft_ver3, type='l', col = 'red4')

maxDiff = max(abs(ft_ref - ft_ver3))
maxDiff
stopifnot(maxDiff < 1e-5)


#---------------
# OldButt4_20_300_1000Hz

prepSos = 
    list(
        g = 0.387132401864575, 
        sos = t(matrix(c(1, 2, 1, 1, 0.365079803471447, 0.197369861882712,1, -2, 1, 1, -1.98222897911192, 0.982386928837799), ncol = 2)))

prepSos$sos

zi = sos_zi(prepSos) 

x = 1564 + 1:5000 %% 100
x = x*1e-9

ft = sosFilter(prepSos, x)
ft_zi = sosFilter(prepSos, x, zi*x[1])


plot(ft, type = 'l', col = 'red')
lines(ft_zi, type = 'l', col = 'blue')

paste(zi[1,], collapse = ', ')
paste(zi[2,], collapse = ', ')

zi*prepSos$g

c(	9.448075803990889E-07,
	2.9960740126759604E-07,
	-1.550421448895139E-06,
	1.5504214488951396E-06) / 1.56435851E-06

# wersja w C# może mieć wzmocnienie na wejściu dlatego ew. trzeba zi pomnożyć przez g
sos = tf2sos(butter(2, c(0.1, 2.1) / (1000 / 2), type = "pass"))
sos_zi(sos) * sos$g

sos = tf2sos(butter(2, c(2, 300) / (1000 / 2), type = "pass"))
sos_zi(sos) * sos$g

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
print(getwd())
library(signal)
source("sos.R")


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




sectionFormII = function(state, section, value)
{
    b = section[1:3]
    a = section[4:6]
    state[1] = value - a[2] * state[2] - a[3] * state[3]
    res = sum(state * b)
    state[2:3] = state[1:2]
    list(res = res, state = state)
}


sectionFormIItrans = function(state, section, value)
{
    b = section[1:3]
    a = section[4:6]

    res =      b[1] * value + state[2];
    state[2] = b[2] * value - a[2] * res + state[3];
    state[3] = b[3] * value - a[3] * res;

    list(res = res, state = state)
}
# debugonce(sectionFormIItrans)



sosFilter = function(sos, x, zi)
{
   nsec = nrow(sos$sos)
   
   state = matrix(rep(0, 3 * nsec), nrow = nsec)
   
   if(!missing(zi))
   {
       state[, 2:3] = zi
   }

   y = x * 0
   for(i in 1:length(x))
   {          
       value = x[i]
       for (r in 1:nsec)
       {
           cmp = sectionFormIItrans(state = state[r,], section = sos$sos[r,], value)
           value = cmp$res
           state[r,] = cmp$state
       }
       y[i] = value
   }
   #print(state)
   y * sos$g
}


# wg. scipy.signal.lfilter_zi(b, a)
biquad_zi = function(b, a)
{
    companT = matrix(c(-a[2:3], c(1, 0)), nrow = 2)
    ImA = diag(2) - companT
    B = b[2:3] - a[2:3] * b[1]
    solve(ImA, B)
}



sos_zi = function(sos)
{
    nsec = nrow(sos$sos)
    zi = matrix(rep(0, 2 * nsec), nrow = nsec)

    scale = 1
    
    for (r in 1:nsec)
    {
        section = sos$sos[r, ]
        b = section[1:3]
        a = section[4:6]
        zi[r, ] = scale * biquad_zi(b, a)
        scale = scale * sum(b) / sum(a)
    }
    zi
}

hackFilter = function(filt, a, x, init.x=c())
{
    x <- stats::filter(c(init.x, x), filt / a[1], sides = 1)
    x <- na.omit(x, filt / a[1] , sides = 1)
    x <- stats::filter(x, -a[-1] / a[1], method = "recursive")
    x
}



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













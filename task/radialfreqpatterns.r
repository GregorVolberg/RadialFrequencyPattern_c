r0 = 1
A  = 0.1
phi = 0
frq = 10
theta = seq(-pi, pi, by = pi/180)
R = r0 * (1+A*sin(frq*theta+phi))
x = R*cos(theta)
y = R*sin(theta)
plot(x, y, type = 'l', xlim = c(-2, 2), ylim = c(-2, 2))
plot(R, type = 'l')


# mit spikes, siehe DOI: 10.1038/srep26681
# oder auch https://en.wikipedia.org/wiki/Triangle_wave, section Harmonics
# phi nicht setzbar
r0 = 1
A  = 0.1
frq = 10
inc = pi/360
theta = seq(-pi, pi, by = inc)
nodd = 5 #(number of odd harmonics) (converges quickly)

erg = numeric(length(theta))
  for (k in 0:nodd){
  tmp = (-1)^k * (sin((2*k+1)*frq*theta) / (2*k+1)^2)
  erg = erg + tmp}

R = r0 * (1+A*erg)

x = R*cos(theta)
y = R*sin(theta)

plot(x, y, type = 'l', xlim = c(-2, 2), ylim = c(-2, 2))
#plot(R, type='l')

# try a different base figure (other than circle)

r0 = 1
A  = 0.1
phi = 0
frq = 3
inc = pi/360
theta = seq(-pi, pi, by = inc)
Rbase = r0 * (1+A*sin(frq*theta+phi))
#x = Rbase*cos(theta)
#y = Rbase*sin(theta)
#plot(x, y, type = 'l', xlim = c(-2, 2), ylim = c(-2, 2))


Rn = Rbase*(1+A*erg)
x  = Rn*cos(theta)
y  = Rn*sin(theta)
plot(x, y, type = 'l', lwd = 2, xlim = c(-2, 2), ylim = c(-2, 2))

# try window function for triang/ circ
As   = c(0.2, 0.5, 1.3)
frqs = c(2, 4, 8)
phis = c(pi/2, pi, 0)
inc = pi/360
theta = seq(-pi, pi, by = inc)

erg2 = numeric(length(theta))
for (j in 1:3){
  tmp = As[j]*sin(frqs[j]*theta+phis[j])
  erg2 = erg2+tmp
    }
eerg = erg*erg2

Rn = Rbase*(1+A*eerg)
x  = Rn*cos(theta)
y  = Rn*sin(theta)
plot(x, y, type = 'l', lwd = 2, xlim = c(-2, 2), ylim = c(-2, 2))

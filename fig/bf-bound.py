#!/usr/bin/python2.7

from pylab import *

# create an array of 256 points (x), from -pi to +pi
points = 100

q = 1000
m = 1<<12
n = 100
k = 28
r = 1
salt = 128
M = np.linspace(1, m, points)

def Mitz(k, n, m):
    return (1 - np.exp(-1 * float(k*n) / float(m)))**k

Ym = np.cos(M)

# compute the plot
plot(M, map(lambda X:Mitz(k,n,X), M))

# Show it in a new window
savefig('bf-bound.pdf')

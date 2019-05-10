#!/usr/bin/python2.7

from pylab import *

# create an array of 256 points (x), from -pi to +pi
points = 100

def Mitz(k, n, m):
    return (1 - np.exp(-1 * float(k*n) / float(m)))**k

def Us(r, q, k, n, m):
    p = float(Mitz(k, n, m))
    x = p * q
    if r <= x:
        return 1.0
    return ((x / r)**r) * np.exp(r - x)

q = 1<<30
n = 100
k = 16
r = 1

start = 1<<12
end = (1<<16)
M = np.linspace(start, end, (end-start)/8)

N = range(100,1100,100)

# compute the plot
lines = []
for n in N:
    line, = plot(M, map(lambda X:Us(r,q,k,n,X), M))
    lines.append(line)
legend([lines[0], lines[-1]], ['$n=%d$' % N[0], '$n=%d$' % N[-1]])

# Show it in a new window
savefig('bf-bound.pdf')

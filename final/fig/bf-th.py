#!/usr/bin/python2.7

from pylab import *

# Fonts for the plot
rc('font', **{
        'family' : 'serif',
        #'serif': ['litertime']
        'size' : 7
    }
)
rc('text', usetex=True)

# How many bits in a KB
kb = 8192.0
mb = 1024.0*kb
gb = 1024.0*mb

# Kirsch-Mitzenmacher
def Mitz(k, n, m):
    return (1 - np.exp(-1 * float(k*n) / float(m)))**k


def Thresh(q, r, k, l, m):
    p = (float(l+k)/float(m))**k
    x = p * q
    if r <= x:
        return 1.0
    return ((x / r)**r) * np.exp(r - x)

# Our \Zeta function
def Us(q, r, k, n, m):
    p = float(Mitz(k, n, m))
    x = p * q
    if r <= x:
        return 1.0
    return ((x / r)**r) * np.exp(r - x)



# Fixed parameters
qorder = 32
q = 1<<qorder
n = 10**9
k = 16
r = 1

# Range of filter length
start = 7*gb
end   = 15*gb
M = np.linspace(start, end, 10)

# Which parameter do you want to plot on the x-axis?
N = range(10**8,10**9+10**8,10**8)
R = range(1,11,1)
param = R

figure(figsize=(3.6,1.5))
yscale('log', basey=2)

r = 1
k = 16
l = (n * k)
th, = plot(M/gb, map(lambda X:Thresh(q,r,k,l,X), M), ':', lw=0.8)
bf, = plot(M/gb, map(lambda X:Us(q,r,k,n,X), M), lw=0.8)
print Thresh(q,r,k,l,start)
legend([bf, th], ['$n$-capped', '$\ell$-thresholded'])
xticks(np.arange(0,int(end/gb+1),1))
xlabel('Filter length (GB)')
xlim(start/gb, end/gb)
ylim(10**-6, 1)
grid(lw=0.25)
tight_layout()

# Threshold

# Show it in a new window
savefig('bf-th.pdf')

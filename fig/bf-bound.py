#!/usr/bin/python2.7

from pylab import *

rc('font', **{
        'family' : 'serif',
        'serif': ['litertime']
        #'size' : 7
    }
)
rc('text', usetex=True)

def Mitz(k, n, m):
    return (1 - np.exp(-1 * float(k*n) / float(m)))**k

def Us(q, r, k, n, m):
    p = float(Mitz(k, n, m))
    x = p * q
    if r <= x:
        return 1.0
    return ((x / r)**r) * np.exp(r - x)

qorder = 64
q = 1<<qorder
n = 100
k = 16
r = 1

kb = 8192.0
start = 19500 / kb
end   = (start + 40000) / kb
M = np.linspace(start*kb, end*kb, (end*kb-start*kb)/10)

N = range(100,1100,100)
R = range(1,42,3)

# compute the plot
param = R

figure(figsize=(6,2.25))
lines = []
yscale('log', basey=2)
for (i, r) in enumerate(param):
    a = 0.3 + (0.7* (i+1) / len(param))
    line, = plot(M/kb, map(lambda X:Us(q,r,k,n,X), M),
                 color=(65/256.0,115/256.0,150/256.0), alpha=a, lw=1)
    lines.append(line)
legend([lines[0], lines[-1]], ['$r=%d$' % param[0], '$r=%d$' % param[-1]])
xticks(np.arange(0,int(end+1),0.5))
#ylabel('$\zeta(q,r,k,m,n)$')
xlabel('Filter length (KB)')
title('$\zeta$ value for $q=2^{%d}$, $k=%d$, $n=%d$ and varying $r$.' % (
        qorder, k, n))
xlim(start, end)
ylim(10**-6, 1)
tight_layout()

# Show it in a new window
savefig('bf-bound.pdf')

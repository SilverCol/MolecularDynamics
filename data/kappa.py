from matplotlib import pyplot as plt
import numpy as np
from scipy import optimize as opt

k1 = []
dk1 = []
k2 = []
dk2 = []
N = []

for i in range(300):
    try:
        data1 = np.fromfile('fluxFrame' + str(i) + 'a.bin')
        data2 = np.fromfile('mFluxFrame' + str(i) + 'a.bin')
    except:
        continue
    N.append(1/i)
    k1.append(np.mean(data1))
    dk1.append(np.std(data1))
    k2.append(np.mean(data2))
    dk2.append(np.std(data2))

plt.rcParams.update({'font.size': 20})
fig = plt.figure()
yMax = .05
yMin = -.01
xMin = 0
xMax = 1.1*max(N)
ax = fig.subplots(subplot_kw=dict(aspect='auto', autoscale_on=False, xlim=(xMin, xMax), ylim=(yMin, yMax)))
ax.grid()


def linear(x, k):
    return k*x


popt1, pcov1 = opt.curve_fit(linear, N, k1, sigma=dk1, absolute_sigma=True)
print(popt1)
print(pcov1)
popt2, pcov2 = opt.curve_fit(linear, N, k2, sigma=dk2, absolute_sigma=True)
print(popt2)
print(pcov2)
xes = np.linspace(xMin, xMax)

ax.errorbar(N, k1, dk1, fmt='.', barsabove=True, color='blue')
line1, = ax.plot(xes, linear(xes, popt1[0]), color='blue')
ax.errorbar(N, k2, dk2, fmt='.', barsabove=True, color='orange')
line2, = ax.plot(xes, linear(xes, popt2[0]), color='orange')

ax.legend((line1, line2), ('Nose-Hoover: $\\kappa = %.2f \\pm %.2f$' % (popt1[0], pcov1[0][0]),
                           'Maxwell: $\\kappa = %.2f \\pm %.2f$' % (popt2[0], pcov2[0][0])), loc=2)
ax.set_ylabel('$\\langle J_j \\rangle$')
ax.set_xlabel('$1/N$')

plt.show()

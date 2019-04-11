from matplotlib import pyplot as plt
import numpy as np
import matplotlib.animation as animation

file = 'contact'
txt = input('Number? ')
file += txt
txt = input('Anharmonic? ')
if txt == 'y':
    file += 'a'
file += '.bin'

data = np.fromfile(file)

t = []
j1 = []
j2 = []
for n, entry in enumerate(data):
    stage = n % 3
    if stage == 0:
        t.append(entry)
    elif stage == 1:
        j1.append(entry)
    elif stage == 2:
        j2.append(entry)

plt.rcParams.update({'font.size': 20})
fig = plt.figure()
yMax = .1
yMin = -.05
xMin = min(t)
xMax = max(t)
ax = fig.subplots(subplot_kw=dict(aspect='auto', autoscale_on=False, xlim=(xMin, xMax), ylim=(yMin, yMax)))
ax.grid()

line1, = ax.plot(t, j1, '-')
line2, = ax.plot(t, j2, '-')
ax.set_xscale('log')
ax.set_ylabel('$J$')
ax.set_xlabel('$\\tau$')
ax.legend((line1, line2), ('$J_1$', '$J_N$'))

plt.show()

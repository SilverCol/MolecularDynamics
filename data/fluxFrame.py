from matplotlib import pyplot as plt
import numpy as np
import matplotlib.animation as animation

txt = input('Maxwell? ')
if txt == 'y':
    file = 'mFluxFrame'
else:
    file = 'fluxFrame'
txt = input('Number? ')
file += txt
txt = input('Anharmonic? ')
if txt == 'y':
    file += 'a'
file += '.bin'

data = np.fromfile(file)

fig = plt.figure()
yMax = 1.1*max(data)
yMin = 0
xMin = 0
xMax = len(data) - 1
ax = fig.subplots(subplot_kw=dict(aspect='auto', autoscale_on=False, xlim=(xMin, xMax), ylim=(yMin, yMax)))
ax.grid()

m = np.mean(data)
m = [m for val in data]

ax.plot(data, '-o')
ax.plot(m, '-.')

plt.show()

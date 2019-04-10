from matplotlib import pyplot as plt
import numpy as np
import matplotlib.animation as animation

txt = input('Maxwell? ')
if txt == 'y':
    file = 'mprofile'
else:
    file = 'profile'
txt = input('Number? ')
file += txt
txt = input('Anharmonic? ')
if txt == 'y':
    file += 'a'
file += '.bin'

data = np.fromfile(file)

fig = plt.figure()
yMax = 2.1
yMin = .9
xMin = 0
xMax = len(data) - 1
ax = fig.subplots(subplot_kw=dict(aspect='auto', autoscale_on=False, xlim=(xMin, xMax), ylim=(yMin, yMax)))
ax.grid()


ax.plot(data, '-o')

if txt != 'y':
    m = np.mean(data)
    m = [m for val in data]
    ax.plot(m, '-.')

plt.show()

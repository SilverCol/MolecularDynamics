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

plt.rcParams.update({'font.size': 20})
fig = plt.figure()
yMax = 2.1
yMin = .9
xMin = 0
xMax = len(data) - 1
ax = fig.subplots(subplot_kw=dict(aspect='auto', autoscale_on=False, xlim=(xMin, xMax), ylim=(yMin, yMax)))
ax.grid()


ax.plot(data, '-o')

if txt != 'y':
    mn = np.mean(data[1:len(data)-1])
    dev = np.std(data[1:len(data)-1])
    m = [mn for val in data]
    line, = ax.plot(m, '-.')
    ax.legend((line,), ('$%.2f \\pm %.2f$' % (mn, dev),))

ax.set_ylabel('$\\langle p_j \\rangle$')
ax.set_xlabel('$j$')

plt.show()

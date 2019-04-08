from matplotlib import pyplot as plt
import numpy as np
import matplotlib.animation as animation


def myMax(array2d):
    c = -10e100
    for array in array2d:
        if max(array) > c:
            c = max(array)
    return c


def myMin(array2d):
    c = 10e100
    for array in array2d:
        if min(array) < c:
            c = min(array)
    return c


txt = input('Maxwell? ')
if txt == 'y':
    file = 'mTProfile'
else:
    file = 'tProfile'
txt = input('Number? ')
file += txt
txt = input('Anharmonic? ')
if txt == 'y':
    file += 'a'
file += '.bin'

data = np.fromfile(file)
binDelimiter = -1234567891

T = []
t = []
time = []

extractTime = True
for entry in data:
    if extractTime:
        time.append(entry)
        extractTime = False
        continue
    elif int(entry) == binDelimiter:
        T.append(t)
        t = []
        extractTime = True
        continue
    t.append(entry)

fig = plt.figure()
yMax = myMax(T)
yMin = myMin(T)
xMin = 0
xMax = len(T[0]) - 1
ax = fig.subplots(subplot_kw=dict(aspect='auto', autoscale_on=False, xlim=(xMin, xMax), ylim=(yMin, yMax)))
ax.grid()

image = ax.plot([], [], '.')[0]
text = ax.text(.02, .95, '', transform=ax.transAxes)

indices = []
for n, t in enumerate(T[0]):
    indices.append(n);


def animate(i):
    image.set_data(indices, T[i])
    text.set_text('$time$ = %.2f' % time[i])
    return image, text,


ani = animation.FuncAnimation(fig, animate, frames=len(T), interval=40, blit=True)

# ani.save('../report/RENAME.gif', writer='imagemagick', fps=25)
plt.show()

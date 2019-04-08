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
xMin = min(time)
xMax = max(time)
ax = fig.subplots(subplot_kw=dict(aspect='auto', autoscale_on=False, xlim=(xMin, xMax), ylim=(yMin, yMax)))
ax.grid()

T = np.transpose(T)
Tm = np.mean(T, axis=0)

lines = []
for t in T:
    line = ax.plot(time, t, '-.', lw=.5)[0]
    lines.append(line)
lines.append(ax.plot(time, Tm, '-k')[0])

indices = []
for n, t in enumerate(T[0]):
    indices.append(n);

plt.show()

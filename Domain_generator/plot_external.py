import numpy as np
import matplotlib.pyplot as plt

with open('../Domain/data/arctic_crude_external.txt', 'rt') as f:
    next(f)
    data = f.readlines()

fig = plt.figure(figsize=(20, 20))
xs = np.empty((0, 0))
ys = np.empty((0, 0))
counter = 0
colors = [
    'tab:blue',
    'tab:orange',
    'tab:green',
    'tab:red',
    'tab:purple',
    'tab:brown',
    'tab:pink',
    'tab:gray',
    'tab:olive',
    'tab:cyan'
]

style = ['-', '--', ':']

for i, line in enumerate(data):
    if line[0] == '#':
        plt.plot(xs, ys, linestyle=style[counter%3], c=colors[counter%10], label=str(counter))
        counter += 1
        xs = np.empty((0, 0))
        ys = np.empty((0, 0))
    else:
        x, y = [float(x) for x in line.split(' ')]
        xs = np.append(xs, x)
        ys = np.append(ys, y)
plt.plot(xs, ys, linestyle=style[counter%3], c=colors[counter%10], label=str(counter))
plt.grid(True)
plt.legend(loc='best')
fig.savefig('../Domain/pictures/external.png')

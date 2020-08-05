import numpy as np
import matplotlib.pyplot as plt
import matplotlib._color_data as mcd

colors_num = 21
with open('../../Domain/data/arctic_crude_coarsened_external.txt', 'rt') as f_ext:
    next(f_ext)
    data_ext = f_ext.readlines()

with open('../../Domain/data/arctic_crude_coarsened_islands.txt', 'rt') as f_isl:
    data_isl = f_isl.readlines()

data = data_ext + data_isl
fig = plt.figure(figsize=(20, 20))
xs = np.empty((0, 0))
ys = np.empty((0, 0))
counter = 0
colors = [name for name in mcd.CSS4_COLORS
           if "xkcd:" + name in mcd.XKCD_COLORS][4:5+colors_num]
style = ['-', '--', '-.']

for i, line in enumerate(data):
    if line[0] == '#':
        plt.plot(xs, ys,
		 linestyle=style[counter%len(style)],
		 c=colors[counter%len(colors)],
		 label=str(counter))
        counter += 1
        xs = np.empty((0, 0))
        ys = np.empty((0, 0))
    else:
        x, y = [float(x) for x in line.split(' ')]
        xs = np.append(xs, x)
        ys = np.append(ys, y)
plt.plot(xs, ys,
	 linestyle=style[counter%len(style)],
	 c=colors[counter%len(colors)],
	 label=str(counter))
plt.grid(True)
plt.legend(loc='best')

fig.savefig('../../Domain/pictures/full_coarsened.png')

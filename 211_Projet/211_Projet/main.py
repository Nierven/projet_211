from fem_resolution import *
from source import *

from time import time
from math import *

from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import lines
import numpy as np

# Variables & Constants
# =====================

# map(width, dx, k, m)
m = map(10, 0.2, 0.1, 0.01)
sources = []

mult = 1.0
map_range = np.arange(0, 10*mult)

start_time = 2
stop_time = 30
time_speed = 1       # Time speed multiplier
fps = 10             # Frames per second
dt = time_speed/fps

times = np.arange(0, stop_time, dt)

# Model & Resolution
# ==================

def acoustic_level(p):
    sum = 0
    for src in sources:
        # sum distances relative to speaker position
        sum += src.level / (1 + src.p.distance(p))**2
    return sum

#resolve_acoustics_problem()

# Problem optimization
# ====================

sources.append(source(point(3, 5), 90))
sources.append(source(point(7, 5), 50))

positions_x = [src.p.x * mult for src in sources]
positions_y = [src.p.y * mult for src in sources]

m.sources.append(m.boules[int(len(m.boules)/2)][int(len(m.boules)/2)])

# Solutions display
# =================

# Create a new plot
fig_plan = plt.figure(figsize=[8, 8])
ax = fig_plan.gca()

#ax.grid(True)
ax.set_xlim(-2*mult, 12*mult)
ax.set_ylim(-2*mult, 12*mult)

#c = ax.pcolor([[acoustic_level(point(x/mult, y/mult)) for x in map_range] for y in map_range], vmin=0, vmax=100, cmap='bwr')
#ax.scatter(positions_x, positions_y, s=20)
#fig_plan.colorbar(c, ax=ax)

# Draw the map
#for seg in segments:
#    ax.add_line(lines.Line2D([seg[0][0], seg[1][0]], [seg[0][1], seg[1][1]]))

boule_scatter = ax.scatter([], [], s=15)
source_scatter = ax.scatter([], [], s=20, c='r')
b_src = m.sources[0]
#b_src.p.x += 0.0005
f = 1
a = 0.5

global avg_ctime
global cpt

avg_ctime = 0
cpt = 0

def animate(t):
    global avg_ctime
    global cpt

    # 2D animation
    if t >= start_time:
        # Get next source state
        if t < start_time + 0.5/f:
            b_src.p.x += a*sin((t-start_time) * 2*pi*f)
        else:
            m.sources.clear()
        #b_src.p.y += a*sin(t * 2*pi*f)

        # Solve step time
        tstart = time()
        m.resolution(dt)
        tstop = time()

        avg_ctime += tstop - tstart
        cpt += 1
    
    boule_scatter.set_offsets(np.c_[[b.p.x for (i, j), b in np.ndenumerate(m.boules)], [b.p.y for (i, j), b in np.ndenumerate(m.boules)]])
    source_scatter.set_offsets(np.c_[[b.p.x for b in m.sources], [b.p.y for b in m.sources]])

    #return everything that must be updated
    return (boule_scatter, source_scatter)

anim = animation.FuncAnimation(fig_plan, animate, frames=times, interval = int(1000/fps), blit=False)

## Show the plot to the screen
plt.show()

print('average compute time: {} s'.format(avg_ctime / cpt))
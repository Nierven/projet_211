from sound import Map, Particle, Source
from utilities import Point, print_header, print_progress_bar

from time import time
from math import pi, sin
from copy import deepcopy

from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

verbose = True
deferred_mode = True

# Variables & Constants
# =====================

# Geometry properties
room_width = 10            # m

# Physical properties
dx = 0.2                   # m
rho = 1.2e-3               # kg/m**3
k = 0.0004                 # mN/m
mass = dx**3 * rho         # ng

# ACCELERATION = m/s**2 * 1e-12

# Map properties
m = Map(room_width, dx, k, mass)
color_mult = 1.0
colormap_range = np.arange(0, room_width*color_mult)

# Time properties
start_time = 1             # Time to start simulation [s]
stop_time = 10              # Time to stop simulation  [s]
time_speed = 1             # Time speed multiplier    [1]
fps = 30                   # Frames per second        [frame/s]
dt = time_speed/fps        # Time interval            [s]

times = np.arange(0, stop_time, dt)

# Model & Resolution
# ==================

watching_particle = m.particles[len(m.particles) // 2][len(m.particles) - 5]
def compute_velocity(dt):
    """Compute velocity of the wave"""
    if 't' not in compute_velocity.__dict__: compute_velocity.t = 0
    if 'done' not in compute_velocity.__dict__: compute_velocity.done = False

    compute_velocity.t += dt; tolerance = 1e-10
    if watching_particle.p.distance(watching_particle.p_0) >= tolerance and not compute_velocity.done:
        speed = watching_particle.p_0.distance(b_src.particle.p_0) / compute_velocity.t
        #print("Wave velocity: {}".format(speed)) # TODO afficher à la fin de simu
        compute_velocity.done = True

# Problem optimization
# ====================

m.add_source(room_width/2, room_width/2, 90)
b_src = m.sources[0]

f = 1
a = dx*4

# Solutions display
# =================

# Create a new figure
fig_plan = plt.figure(figsize=[8, 8])
ax = fig_plan.gca()

# Set figure limits and properties
ax.set_xlim(-0.2*room_width*color_mult, 1.2*room_width*color_mult)
ax.set_ylim(-0.2*room_width*color_mult, 1.2*room_width*color_mult)

# Create scatter plots
particles_scatter = ax.scatter([], [], s=10, marker='s')
source_scatter = ax.scatter([], [], s=15, marker='s', c='r')

#c = ax.pcolor([[acoustic_level(point(x/mult, y/mult)) for x in map_range] for y in map_range], vmin=0, vmax=100, cmap='bwr')
#ax.scatter(positions_x, positions_y, s=20)
#fig_plan.colorbar(c, ax=ax)
    
global avg_ctime; avg_ctime = 0
global cpt; cpt = 0
particles_history = []
sources_history = []

def compute_all():
    print_progress_bar(0, times[-1], prefix = 'Compute progress:', suffix = 'Complete', length = 50)
    for t in times:
        if t >= start_time:
            compute_step(t)
            print_progress_bar(t, times[-1], prefix = 'Compute progress:', suffix = 'Complete', length = 50)

def compute_step(t):
    global avg_ctime
    global cpt

    # Get next source state
    if t < start_time + 0.5/f:
        b_src.particle.p.x = b_src.particle.p_0.x + a/2*(1 + sin((t-start_time) * 2*pi*f))
    else:
        compute_velocity(dt)
        m.clear_sources()

    # Solve step time
    tstart = time()
    m.resolution(dt)
    tstop = time()

    # Add frame to history
    if deferred_mode:
        particles_history.append(deepcopy(m.particles))
        sources_history.append(deepcopy(m.sources))
        
    # Update average time
    avg_ctime += tstop - tstart
    cpt += 1

# Figure animation live
def animate_live(t):
    global avg_ctime
    global cpt
    
    if t >= start_time:
        compute_step(t)

    # Update all points positions
    particles_scatter.set_offsets(np.c_[[b.p.x for (i, j), b in np.ndenumerate(m.particles)], [b.p.y for (i, j), b in np.ndenumerate(m.particles)]])
    source_scatter.set_offsets(np.c_[[src.particle.p.x for src in m.sources], [src.particle.p.y for src in m.sources]])

    # Return everything that must be updated
    return (particles_scatter, source_scatter)

# Figure animation deferred
def animate_deferred(t):
    global avg_ctime
    global cpt

    # Update all points positions
    particles_scatter.set_offsets(np.c_[[b.p.x for (i, j), b in np.ndenumerate(particles_history[t])], [b.p.y for (i, j), b in np.ndenumerate(particles_history[t])]])
    source_scatter.set_offsets(np.c_[[src.particle.p.x for src in sources_history[t]], [src.particle.p.y for src in sources_history[t]]])

    # Return everything that must be updated
    return (particles_scatter, source_scatter)

if deferred_mode:
    compute_all()
    # Configure animation with blit and repeating
    anim = FuncAnimation(fig_plan, animate_deferred, frames=range(len(particles_history)), interval = 1000//fps, blit=True, repeat=True, repeat_delay=1000)
else:
    # Configure animation with blit and not repeating
    anim = FuncAnimation(fig_plan, animate_live, frames=times, interval = 1000//fps, blit=True, repeat=False)

# Show the plot to the screen
plt.show()

# Show debug outputs
if verbose:
    avg_ctime /= cpt
    print('Average step compute time: {} s'.format(avg_ctime))
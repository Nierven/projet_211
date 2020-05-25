from sound import Map, Particle, Source
from utilities import Point, print_header, print_progress_bar

from time import time
from math import pi, sin
from copy import deepcopy

from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

display = True
verbose = True
deferred_mode = True
gradient = False

# Variables & Constants
# =====================

# Geometry properties
room_width = 10            # m

# Physical properties
dx = 0.25                  # m
rho = 1.2e-3               # kg/m**3
k = 0.05                   # mN/m
mass = dx**3 * rho         # ng

# ACCELERATION = m/s**2 * 1e-12

# Map properties
global m
m = Map(room_width, dx, k, mass)
color_mult = 1.0
colormap_range = np.arange(0, room_width*color_mult)

# Time properties
start_time = 1             # Time to start simulation [s]
stop_time = 3              # Time to stop simulation  [s]
dt = 0.01                 # Time interval            [s]
fps = 30                   # Frames per second        [frame/s]

times = np.arange(0, stop_time, dt)

# Model & Resolution
# ==================

global speed
global b_src

m.add_source(room_width/2, room_width/2, 90)
b_src = m.sources[0]

h  = []
hv = []
ha = []
ht = []

def compute_velocity(dt):
    """Compute velocity of the wave"""
    global speed
    global b_src
    
    watching_particle = m.particles[len(m.particles) // 2][len(m.particles) - 4]
    if 't' not in compute_velocity.__dict__: compute_velocity.t = 0
    if 'max_d' not in compute_velocity.__dict__: compute_velocity.max_d = 0
    compute_velocity.t += dt
    
    d = watching_particle.p.distance(watching_particle.p_0)
    h.append(d)
    hv.append(watching_particle.v.x)
    ha.append(watching_particle.a.x)
    ht.append(compute_velocity.t)

    if d >= compute_velocity.max_d:
        compute_velocity.max_d = d
        speed = b_src.particle.p_0.distance(watching_particle.p_0) / compute_velocity.t
        return False
    elif watching_particle.p.x - watching_particle.p_0.x < 1e-9:
        return False
    else:
        return True
        
global avg_ctime; avg_ctime = 0
global cpt; cpt = 0
particles_history = []
sources_history = []

def compute_all():
    compute_velocity.t = 0
    compute_velocity.max_d = 0
    print_progress_bar(0, times[-1], prefix = 'Compute progress:', suffix = 'Complete', length = 50)
    for t in times:
        if t >= start_time:
            max_reached = compute_step(t)
            if max_reached:
                print_progress_bar(times[-1], times[-1], prefix = 'Compute progress:', suffix = 'Complete', length = 50)
                break
            print_progress_bar(t, times[-1], prefix = 'Compute progress:', suffix = 'Complete', length = 50)

def compute_step(t):
    global avg_ctime
    global cpt
    global b_src

    # Get next source state
    if t < start_time + 0.5/f:
        b_src.particle.p.x = b_src.particle.p_0.x + a/2*(1 + sin((t-start_time) * 2*pi*f))
    else:
        m.clear_sources()

    # Solve step time
    tstart = time()
    m.resolution(dt)
    tstop = time()

    # Add frame to history
    if deferred_mode and display:
        particles_history.append(deepcopy(m.particles))
        sources_history.append(deepcopy(m.sources))
        
    # Update average time
    avg_ctime += tstop - tstart
    cpt += 1
    
    max_reached = compute_velocity(dt)
    return max_reached

def get_speed(k):
    global b_src
    global speed
    global m
    
    #m.k = k
    mass = dx**3 * rho
    m = Map(room_width, dx, k, mass)
    m.add_source(room_width/2, room_width/2, 90)
    b_src = m.sources[0]  
    compute_all()
    
    return speed
    
# Problem optimization
# ====================

f = 1
a = 1

def cost(k):
    return (get_speed(k)-340)**2

t0, t1 = 0.01, 100
learning_rate = lambda n: t0 / (n + t1) #0.000001
initial_k = k
nber_iterations = 10
step_k = 0.005  #0.03

def partial_derivatives(k):
    new_k = k + step_k
    error_k = cost(k)
    new_error_k = cost(new_k)
    derivate_cost_k = (new_error_k - error_k)/step_k

    return derivate_cost_k
 
def gradient_descent():
    global speed
    temporary_k = initial_k
    for i in range(nber_iterations):
        new_k = temporary_k - learning_rate(i) * partial_derivatives(temporary_k)
        temporary_k = new_k
        print("k: {0}, speed: {1}".format(temporary_k, speed))
    return temporary_k       
 
if gradient:
    final_k = gradient_descent() 
    print("After {0} iterations k = {1}".format(nber_iterations, final_k))
    

# Solutions display
# =================

if display:
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
    print("Wave velocity: {}, dx: {}, k: {}".format(speed,dx,k))

    plt.plot(ht, h)
    #plt.plot(ht, hv)
    #plt.plot(ht, ha)
    plt.show()
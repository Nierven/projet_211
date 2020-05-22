import numpy as np
from copy import deepcopy

from utilities import Point, get_threads_count
from threading import Thread

class Task(Thread):
    def __init__(self, batch, start_index, columns):
        Thread.__init__(self)

        self.batch = batch
        self.c = columns
        self.start_index = start_index
        self.func = None

    def run(self):
        for i, particle in enumerate(self.batch):
            if not particle.is_wall and not particle.is_src:
                index = self.start_index + i
                self.func(particle, index//self.c, index%self.c)

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        return result

class Map:
    def __init__(self, width, dx, k, m):
        r = np.arange(0, width + dx, dx)
        self.particles = [[Particle(Point(x, y), m) for x in r] for y in r]
        self.particles = np.array(self.particles)
        self.w = len(self.particles) - 1

        # Handle sources & walls
        self.sources = []
        for (i, j), b in np.ndenumerate(self.particles):
            b.is_wall = i == 0 or j == 0 or i == self.w or j == self.w

        # Constants
        self.k = k
        self.e = dx

        # Creation of all threads
        self.compute_threads = []
        number_of_threads = get_threads_count()
        particles_per_thread = self.particles.size // number_of_threads
        columns = self.particles.shape[1]

        particles_list = self.particles.ravel()

        for i in range(number_of_threads):
            nb_i = i * particles_per_thread
            nb_f = (i + 1) * particles_per_thread - 1
            if self.particles.size - nb_f < particles_per_thread:
                stop = self.particles.size - 1

            batch = particles_list[nb_i:nb_f+1]

            self.compute_threads.append(Task(batch, nb_i, columns))

    def add_source(self, x, y, level):
        """Add source using class Source"""

        # Get local coordinates & source particle
        i = int(y / self.e)
        j = int(x / self.e)
        p = self.particles[i][j]
        p.is_src = True

        # Add corresponding particle as source
        self.sources.append(Source(p, level))

    def clear_sources(self):
        for src in self.sources:
            src.particle.is_src = False
        self.sources.clear()

    def resolution(self, dt):
        """Resolve the current timestep problem"""
        k = self.k
        e = self.e
        w = self.w

        def compute_positions(b, i, j):
            """Compute and assign particle acceleration then update next positions"""
            # Neighbors
            bl = self.particles[i][j - 1]
            bt = self.particles[i + 1][j]
            br = self.particles[i][j + 1]
            bb = self.particles[i - 1][j]

            # Distances with neighbors
            bld = b.p - bl.p
            btd = b.p - bt.p
            brd = b.p - br.p
            bbd = b.p - bb.p
                
            # Forces of springs
            fl_x = -k * (bld.x - e)
            fl_y = -k * (bld.y)
            ft_x = -k * (btd.x)
            ft_y = -k * (btd.y + e)
            fr_x = -k * (brd.x + e)
            fr_y = -k * (brd.y)
            fb_x = -k * (bbd.x)
            fb_y = -k * (bbd.y - e)
                
            # Update of ball acceleration & positions
            b.a.x = (fl_x + ft_x + fr_x + fb_x) / b.m
            b.a.y = (fl_y + ft_y + fr_y + fb_y) / b.m
            b.v.x += b.a.x * dt
            b.v.y += b.a.y * dt
            b.p_1.x = b.p.x + b.v.x * dt
            b.p_1.y = b.p.y + b.v.y * dt

        def update_new_positions(b, i, j):
            """Update position values using p_1 values"""
            b.p.x = b.p_1.x
            b.p.y = b.p_1.y

        def launch_threads(func):
            # Assign function used for tasks then launch them
            for thread in self.compute_threads:
                thread.func = func
                if thread._started._flag: thread.run()
                else: thread.start()
            
            for thread in self.compute_threads:
                thread.join()
            
        launch_threads(compute_positions)
        launch_threads(update_new_positions)

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        return result

    def __str__(self):
        return repr(self.particles)

    def __repr__(self):
        return '<map>: ' + str(self)

class Particle:
    def __init__(self, p_0, m):
        self.p_0 = Point(p_0.x, p_0.y)
        self.p   = Point(p_0.x, p_0.y)
        self.p_1 = Point(p_0.x, p_0.y)

        self.a = Point(0, 0)
        self.v = Point(0, 0)

        self.is_wall = False
        self.is_src = False

        self.m = m

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        return result

    def __str__(self):
        return str(self.p)

    def __repr__(self):
        return '<particle>: ' + str(self)

class Source:
    def __init__(self, particle, level):
        self.particle = particle
        self.level = level

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        return result

    def __str__(self):
        return '{}, {} dB'.format(self.particle.p, self.level)

    def __repr__(self):
        return '<source>: ' + str(self)
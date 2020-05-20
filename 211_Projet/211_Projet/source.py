import numpy as np

class source:
    def __init__(self, p, level):
        self.p = p
        self.level = level

class map:
    def __init__(self, width, dx, k, m):
        r = np.arange(0, width + dx, dx)
        self.boules = [[boule(point(x, y), m, k) for x in r] for y in r]
        self.boules = np.array(self.boules)
        w = len(self.boules) - 1

        # Handle walls
        for (i, j), b in np.ndenumerate(self.boules):
            is_wall = i == 0 or j == 0 or i == w or j == w
            if is_wall:
               b.m = 1e9
               b.k = k * 10

        # Handle sources
        self.sources = []

        # Constants
        self.k = k
        self.e = dx

    def resolution(self, dt):
        """Resolve the current timestep problem"""
        k = self.k
        e = self.e
        w = len(self.boules) - 1

        for (i, j), b in np.ndenumerate(self.boules):

            # Wall handling
            is_wall = i == 0 or j == 0 or i == w or j == w
            is_src  = False
            for src in self.sources:
                if b == src:
                    is_src = True

            if not is_wall and not is_src:
                # Neighbors
                bl = self.boules[i][j - 1]
                bt = self.boules[i + 1][j]
                br = self.boules[i][j + 1]
                bb = self.boules[i - 1][j]
                
                fl_x = -k * ((b.p - bl.p).x - e)
                fl_y = -k * ((b.p - bl.p).y)
                ft_x = -k * ((b.p - bt.p).x)
                ft_y = -k * ((b.p - bt.p).y + e)
                fr_x = -k * ((b.p - br.p).x + e)
                fr_y = -k * ((b.p - br.p).y)
                fb_x = -k * ((b.p - bb.p).x)
                fb_y = -k * ((b.p - bb.p).y - e)
                
                # Update of ball acceleration
                b.a.x = (fl_x + ft_x + fr_x + fb_x) / b.m
                b.a.y = (fl_y + ft_y + fr_y + fb_y) / b.m
            else:
                b.a.x = 0
                b.a.y = 0

        for (i, j), b in np.ndenumerate(self.boules):
            # Update of integred values
            b.v.x += b.a.x * dt
            b.v.y += b.a.y * dt
            b.p.x += b.v.x * dt
            b.p.y += b.v.y * dt

    def __str__(self):
        return repr(self.boules)

    def __repr__(self):
        return str(self)

class boule:
    def __init__(self, p_0, m, k):
        self.p_0 = point(p_0.x, p_0.y)
        self.p = point(p_0.x, p_0.y)
        self.p_1 = point(p_0.x, p_0.y)

        self.a = point(0, 0)
        self.v = point(0, 0)

        self.m = m
        self.k = k

    def __str__(self):
        return repr(self.p)

    def __repr__(self):
        return str(self)

class point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __sub__(self, p):
        return point(self.x - p.x, self.y - p.y)

    def distance(self, p):
        dist_x = p.x - self.x
        dist_y = p.y - self.y
        return np.sqrt(dist_x**2 + dist_y**2)

    def angle(self, p):
        dist_x = p.x - self.x
        dist_y = p.y - self.y
        return np.arctan2(dist_y, dist_x)

    def __str__(self):
        return '({};{})'.format(self.x, self.y)

    def __repr__(self):
        return str(self)
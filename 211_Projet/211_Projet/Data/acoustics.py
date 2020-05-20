filename_mesh = 'room.mesh'

v_n = 1.0 # m/s
w = 1000.0
c = 343.0 # m/s
rho = 1.55 # kg/m^3

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'ts' : 'ts',
    
    'save_times' : 'all',
    'active_only' : False
}

materials = {
    'one' : ({'one' : 1.0},),
}

regions = {
    'Omega' : 'all',
    'Gamma_in' : ('vertices in (x < 0.01)', 'facet'),
    'Gamma_out' : ('vertices in (x > 0.99)', 'facet'),
}

fields = {
    'accoustic_pressure' : ('complex', 1, 'Omega', 1),
}

variables = {
    'p'   : ('unknown field', 'accoustic_pressure', 0),
    'q'   : ('test field',    'accoustic_pressure', 'p'),
}

# Essential boundary conditions
ebcs = {
    'p1': ('Gamma_in', [(1.0, 2.0)], {'p.0' : 2})
}

# Integrals specify which numerical scheme to use.
# Here we are using a 2nd order quadrature over a 2 dimensional space.
integrals = {
    'i' : 2,
}

equations = {
    'Acoustic pressure' :
    """dw_laplace.i.Omega( one.one, q, p )
    = 0"""
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max'      : 1,
        'eps_a'      : 1e-1,
        'eps_r'      : 1.0,
        'macheps'   : 1e-16,
        'lin_red'    : 1e-1, # Linear system error < (eps_a * lin_red).
        'ls_red'     : 0.1,
        'ls_red_warp' : 0.001,
        'ls_on'      : 1.1,
        'ls_min'     : 1e-5,
        'check'     : 0,
        'delta'     : 1e-6,
    }),
    'ts' : ('ts.simple', {
        't0'     : 0.0,
        't1'     : 4.0,
        'dt'     : None,
        'n_step' : 5, # has precedence over dt!

        'quasistatic' : True,
        'verbose' : 1,
    })
}
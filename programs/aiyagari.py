
import numpy as np
import math
from scipy.optimize import fminbound, brentq
from lininterp import LinInterp

# Parameters, a la Aiyagari 
w = 1.3712
r =  0.0129
beta = 0.96
mu = 2
sigma = 0.4
m = 1 - mu

m = 1 - mu
def U(c): 
    return ((c + 1e-10)**m - 1) / m

def Up(c): # derivative
    return (c + 1e-10)**(-mu)

# The L shock has 3 states: {smin, 1, smax}, and is uniform.  We
# assume that |smin -1| = |smax - 1| = d, so that E L = 1, and 
# Var L = (2/3) * d**2.  We choose d so that Var L = sigma**2
d = math.sqrt(3.0 / 2.0) * sigma
smin, smax = 1 - d, 1 + d

gridmin, gridmax, gridsize = w * smin, 14, 150
grid = np.linspace(gridmin, gridmax, gridsize)

def maximizer(h, a, b):
    return fminbound(lambda x: -h(x), a, b)

# The objective function for the v-greedy problem
def obj_fun(z, a, v):
    t = v(w * smin + (1 + r) * a) \
                + v(w + (1 + r) * a) + v(w * smax + (1 + r) * a)
    return U(z - a) + (beta / 3) * t

def bellman(v):
    """
    The approximate Bellman operator.
    Parameters: v is a vectorized function (i.e., a callable object which 
        acts pointwise on arrays).
    Returns: An instance of LinInterp.
    """
    vals = []
    for z in grid:
        h = lambda a: obj_fun(z, a, v)
        vals.append(h(maximizer(h, 0, z)))
    return LinInterp(grid, vals)


def find_zb(policy_function):
    h = lambda z: - Up(z) + (beta * (1 + r) / 3) * (\
              Up(w * smin - policy_function(w * smin))  \
            + Up(w - policy_function(w))  \
            + Up(w * smax - policy_function(w * smax)) )
    return brentq(h, gridmin, gridmax)


def compute_optimal_policy(read_from_file=True):

    if read_from_file:
        policy = LinInterp(np.loadtxt('grid_dat.txt'), np.loadtxt('pol_dat.txt'))

    else:
        v = U          # Initial condition
        tol = 0.1      # Error tolerance 

        while 1:
            new_v = bellman(v)
            err = np.max(np.abs(new_v(grid) - v(grid)))
            print err
            if err < tol:
                break
            v = new_v

        policy_on_grid = []
        for z in grid:
            h = lambda a: obj_fun(z, a, v)
            policy_on_grid.append(maximizer(h, 0, z))

        policy = LinInterp(grid, policy_on_grid)

        # for z < zb, set policy(z) = 0.0
        zb = find_zb(policy)
        print "zb = ", zb
        for i in range(gridsize):
            if grid[i] < zb:
                policy_on_grid[i] = 0.0

        np.savetxt('grid_dat.txt', grid, fmt='%.10f')
        np.savetxt('pol_dat.txt', policy_on_grid, fmt='%.10f')

        policy.Y = policy_on_grid

    return policy



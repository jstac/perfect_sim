
from aiyagari import *
import random

def draw_shock():
    return random.choice((smin, 1, smax))

increment = 50  # Test at multiples of this value

policy = compute_optimal_policy(read_from_file=False)
zb = find_zb(policy)

def update(z, L):
    return w * L + (1 + r) * policy(z)

# Compute upper bound of state space
#h = lambda z: update(z, smax) - z
#zmax = brentq(h, gridmax * 0.5, gridmax * 1.5)
zmax = gridmax

def compute_yt(shock_path):
    z = zmax
    regenerated = False
    indexes = range(len(shock_path))
    indexes.reverse()
    for t in indexes:
        z = update(z, shock_path[t])
        if z < zb and t > 0:
            regenerated = True
    return z, regenerated

def draw_perfect():
    shock_path = [draw_shock() for i in range(increment)]
    while 1:
        z, regenerated = compute_yt(shock_path)
        if regenerated:
            break
        else:
            shock_path.extend([draw_shock() for i in range(increment)])
    return z, len(shock_path)


# psi_solve/potentials.property
import numpy as np

# --- Potential Play Grounds ---
Last_k_value = None

def constant(x,c):
    """Adds a flat baseline potential."""
    return np.ones_like(x)*c

def harmonic(x,k,center=L/2):
    """A Parabola"""
    global Last_k_value
    Last_k_value = k
    
    constant = 1 #m*L**2
    potential = 0.5*k*(x - center)**2 # Parabola
    return constant*potential

def gaussian_well(x, center=L/2, width=1.0, depth=50): # A new Potnetial
    """A dip in the potential (finite well)."""
    return -depth * np.exp(-(x - center)**2 / (2 * width**2))

def inf_sqaure_well(x,lower_bound,upper_bound):
    "Gives you an Infinite square well of Length L"
    "Takes the Grind Points array x and creates an array of Zeros in the same shape"
    ""
    #walls_array.append(lower_bound)
    #walls_array.append(upper_bound)

    V = np.zeros_like(x) 
    V[x<lower_bound] = np.inf # Potnetial of infinite below Lower Bound
    V[x>upper_bound] = np.inf # Potnetial of infinte above upper Bound
    # Effetively, all untouced grind points have potential Zero

    return V

def inf_wall(x,side,bound):
    """
    Places an 'Infinite' wall.
    Instead of np.inf, we use 1e10 (10 billion). This is called Penalty Method
    This forces psi -> 0 without breaking the solver.
    """
    V = np.zeros_like(x)
    HUGE_NUMBER = 9e10 # Replacment for np.inf

    #walls_array.append(bound)
    side = side.strip(', . ').lower() # Becase Strings are immutable and .lower() is a method

    if side =='left':
        # Wall exists everywhere to the left of 'bound'
        V[x<bound] = HUGE_NUMBER

    elif side =='right':
        # Wall exists everywhere to the right of 'bound'
        V[x>bound] = HUGE_NUMBER

    return V


def finite_barrier(x, center, width, height):
    """A square block in the middle.""" # Maybe we can see Tunneling??
    V = np.zeros_like(x)
    mask = (x > (center - width/2)) & (x < (center + width/2))
    V[mask] = height
    return V


def V_double_well(x, depth=20, separation=1):
    V = depth * ( (x-L/2)**2 - separation )**2
    return V


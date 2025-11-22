# psi_solve/constants.py

import numpy as np

# Grid Parameters

"Limitation: only 0 to L now, not -a to a"
L = 50

N = 2000 # Number of Internal Grid Points
#x = np.linspace(0,L,N+1) # N+1, because if you need N intervals we need N+1 end points |..|..|
x = np.linspace(-L/2,L/2,N+2) # For also Wall # N+1, because if you need N intervals we need N+1 end points |..|..|

# x0 , x1 ... xN (Thats N+1 end Points)
dx = x[1] - x[0] # a Constant Grid Space


# Constants in AMU
hbar = 1.0 # Reduced Planks Constant
m = 1.0


import numpy as np
from .constants import hbar, m, dx, N


def Hamiltonian(V):

    N
    # The Differntial operator
    main_diagonal = (1/dx**2)*np.diag(-2*np.ones(N))
    off_diagonal1 = (1/dx**2)*np.diag(np.ones(N-1),-1)
    off_diagonal2 = (1/dx**2)*np.diag(np.ones(N-1),1)
    D2 =  (main_diagonal + off_diagonal1 + off_diagonal2)
    
    # Hamiltonian
    #T = (-(hbar**2/2*m ) * D2 ) # Recall D is the Differntial Operator 

    # Correcting for Python Left to Right Product order
    T = (-(hbar**2/(2*m) ) * D2 ) # Recall D is the Differntial Operator 

    H = T + np.diag(V[1:-1])

    return H


def normalize(psi):
    """Normalizes the eigenstates using the Riemann sum approximation."""
    for i in range(psi.shape[1]):
        # 1. Find the integral of probability density
        # sum( |psi|^2 * dx )
        norm_factor = np.sum(psi[:, i]**2) * dx
        
        # 2. Divide the vector by the square root of that integral
        psi[:, i] = psi[:, i] / np.sqrt(norm_factor)


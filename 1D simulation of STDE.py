# 1 Demension Time dependent Simulation of Schrödinger Equation

# Background Infromation:
# The Time dependent Schrödinger Deals with the Matter Waves ( Also known 
# as de-Broglie Waves) ( λ = h/p)

# At all scales where measurements have been practical, matter exhibits wave-like behavior. 
# For example, a beam of electrons can be diffracted just like a beam of light 
# or a water wave.The concept that matter behaves like a wave was proposed by 
# French physicist Louis de Broglie. Hence the name

# In λ = h/p, The de-Broglie Wavelength λ is associated with a particle of Momentum p through 
# the plan constant h

# The Wave Equation ( The one below?) is important becase it explains why an electron does not 
# simply orbit into the nucleus of the atom. This is the foundation of QM.The answer lies in thinking
# of the electron not as a particle but also as a wave, hence we need a wave equation that
# represents this. 

# Ψ is the Wave Function that represents the de-Broglie / matter waves.
# 
# For a freely moving particle moving in +x direction, the equation that represents it is:
# Ψ = A e^-i(wt-kx). As Time is postive, x is postive hence particel movies right. 


# Side Note: Remember that Ψ = e^(iαx) = cos(αx) + i sin(αx)



# Remember the Formula for TOTAL Energy of a Neutral Ball/ or any Human scale Obejct?
# It has both Kinetic and Potnetial energy in its simplest form 
# Total Energy = 1/2 mv^2 - mgh (So that the total energy at equilibrium is Zero)
# Total Energy = p^2/2m - mgh = -p^2/2m + mgh  

# || Similarly, the Total Energy For a single non-relativistic spinless particle of mass, 
# in a time-invariant potential is: 

# = 1/2m ∂^2/∂x^2 + V(x)
# (ASK FEW QUESTION, Velocity = dx/dt, this is differnet. 
# why is m in denominator? Because this is in the Momentum Form where p = mv 
# Why KE part looks like it needs to be acted on soemthing while V(x) looks complete ?
# It looks like that but, both need to act on something. BOTH are operators, that when acting
# on an operand like the wave equation give the total energy. 
# This new operator (sum of other operators) is called the Hamiltonian

# H = 1/2m ∂^2/∂x^2 + V(x) or
# H = 1/2m ∇^2 + V(x). Where ∇^2 = ∂^2/∂x^2, The Laplace Operator. 

# Ok now By Schrodinger we have the Wave Equation: 
# ∂Ψ(x,t)/∂t = -iH Ψ(x,t) , Using Time Evolution and hence
# by using a U Time Evalution Operator and knowing that Ψ(0) initial State is Arbitrary, we
# can find out what the Wave Function is:
# Ψ(x,t) = e^(-iH) Ψ(x,0)


# The Computational Challenge, Becase Computers work For discreate operations. We need to convert
# THis COntinous Function into a discrete one. 
# We do this by splitting the entier function into N parts with equally spaced ∂x interval. This
# Process is called Discretization. The start/0th point is x0 and the last point is xN. 
# And there are N spaces of ∂x

# Hence we can COnvert this contnous line into a vector having N slots/ elements, with each element 
# having the corresponding value of the wave function at that point on the grid. 


# When we convert this into a the Grid, we get: We find the Hamiltonian to be 
# H = -1/2Δx^2 * [a Tri Diagonal Teoplitz Sparse Matrix on NxN (Discrete Laplace Operator)] + 
# [a Diagonal Matrix with Values V0 to VN], 

# Here the Discrete Laplace Operator, DLO is a Tri Diagonal Teoplitz Sparse Matrix on NxN with 
# 1 on the Top , 2 on the Center Diagonal and 1 on the Below Center Diagonal. Again this is 
# Teoplits
# Refer to Notes for more Information

# Let's Construct this Discrete Hamiltonian Function:

import scipy.sparse
import numpy as np

def hamiltonian(N,dx,V=None): 
    """
    

    Args:
        N: Number of Slots you want to divide into
        dx: Spacing of Each Slot
        V: An array containing the Potnetioal at each point N, 
            Must be (N,) format
            Keyword argument V is set to None to represent the Case of No Potentiol Energy
    
    Return:
        Discrete Hamiltonian Matrix of Dimention NxN 

    """

    DLO = scipy.sparse.diags([1,-2,1], offsets=[1,0,-1],shape=(N,N))
    KE = (-1/(2*(dx**2)))* DLO

    H = KE

    # The truth Value of an Array with more than 1 element is apparenlty Ambigious
    if V is not None: # Recall bool(None) = False

        PE = scipy.sparse.spdiags(V,0,N,N)

        H += PE
    
    
    return H.tocsc() # COnverts the Sparse Matrix to Compressed Sparse Matrix (csc)
    # Effecient for efficient for column slicing and for solving linear systems


# Create a Time Evolution Operator

def time_evolution(H, dt):
    H_dense = H.toarray()  # convert sparse CSC to a dense NumPy array
    return scipy.linalg.expm(-1j * H_dense * dt)



# Generate a Simulator Function that gives the Wave FUcntion 



def simulator(psi,H,dt): # Once
    U = time_evolution(H,dt)

    while True:
        yield psi # Get the Last Psi Variable Not Input one
        psi = U @ psi # Matrix Multiplication
        # psi = expm_multiply(U, psi) # More Effecient Multiplication



# Now the Wave Function Ψ(x,t) has both complex and Real part. This can not be used to find 
# the probability of a particle at a particlar point and at a particular time
# Instead we need the Probability Density Function. WHich gives us the probability
# Of finding a particlar Particle in that point in space. This probablity Density FUnction 
# Does a Modulas Sqare |Ψ(x,t)|**2 of the FUnction, which is essentially the sum of the sqared real
# and Imaginary Parts

# Lets make a probability density function for a given Psi


def pdf(psi):
    return (psi.real)**2 + (psi.imag)**2

# Now Lets Finnaly get our Inital Psi, we want a particle that mostly starts and is concenrtaed 
# at a point. We hence use the well known, localized Gaussian wave packet centered 
# at x0 with average initial momentum, p0, and Gaussian width, σ0. 

# Here we can find the distance the particle traveles in time t, by v0*t. 
# Since we set mass to 1, p0 = 1*v0 ->, v0=p0 
# and hence the distance can be found by p0*t. 

# The uncertinity in the particles position is given by the Gaussian Width, Δx = σ0

# Since a the Gaussian we have choosen Saturates/ takes minimum value in the
# position-momentum uncertinity, Δp = 1/(2*σ0)
# As: ΔxΔp = 1/2, σ0*Δp = 1/2 hence. Δp = 1/(2*σ0)
# Its a Minimum uncertinity Wave Packet
# ENergy of this particle is given by E_kinetic = p**2/2m + V
# Since, m is 1 and V for a free particle is 0, here E = p0**2/2 

# We know from experimentation and Schrodinger equation derivation that This 
# gaussian's formula is :

def gaussian(x,x0,sigma0,p0):
    return (1/(2*np.pi*sigma0**2))**(1/4) * np.exp( -(x-x0)**2/(4*sigma0**2) + 1j * p0 *x )



# Now lets Simulate this Free Particle (V=0) 
# A Gaussian wave packet centered at x0=0 with initial momentum 
# of p0=1 and Gaussian width of σ0=5. In the absence of a potential, 
# this will behave like a free particle.

# Lets say we want 2000 Steps and want the simulation to start from positions 0 and go till
# 200. 
# First we need to create an array that goes from 0 to 200 and has 2000 steps inbetween so,
'''
N = 2000

start_x = -10
stop_x = 200
x = np.linspace(start_x,stop_x,N)
'''

# But it cant stop here, You have told it to be evenly spaced and have N points. But you dont 
# know what that spacing dx is going to be between each point. To get this we must say, 
# retstep=True. Hence


N = 2000

#start_x = 0 # Very Intresting, if we start from a wall. x = 0. it immedielty interfeares with itslef. How does it know x = 0 is a wall? idk???
start_x = -50
stop_x = 100
x,dx = np.linspace(start_x,stop_x,N,retstep=True)

# Great now you have an array of all points you want the simulation to run through next.

# Generat an initial Wave Function at x0=0 with initial momentum 
# of p0=1 and Gaussian width of σ0=5.

# x the array will go into this and the output would be for each point in space. 
# Hence this is the Discrete Wave function that represnets the inital condition at each point in 
# space (As much as a discrete system can)
psi0 = gaussian(x,x0=0,sigma0=5,p0=1)

# COntruct our Hamiltonian from the fucntion

H = hamiltonian(N,dx)



###################################################################################


import matplotlib.pyplot as plt
import matplotlib.animation as animation

# --- helper to (softly) renormalize the wavefunction to avoid drift ---
def renormalize(psi, dx):
    norm = np.sqrt(np.sum(np.abs(psi)**2) * dx)
    if norm == 0:
        return psi
    return psi / norm


# Animation Section


# Set up the figure and a single line for probability density
fig, ax = plt.subplots(figsize=(8,4.5))
pen, = ax.plot([], [], lw=2) # draws an empty line on the axes, Initilize the pen????
# create  horizontal line once to show Max peak
peak_line = ax.axhline(y=0.0, color="red", linestyle="dashed", linewidth = 2)  # Line2D
width_line = ax.axhline(y=0.0, color="b", linewidth = 2 )
## ax.hlines does not have Set ydata. 


# Set Canvas Axis 
ax.set_xlim(start_x, stop_x) 
ax.set_ylim(0, 1.1*np.max(pdf(psi0)))  # initial guess; we’ll update on first frame
#x.set_ylim(0, 0.25)  # initial guess; we’ll update on first frame

ax.set_xlabel("x")
ax.set_ylabel(r"$|\Psi(x,t)|^2$")
title = ax.set_title("t = 0.0")


# We will update this psi from the generator
current_psi = None
t = 0.0

# Animation parameters
n_frames = 400          # how many time steps to show
dt = 0.1                # must match the dt used in simulator(...)

#simulation1 = simulator(psi0, H, dt=0.1)

sim = simulator(psi0, H, dt=dt)  # use_denseU=True builds a dense U (heavy)

def init():
    global current_psi, t # Python forces you to declare if you want to edit it, otehrwise it thinks you are making a local variable
    # start from t=0 state from the generator
    current_psi = next(sim)
    current_psi = renormalize(current_psi, dx)
    dens = pdf(current_psi)
    
    pen.set_data(x, dens)
    peak_line.set_ydata([dens.max(), dens.max()])  # set horizontal line at peak
    width_line.set_ydata([0.5*dens.max(),0.5*dens.max()])
    width_line.set_xdata([])
    # How does this work??


    t = 0.0
    title.set_text(f"t = {t:.2f}")
    
    return pen, title , peak_line

def update(n_frames):
    global current_psi, t
    # advance one step in time
    current_psi = next(sim) # This handles the next Data Simulation. No need of frames into update
    current_psi = renormalize(current_psi, dx)
    dens = pdf(current_psi)
    
    pen.set_data(x, dens)
    peak_line.set_ydata([dens.max(), dens.max()])  # update peak line

    t += dt
    title.set_text(f"t = {t:.2f}")
    
    return pen, title,peak_line # With blit= True, Matplotlib expects you to tell what artist to redraw each frame. Here we want only the 

ani = animation.FuncAnimation(
    fig,
    update,               # called once per frame after init()
    frames=n_frames,      # how many frames to render # Un used, kept jsut to supress errors
    init_func=init,       # one-time setup function
    blit=True,            # only re-draw changed artists (faster)
    interval=5,          # ~5 ms between frames (display rate)
    repeat=False
)


#plt.tight_layout()
#plt.show()

''' For Jupyter
from IPython.display import HTML
plt.rcParams["animation.html"] = "jshtml"
HTML(ani.to_jshtml())
'''

from matplotlib.animation import HTMLWriter
import webbrowser
import os

writer = HTMLWriter(embed_frames=True, fps=50)
ani.save("Free_particle.html", writer=writer)
# Now open schrodinger.html in your browser → play/pause slider appears

print("Saving in:", os.getcwd())
webbrowser.open("Free_particle.html")

# New Version
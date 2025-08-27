import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ---------- Parameters ----------
Nx, Ny = 220, 220              # grid
Lx, Ly = 1.0, 1.0
dx = Lx / Nx
dy = Ly / Ny

DA = 1.0e-3                    # diffusion of acid
DB = 1.0e-3                    # diffusion of base
D_dye = 7.5e-4                 # diffusion of dye
dt = 0.20 * dx*dx / (4*max(DA, DB))  # stable explicit step

k_diss = 2.5e-2                # dissolution rate of solids
a_sat = 1.0; b_sat = 1.0       # "solubility" scale (nondimensional)
k_rxn = 4.0                    # acid-base reaction rate
co2_yield = 1.0                # CO2 per unit reaction
dye_yield = 0.10               # dye released per reaction

# Bubble parameters (toy physics)
CO2_spawn_thresh = 0.05
bubble_spawn_cooldown = 10     # frames between spawns per cell
bubble_max = 800
bubble_growth = 0.015
bubble_base_speed = 0.12       # upward speed ~ bubble_base_speed * sqrt(R)

# ---------- Grids ----------
y, x = np.mgrid[0:Ny, 0:Nx]
cx, cy = Nx//2, Ny//2
R0 = int(0.18*Nx)
mask_bomb = (x-cx)**2 + (y-cy)**2 <= R0**2

# Solid reservoirs (dimensionless mass)
sA = np.zeros((Ny, Nx)); sA[mask_bomb] = 1.0
sB = np.zeros((Ny, Nx)); sB[mask_bomb] = 1.0

# Aqueous concentrations
A = np.zeros((Ny, Nx), dtype=np.float32)
B = np.zeros((Ny, Nx), dtype=np.float32)
DYE = np.zeros((Ny, Nx), dtype=np.float32)
CO2 = np.zeros((Ny, Nx), dtype=np.float32)

# Spawn cooldown map
cooldown = np.zeros((Ny, Nx), dtype=np.int16)

# ---------- Helpers ----------
def laplacian(U):
    return (np.roll(U,1,0) + np.roll(U,-1,0) + np.roll(U,1,1) + np.roll(U,-1,1) - 4*U) / (dx*dx)

# Bubble store: positions (float), radii
b_pos = []  # list of [y, x]
b_R = []    # list of radii

rng = np.random.default_rng(0)

# ---------- Visualization setup ----------
fig, ax = plt.subplots(figsize=(6,6))
im = ax.imshow(DYE, origin='lower', vmin=0, vmax=0.8, interpolation='bilinear')
ax.set_xticks([]); ax.set_yticks([])
ax.set_title("Bath bomb reaction–diffusion + bubbles (toy)")

bubble_scat = ax.scatter([], [], s=[], edgecolors='white', facecolors='none', linewidths=1.0)

# ---------- Time stepping ----------
def step():
    global A, B, DYE, CO2, sA, sB, cooldown, b_pos, b_R

    # 1) Dissolution (only where there is solid)
    #    Driving force: undersaturation (1 - A/a_sat), clipped to [0,1]
    undersatA = np.clip(1.0 - A/a_sat, 0.0, 1.0)
    undersatB = np.clip(1.0 - B/b_sat, 0.0, 1.0)
    JA = k_diss * sA * undersatA
    JB = k_diss * sB * undersatB

    # solid depletion
    sA = np.maximum(sA - JA*dt, 0.0)
    sB = np.maximum(sB - JB*dt, 0.0)

    # 2) Reaction (A + B -> CO2 + salt), rate = k_rxn * A * B
    R = k_rxn * A * B
    # Do not consume more than available in a single step
    R = np.minimum(R, np.minimum(A/dt + JA, B/dt + JB))  # conservative cap

    # 3) Diffusion + source/sink updates
    A = A + dt*(DA*laplacian(A)) + dt*JA - dt*R
    B = B + dt*(DB*laplacian(B)) + dt*JB - dt*R
    DYE = DYE + dt*(D_dye*laplacian(DYE)) + dye_yield*dt*R
    CO2 = CO2 + co2_yield*dt*R

    # Non-negative clamp
    A = np.clip(A, 0, None)
    B = np.clip(B, 0, None)
    DYE[:] = np.clip(DYE, 0, 1.0)

    # 4) Bubble spawning (throttled)
    cooldown = np.maximum(cooldown-1, 0)
    spawn_mask = (CO2 > CO2_spawn_thresh) & (cooldown == 0)
    # Sparse random thinning so we don't spawn too many
    if spawn_mask.any() and len(b_pos) < bubble_max:
        ys, xs = np.where(spawn_mask)
        take = rng.choice(len(ys), size=min(50, len(ys)), replace=False)
        for yi, xi in zip(ys[take], xs[take]):
            b_pos.append([float(yi), float(xi)])
            b_R.append(0.8 + 0.4*rng.random())  # small initial radius
            cooldown[yi, xi] = bubble_spawn_cooldown
            # Local CO2 consumption into bubbles (visual bookkeeping)
            CO2[yi, xi] *= 0.5

    # 5) Bubble update: rise, grow a bit, cull at top
    if b_pos:
        for i in range(len(b_pos)):
            Rb = b_R[i]
            speed = bubble_base_speed * np.sqrt(max(Rb, 1e-3))
            b_pos[i][0] += speed    # y increases upward in our plot? (origin='lower') -> increasing y = up
            b_R[i] = Rb + bubble_growth*(0.5 + 0.5*np.tanh(CO2[int(b_pos[i][0])%Ny, int(b_pos[i][1])%Nx]))
        # Remove bubbles that left domain or got too big
        keep = []
        for i,(yy,xx) in enumerate(b_pos):
            if yy < Ny-1 and b_R[i] < 12.0:
                keep.append(i)
        b_pos = [b_pos[i] for i in keep]
        b_R  = [b_R[i]  for i in keep]

def animate(_):
    step()
    im.set_data(DYE)
    if b_pos:
        Y = np.array([p[0] for p in b_pos])
        X = np.array([p[1] for p in b_pos])
        S = (np.array(b_R)*2.0)**2  # marker size ~ radius^2
        bubble_scat.set_offsets(np.c_[X, Y])
        bubble_scat.set_sizes(S)
    else:
        bubble_scat.set_offsets(np.empty((0,2)))
        bubble_scat.set_sizes([])
    return (im, bubble_scat)

ani = FuncAnimation(fig, animate, interval=30, blit=False)
plt.show()

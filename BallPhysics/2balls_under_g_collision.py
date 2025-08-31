import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, HTMLWriter
from matplotlib.patches import Circle
import webbrowser, os


wall_below = 0 #m 
wall_right = 25 #m
wall_left = -5

T = 40 # s, total time
dt = 0.01 # s, Sampling time

N = int(T/dt) # How many times the Function should run

# Modify drop_ball to include two ball collisions
# Lets start with 1D Collison in y

# Lets Add Horizontal Velocty and Make this Guys 2D and add walls on side

A_y0 = 5 # Start at height above origin
A_ys = [A_y0]
A_x0 = 0 # Start x at origin
A_xs = [A_x0]
A_radius = 0.5 #m , Radius of Ball


B_y0 = 10 # Start at height above origin
B_ys = [B_y0]
B_x0 = 0 # Start x at origin
B_xs = [B_x0]
B_radius = 0.5 #m , Radius of Ball

A_vy0 = 0 # m/s down
A_a = -9.81 #m/s^2, Gravity on Earth
A_vx0 = -2 # m/s

A_r = 0.9 # coefficient of restitution, 
# r=1: perfectly elastic → no energy loss → bounces forever.
# r<1: inelastic → each bounce smaller.
# r=0: perfectly inelastic → the ball sticks to the ground
B_vy0 = 0 # m/s down
B_a = -9.81 #m/s^2, Gravity on Earth

B_vx0 = 20 #m/s right
B_r = 0.9 # coefficient of restitution, 
# r=1: perfectly elastic → no energy loss → bounces forever.
# r<1: inelastic → each bounce smaller.
# r=0: perfectly inelastic → the ball sticks to the ground

collision = False


for i in range(N):
    A_ds = A_vy0*dt + (1/2)*A_a*(dt)**2
    A_y0 = A_y0 + A_ds

    A_x0 = A_x0 + A_vx0*dt

# Seems like i have to update the velocities as well because otherwise, 
# there just seems to a constant Velocity kinda of effect from acceleration
    A_vf = A_vy0 + A_a*dt
    A_vy0 = A_vf

    A_ys.append(A_y0)
    A_xs.append(A_x0)


    B_ds = B_vy0*dt + (1/2)*B_a*(dt)**2
    B_y0 = B_y0 + B_ds

    B_x0 = B_x0 + B_vx0*dt

    B_vf = B_vy0 + B_a*dt
    B_vy0 = B_vf

    B_ys.append(B_y0)
    B_xs.append(B_x0)

    # if I ask the radius positins to be eqaul, that never happenes becase in the
    # Time step, when the quantaties are a float, they will never exactly equal
    # Each other

    # Lets Shorted this logic using net distance < Combined Radius
    # Then Lets Make this a 2D collision

    net_distance = np.sqrt((A_y0 - B_y0)**2 + (A_x0 - B_x0)**2)
    if net_distance < A_radius+B_radius:
        A_vy0 = -A_r * A_vy0
        B_vy0 = -B_r * B_vy0

        A_vx0 = -A_r * A_vx0
        B_vx0 = -B_r * B_vx0



    if (A_y0-A_radius) <= wall_below:
        collision = True
        A_vy0 = -A_r * A_vy0

    if (B_y0-B_radius) <= wall_below:
        collision = True
        B_vy0 = -B_r * B_vy0

    if (A_x0 + A_radius) >= wall_right:
        A_vx0 = -A_r * A_vx0

    if (B_x0 + B_radius) >= wall_right:
        B_vx0 = -B_r * B_vx0

    if (A_x0 - A_radius) <= wall_left:
        A_vx0 = -A_r * A_vx0

    if (B_x0 - B_radius) <= wall_left:
        B_vx0 = -B_r * B_vx0

    
"""
A_xs,A_ys = drop_ball(x0=A_x0,y0=A_y0,v0=A_vy0,radius=A_radius,a=A_a,r=A_r)

B_xs,B_ys = drop_ball(x0=B_x0,y0=B_y0,v0=B_vy0,radius=B_radius,a=B_a,r=B_r)

"""


# ---------------- Figure: draw a real ball ----------------
fig, ax = plt.subplots(figsize=(8, 7))
#ax.set_xlim(min(A_xs[0],B_xs[0])-5 ,max(A_xs[0],B_xs[0])+5)
ax.set_xlim(-10 ,30)

ax.set_ylim(0,20)
ax.set_aspect('equal', adjustable='box') #  if not equal, the Output is sqiushed


# ALl walls
ax.axvline(x=wall_left,color="Black",lw=3)
ax.axvline(x=wall_right,color="Black",lw=3)

# Labels 
ax.set_xlabel("x (m)")
ax.set_yticks([])
ax.set_title("2 Balls bouncing Under Gravity")



A_ball = Circle(xy=[A_xs[0],A_ys[0]],radius=A_radius,color="Blue")
B_ball = Circle(xy=[B_xs[0],B_ys[0]],radius=B_radius,color="Black")
ax.add_patch(A_ball)
ax.add_patch(B_ball)


def init():
    A_ball.center = (A_xs[0],A_ys[0])
    B_ball.center = (B_xs[0],B_ys[0])
    return (A_ball,B_ball)
#Your animation callbacks must return an 
# iterable of artists when blit=True. So not (ball,), not ball

def update(i):
    A_ball.center = (A_xs[i],A_ys[i])
    B_ball.center = (B_xs[i],B_ys[i])
    return (A_ball,B_ball)

ani = FuncAnimation(fig, update, init_func=init, frames=N, interval=1000*dt, blit=True)

# --- Save to standalone HTML and show ---
outfile = "2balls_fall_gravity.html"
ani.save(outfile, writer=HTMLWriter(fps=int(round(1/dt)), embed_frames=True))
plt.close(fig)
print("Saved:", os.path.abspath(outfile))
webbrowser.open("file://" + os.path.abspath(outfile))

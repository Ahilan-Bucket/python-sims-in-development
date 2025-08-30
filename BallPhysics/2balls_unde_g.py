import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, HTMLWriter
from matplotlib.patches import Circle
import webbrowser, os


# Define Parameters



wall_below = 0 #m 

T = 60 # s, total time
dt = 0.02 # s, Sampling time

N = int(T/dt) # How many times the Function should run

def drop_ball(y0,v0,radius,a,r):
    collision = False
    ys = []
    for _ in range(N):
        if (y0-radius) <= wall_below:
            collision = True
            v0 = -r * v0

    # Seems like i have to update the velocities as well because otherwise, 
    # there just seems to a constant Velocity kinda of effect from acceleration

        ds = v0*dt + (1/2)*a*(dt)**2
        y0 = y0 + ds

        vf = v0 + a*dt
        v0 = vf
        ys.append(y0)
    return ys


A_y0 = 10 # Start at height above origin
A_x0 = 0 # Start x at origin
A_v0 = 0 # m/s down
A_radius = 0.5 #m , Radius of Ball
A_a = -9.81 #m/s^2, Gravity on Earth
A_r = 0.9 # coefficient of restitution, 
# r=1: perfectly elastic → no energy loss → bounces forever.
# r<1: inelastic → each bounce smaller.
# r=0: perfectly inelastic → the ball sticks to the ground



A_ys = drop_ball(y0=A_y0,v0=A_v0,radius=A_radius,a=A_a,r=A_r)


B_y0 = 5 # Start at height above origin
B_x0 = 5 # Start x at origin
B_v0 = 0 # m/s down
B_radius = 0.5 #m , Radius of Ball
B_a = -9.81 #m/s^2, Gravity on Earth
B_r = 0.9 # coefficient of restitution, 
# r=1: perfectly elastic → no energy loss → bounces forever.
# r<1: inelastic → each bounce smaller.
# r=0: perfectly inelastic → the ball sticks to the ground

B_ys = drop_ball(y0=B_y0,v0=B_v0,radius=B_radius,a=B_a,r=B_r)


# ---------------- Figure: draw a real ball ----------------
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xlim(min(A_x0,B_x0)-5 ,max(A_x0,B_x0)+5)
ax.set_ylim(0,20)
ax.set_aspect('equal', adjustable='box') #  if not equal, the Output is sqiushed

# Labels 
ax.set_xlabel("x (m)")
ax.set_yticks([])
ax.set_title("2 Balls bouncing Under Gravity")




A_ball = Circle(xy=[A_x0,A_ys[0]],radius=A_radius)
B_ball = Circle(xy=[B_x0,B_ys[0]],radius=B_radius)
ax.add_patch(A_ball)
ax.add_patch(B_ball)


def init():
    A_ball.center = (A_x0,A_ys[0])
    A_ball.center = (B_x0,B_ys[0])
    return (A_ball,A_ball)
#Your animation callbacks must return an 
# iterable of artists when blit=True. So not (ball,), not ball

def update(i):
    A_ball.center = (A_x0,A_ys[i])
    B_ball.center = (B_x0,B_ys[i])
    return (A_ball,B_ball)

ani = FuncAnimation(fig, update, init_func=init, frames=N, interval=1000*dt, blit=True)

# --- Save to standalone HTML and show ---
outfile = "2balls_fall_gravity.html"
ani.save(outfile, writer=HTMLWriter(fps=int(round(1/dt)), embed_frames=True))
plt.close(fig)
print("Saved:", os.path.abspath(outfile))
webbrowser.open("file://" + os.path.abspath(outfile))

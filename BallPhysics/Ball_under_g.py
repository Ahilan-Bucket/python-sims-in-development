import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, HTMLWriter
from matplotlib.patches import Circle
import webbrowser, os


# Define Parameters

y0 = 10 # Start at height above origin
v0 = 0 # m/s down
radius = 0.5 #m , Radius of Ball
a = -9.81 #m/s^2, Gravity on Earth
r = 0.7 # coefficient of restitution, 
# r=1: perfectly elastic → no energy loss → bounces forever.
# r<1: inelastic → each bounce smaller.
# r=0: perfectly inelastic → the ball sticks to the ground

wall_below = 0 #m 

T = 60 # s, total time
dt = 0.02 # s, Sampling time

N = int(T/dt) # How many times the Function should run

collision = False
ys = []
for _ in range(N):
    if (y0-radius) <= wall_below:
        collision = True
        v0 = -r * v0

    if (collision) and (v0 == 0) :
        break

# Seems like i have to update the velocities as well because otherwise, 
# there just seems to a constant Velocity kinda of effect from acceleration

    ds = v0*dt + (1/2)*a*(dt)**2
    y0 = y0 + ds

    vf = v0 + a*dt
    v0 = vf
    ys.append(y0)

#formated = [f"{i:.2f}" for i in xs]

#rint(formated)


# ---------------- Figure: draw a real ball ----------------
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xlim(-5 ,5)
ax.set_ylim(0,20)
ax.set_aspect('equal', adjustable='box') #  if not equal, the Output is sqiushed

# Labels 
ax.set_xlabel("x (m)")
ax.set_yticks([])
ax.set_title("Ball bouncing Under Gravity")




ball = Circle(xy=[0,ys[0]],radius=radius)
ax.add_patch(ball)

def init():
    ball.center = (0,ys[0])
    return (ball,)
#Your animation callbacks must return an 
# iterable of artists when blit=True. So not (ball,), not ball

def update(i):
    ball.center = (0,ys[i])
    return (ball,)

ani = FuncAnimation(fig, update, init_func=init, frames=len(ys), interval=1000*dt, blit=True)

# --- Save to standalone HTML and show ---
outfile = "ball_fall_gravity.html"
ani.save(outfile, writer=HTMLWriter(fps=int(round(1/dt)), embed_frames=True))
plt.close(fig)
print("Saved:", os.path.abspath(outfile))
webbrowser.open("file://" + os.path.abspath(outfile))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, HTMLWriter
from matplotlib.patches import Circle
import webbrowser, os


# Define Parameters

x0 = 0 # Start at Origin
v0 = 5 # m/s right
radius = 0.5 #m , Radius of Ball

wall = 10 #m. 
wall2 = -2 # m

T = 60 # s, total time
dt = 0.02 # s, Sampling time

N = int(T/dt) # How many times the Function should run

xs = []
for _ in range(N):
    if ((x0+radius) >= wall) or ((x0-radius) <= wall2):
        v0 = -v0
    
    x0 = x0 + v0*dt
    xs.append(x0)

#formated = [f"{i:.2f}" for i in xs]

#rint(formated)


# ---------------- Figure: draw a real ball ----------------
fig, ax = plt.subplots(figsize=(8, 2.6))
ax.set_xlim(-5,11)
ax.set_ylim(0,10)

# Draw Walls
ax.axvline(x=wall)
ax.axvline(x=wall2)

ball = Circle(xy=[xs[0],radius],radius=radius)
ax.add_patch(ball)

def init():
    ball.center = (xs[0],radius)
    return (ball,)
#Your animation callbacks must return an 
# iterable of artists when blit=True. So not (ball,), not ball

def update(i):
    ball.center = (xs[i],radius)
    return (ball,)

ani = FuncAnimation(fig, update, init_func=init, frames=N, interval=1000*dt, blit=True)

# --- Save to standalone HTML and show ---
outfile = "ball_bounce_wall.html"
ani.save(outfile, writer=HTMLWriter(fps=int(round(1/dt)), embed_frames=True))
plt.close(fig)
print("Saved:", os.path.abspath(outfile))
webbrowser.open("file://" + os.path.abspath(outfile))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, HTMLWriter
from matplotlib.patches import Circle
import webbrowser, os
import random
from itertools import combinations


wall_below = 0 #m 
wall_right = 25 #m
wall_left = 0 # Making Walls below 0, makes this negative. This makes the code un-necessarliy longer
wall_above = 20
# Defining such that the entier system is on the first quadraant has the same Math as a Larger Space

T = 40 # s, total time
dt = 0.01 # s, Sampling time

N = int(T/dt) # How many times the Function should run

# Modify drop_ball to include two ball collisions
# Lets start with 1D Collison in y

# Lets Add Horizontal Velocty and Make this Guys 2D and add walls on side
nos = 2 # what if i only give it 1?
g = -9.81

balls = []
for i in range(nos):
    id = {
        'x': [random.randint(wall_left,wall_right)], 
        'y': [random.randint(wall_below,30)],
        'vx': [random.randint(-20,20)], 
        'vy': [random.randint(-20,20)],
        'radius': random.uniform(0.5,2),           # radius
        'r': 0.9,            # coefficient of restitution
        # r=1: perfectly elastic → no energy loss → bounces forever.
        # r<1: inelastic → each bounce smaller.
        # r=0: perfectly inelastic → the ball sticks to the ground
        'color': tuple(np.random.rand(3)),
    }
    balls.append(id)


def collision(x0,v0,r,radius,wall):
    # Reverse Velcoity
    v0 = -r * v0

    # Brut Force Guilty, Correct Overshoot of position
    overshoot =  wall - (x0-radius)
    x0 = x0 + overshoot

    return v0,x0

def collision2(x0,v0,r,radius,wall):
    v0 = -r * v0

    overshoot = (x0 + radius) - wall
    x0 = x0 - overshoot 
    x0 = x0
    return v0, x0


for i in range(N):
    for id in range(nos):
        # Ball A
        y_prev = balls[id]['y'][-1]
        vy_prev = balls[id]['vy'][-1]
        

        ds = vy_prev*dt + (1/2)*g*(dt)**2
        y_now = y_prev + ds
        vy_now = vy_prev + g*dt
        balls[id]['vy'].append(vy_now)
        balls[id]['y'].append(y_now)

        
        x_prev = balls[id]['x'][-1]
        vx_prev = balls[id]['vx'][-1]

        x_now = x_prev + vx_prev*dt

        vx_now = vx_prev

        balls[id]['x'].append(x_now)
        balls[id]['vx'].append(vx_now)


    # if I ask the radius positins to be eqaul, that never happenes becase in the
    # Time step, when the quantaties are a float, they will never exactly equal
    # Each other

    # Lets Shorted this logic using net distance < Combined Radius
    # Then Lets Make this a 2D collision


        radius = balls[id]['radius']
        rcoef = balls[id]['r']


        # Add Ball to ball Collisions:
        # Using Simple i-j collsion check, inspired by lukepolson

        ids = np.arange(0,len(balls),1)

        ids_pairs = np.asarray(list(combinations(ids,2)))

        # [[0, 1],[0, 2],...]

        for pair in ids_pairs:
            id_1 = pair[0]
            id_2 = pair[1]

            id_1y = balls[id_1]['y'][-1]
            id_2y = balls[id_2]['y'][-1]

            
            id_1x = balls[id_1]['x'][-1]
            id_2x = balls[id_2]['x'][-1]
            
            net_distance = np.sqrt((id_1y - id_2y)**2 + (id_1x - id_2x)**2)

            id_1_radius = balls[id_1]['radius']
            id_2_radius = balls[id_2]['radius']

            id_1_vy = balls[id_1]['vy'][-1]
            id_2_vy = balls[id_2]['vy'][-1]

            id_1_rcoef = balls[id_1]['r']
            id_2_rcoef = balls[id_2]['r']

            '''
            id_1_vx = balls[id_1]['vx'][-1]
            id_2_vx = balls[id_2]['vx'][-1]
            '''

            if net_distance <= id_1_radius+id_2_radius:
                id_1_vy = -id_1_rcoef*id_1_vy
                id_2_vy = -id_2_rcoef*id_2_vy

                # Append
                balls[id_1]['vy'][-1] = id_1_vy 
                balls[id_2]['vy'][-1] = id_2_vy



                '''
                x coef
                '''


            

        """    
        net_distance = np.sqrt((A_y0 - B_y0)**2 + (A_x0 - B_x0)**2)
        if net_distance < A_radius+B_radius:
            A_vy0 = -A_r * A_vy0
            B_vy0 = -B_r * B_vy0

            A_vx0 = -A_r * A_vx0
            B_vx0 = -B_r * B_vx0
        """


        # Alter the Current 'now' states if there is a collision

        '''
        Right after a bounce, the old state is still exactly on the wall (e.g., floor contact point).

        If you check that old state, you think “still colliding” and flip again—even though the ball is already moving away.

        You keep flipping every frame → looks stuck/jittery.
        '''
        # Balls Stick now for some reason
        if (y_now-radius) <= wall_below:
            vy_now,y_now = collision(y_now,vy_now,rcoef,radius,wall=wall_below)

        if (y_now+radius) >= wall_above:
            vy_now,y_now = collision2(y_now,vy_now,rcoef,radius,wall=wall_above)

        if (x_now - radius) <= wall_left:
            vx_now, x_now = collision(x_now,vx_now,rcoef,radius,wall=wall_left)

        if (x_now + radius) >= wall_right:
            vx_now, x_now = collision2(x_now,vx_now,rcoef,radius,wall=wall_right)

        
        balls[id]['y'][-1]  = y_now
        balls[id]['vy'][-1] = vy_now
        balls[id]['x'][-1]  = x_now
        balls[id]['vx'][-1] = vx_now

# ---------------- Figure: draw a real ball ----------------
fig, ax = plt.subplots(figsize=(11, 8))
ax.set_xlim(wall_left-5 ,wall_right+5)
#ax.set_xlim(-10 ,30)

ax.set_ylim(0,20)
ax.set_aspect('equal', adjustable='box') #  if not equal, the Output is sqiushed


# ALl walls
ax.axvline(x=wall_left,color="Black",lw=3)
ax.axvline(x=wall_right,color="Black",lw=3)
ax.axhline(y=wall_above,color="Black",lw=3)


# Labels 
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")

#ax.set_yticks([])
ax.set_title(f"{nos} Balls bouncing Under Gravity")


# List of Circles:

patches = []
for id in range(nos):
    ball = Circle(xy=[balls[id]['x'][0],balls[id]['y'][0]],
                                  radius=balls[id]['radius'],color=balls[id]['color'])
    ax.add_patch(ball)
    patches.append(ball)
    

def init():
    for id,patch in enumerate(patches):
        patch.center = (balls[id]['x'][0],balls[id]['y'][0])
    return tuple(patches) # Must be Iterable with Artists
#Your animation callbacks must return an 
# iterable of artists when blit=True. So not (ball,), not ball

def update(i):
    for id,patch in enumerate(patches):
        patch.center = (balls[id]['x'][i],balls[id]['y'][i])
    return tuple(patches) # Must be Iterable with Artists

ani = FuncAnimation(fig, update, init_func=init, frames=N, interval=1000*dt, blit=True)

# --- Save to standalone HTML and show ---
outfile = "Nballs_fall_gravity_collision.html"
ani.save(outfile, writer=HTMLWriter(fps=int(round(1/dt)), embed_frames=True))
plt.close(fig)
print("Saved:", os.path.abspath(outfile))
webbrowser.open("file://" + os.path.abspath(outfile))
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import mediapipe as mp
import cv2
import math
import time

# psi_solve2/app.py
# ==========================================
# 0. STREAMLIT SESSION STATE SETUP
# ==========================================
if 'countdown_finished' not in st.session_state:
    st.session_state.countdown_finished = False
if 'V_user_defined' not in st.session_state:
    st.session_state.V_user_defined = None

# ==========================================
# 1. PHYSICS ENGINE (Unchanged)
# ==========================================
hbar = 1
m = 1
L = 50
N_GRID = 2000

def make_grid(L=L, N=N_GRID):
    """Returns x_full (N+2), dx, x_solver (N)"""
    x = np.linspace(-L/2, L/2, N+2)
    dx = x[1] - x[0]
    x_internal = x[1:-1]
    return x, dx, x_internal

def kinetic_operator(N, dx, hbar=hbar, m=m):
    """Builds the Kinetic Energy Matrix (N x N)"""
    main_diagonal = (1/dx**2) * np.diag(-2 * np.ones(N))
    off_diagonal1 = (1/dx**2) * np.diag(np.ones(N-1), -1)
    off_diagonal2 = (1/dx**2) * np.diag(np.ones(N-1), 1)
    D2 = (main_diagonal + off_diagonal1 + off_diagonal2)
    T = (-(hbar**2 / (2*m)) * D2)
    return T

def solve(T, V_full, dx):
    """Solves H psi = E psi"""
    V_internal = V_full[1:-1]
    H = T + np.diag(V_internal)
    E, psi = np.linalg.eigh(H)
    
    # Normalize
    for i in range(psi.shape[1]):
        norm_factor = np.sum(psi[:, i]**2) * dx
        psi[:, i] = psi[:, i] / np.sqrt(norm_factor)
        
    return E, psi

def plot_dead(E, psi, V_full, x_full, nos=5):
    """Returns a Matplotlib Figure object for Streamlit to display."""
    plt.style.use("dark_background")
    fig, (ax_main, ax_bar) = plt.subplots(
        1, 2, figsize=(10, 6), gridspec_kw={'width_ratios': [5, 1]}
    )
    plt.subplots_adjust(bottom=0.25, wspace=0.4)

    # Ensure nos doesn't exceed available states
    states_to_plot = min(nos, len(E))

    # Use first two available states to determine scale
    if len(E) >= 2:
        scale = (E[1] - E[0]) * 0.4
    elif len(E) == 1:
        scale = E[0] * 0.1 # Fallback
    else:
        scale = 1.0

    x_solver = x_full[1:-1]
    
    max_E = E[states_to_plot-1] if states_to_plot > 0 else 10.0
    window_height = max_E * 1.5

    # 1. Plot Wavefunctions
    for n in range(states_to_plot):
        psi_nth = psi[:, n]
        psi_nth = psi_nth / np.max(np.abs(psi_nth))
        vertical_shifted = psi_nth * scale + E[n]
        
        V_internal = V_full[1:-1]
        vertical_shifted[V_internal > 1e5] = np.nan

        color = plt.colormaps['tab20'].colors[n % 20]
        ax_main.plot(x_solver, vertical_shifted, 
                     label=f"n={n+1}", lw=1.2, color=color)

    # 2. Plot Potential
    V_clipped = np.clip(V_full, 0, window_height)
    ax_main.plot(x_full, V_clipped, color='white', lw=2, label="V(x)")

    # Formatting
    ax_main.set_xlabel("x [a.u.]")
    ax_main.set_ylabel("Energy")
    ax_main.set_ylim(0, window_height * 1.1)
    ax_main.set_title("Eigenstates")
    ax_main.legend(loc="upper right", fontsize=8)

    # Bar Chart
    ax_bar.set_title("Spectrum")
    ax_bar.set_xticks([])
    ax_bar.set_ylim(0, np.max(E[:states_to_plot]) * 1.1)
    for n in range(states_to_plot):
        ax_bar.axhline(E[n], c=plt.colormaps['tab20'].colors[n % 20], lw=1)

    return fig

# ==========================================
# 2. COMPUTER VISION (MediaPipe) (Unchanged)
# ==========================================
def process_frame_to_potential(frame):
    """Processes frame to get potential."""
    mp_hands = mp.solutions.hands
    
    with mp_hands.Hands(max_num_hands=2, min_detection_confidence=0.5) as hands:
        h, w, _ = frame.shape
        rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        res = hands.process(rgb)

        if not res.multi_hand_landmarks:
            return None, "No Hands Detected"

        # --- LOGIC: Square Well vs QHO ---
        
        # 1. Square Well (2 Hands)
        if len(res.multi_hand_landmarks) >= 2:
            INDEX_TIP_ID = 8
            x_coords = [lm.landmark[INDEX_TIP_ID].x * w for lm in res.multi_hand_landmarks]
            x_coords.sort()
            
            xL_hand, xR_hand = x_coords[0], x_coords[1]
            well_width = xR_hand - xL_hand
            
            center_screen = w / 2
            centered_L = center_screen - (well_width / 2)
            centered_R = center_screen + (well_width / 2)
            
            x_space = np.linspace(0, w, 400)
            pot_profile = np.ones_like(x_space)
            pot_profile[(x_space > centered_L) & (x_space < centered_R)] = 0
            
            return pot_profile, "Square Well (Captured)"

        # 2. Harmonic Oscillator (1 Hand)
        elif len(res.multi_hand_landmarks) == 1:
            lm = res.multi_hand_landmarks[0]
            THUMB = lm.landmark[4]
            INDEX = lm.landmark[8]
            
            dx = INDEX.x - THUMB.x
            dy = INDEX.y - THUMB.y
            dist = math.sqrt(dx**2 + dy**2)
            
            A = np.interp(dist, [0.05, 0.3], [100.0, 1.0]) 
            
            x_space = np.linspace(-1, 1, 400)
            pot_profile = A * (x_space**2)
            
            pot_profile = np.clip(pot_profile, 0, 100)
            pot_profile = pot_profile / 100.0
            
            return pot_profile, f"Harmonic Oscillator (k={A:.1f})"
            
    return None, "Error"

# ==========================================
# 3. STATIC POTENTIALS
# ==========================================
def generate_static_potential(L, N_solver, potential_type):
    """Generates a fixed potential for comparison."""
    x_full, dx, x_solver = make_grid(L, N_GRID)
    V_physics = np.zeros(N_solver)

    if potential_type == "Static Square Well (Wide)":
        # Square well from -10 to 10
        x_internal = x_solver
        V_physics[np.abs(x_internal) > 10] = 200.0 
    
    elif potential_type == "Static Harmonic Oscillator":
        # QHO V(x) = 0.5 * k * x^2
        x_internal = x_solver
        k = 0.5  # Fixed spring constant
        V_physics = 0.5 * k * x_internal**2
        V_physics = V_physics / np.max(V_physics) * 200.0
        
    # Pad with walls
    V_full = np.pad(V_physics, (1, 1), 'constant', constant_values=1e10)
    return V_full, "Static Potential"

# ==========================================
# 4. STREAMLIT APP LAYOUT & LOGIC
# ==========================================
st.set_page_config(page_title="Quantum Hand Solver", layout="wide")

st.title("Quantum Potential Solver ⚛️")
st.write("Use your hands or static settings to explore the time-independent Schrödinger equation.")

col1, col2 = st.columns([1, 2])

# --- CONTROL PANEL ---
with col1:
    st.header("Control Panel")
    
    # User control for number of eigenstates
    nos_user = st.slider("Number of Eigenstates (n)", 1, 10, 5)

    # Dropdown for potential selection
    potential_mode = st.selectbox(
        "Select Potential Mode",
        ("Static Square Well (Wide)", 
         "Static Harmonic Oscillator", 
         "Hand Gesture (Camera)")
    )
    
    # Reset countdown state when mode changes
    if potential_mode != "Hand Gesture (Camera)":
        st.session_state.countdown_finished = False

    # Initialize variables for the main loop
    V_full_to_solve = None
    status_msg = ""
    run_analysis = False
    
    # --- CAMERA/HAND GESTURE LOGIC ---
    if potential_mode == "Hand Gesture (Camera)":
        
        st.subheader("Hand Gesture Controls")
        st.info("Instructions:\n1. Click **'Start Countdown'**.\n2. Get your **two hands** ready for the Square Well, or **one hand pinch** for the QHO.\n3. Click **'Take a snapshot'** when the countdown finishes.")
        
        # 1. Start Countdown Button
        if st.button("Start Countdown"):
            st.session_state.countdown_finished = False
            
            # Display the countdown and delay
            countdown_placeholder = st.empty()
            
            for i in range(10, 0, -1):
                countdown_placeholder.markdown(f"## 📸 Capture in **{i}** seconds...")
                time.sleep(1)
            
            countdown_placeholder.success("Ready! Hold your pose and click 'Take a snapshot'.")
            st.session_state.countdown_finished = True
            st.rerun() # Rerun to show the success message clearly

        # 2. Camera Input (conditionally visible)
        img_file = None
        if st.session_state.countdown_finished:
            img_file = st.camera_input("Take a snapshot", key="camera_capture")
        else:
            # Show camera for initial setup visibility
            st.camera_input("Take a snapshot", key="camera_initial_setup")

        # 3. Process captured image
        if img_file is not None and st.session_state.countdown_finished:
            
            # Convert uploaded file to OpenCV format
            file_bytes = np.asarray(bytearray(img_file.read()), dtype=np.uint8)
            frame = cv2.imdecode(file_bytes, 1)

            # FIX: Mirror the image for correct x-axis orientation
            frame = cv2.flip(frame, 1) 
            
            # Process Image
            V_vision, msg = process_frame_to_potential(frame)
            
            if V_vision is not None:
                st.success(f"Detected: {msg}")
                # Store V_vision for the solver
                st.session_state.V_user_defined = V_vision
                run_analysis = True
            else:
                st.error(msg)
                st.session_state.countdown_finished = False # Reset for new attempt
                st.rerun()


    else:  # Static potentials
        if potential_mode == "Static Square Well (Wide)":
            st.subheader("Infinite Square Well Settings")

            # Half-width slider (well extends from -W to +W)
            W = st.slider("Well Half-Width (|x| < W)", 2.0, 20.0, 10.0, step=0.5)

            if st.button("Solve Infinite Well"):
                x_full, dx, x_solver = make_grid(L, N_GRID)
                V_physics = np.zeros_like(x_solver)
                V_physics[np.abs(x_solver) > W] = 200.0  # High walls
                V_full_to_solve = np.pad(V_physics, (1,1), constant_values=1e10)
                run_analysis = True
                status_msg = f"Infinite Square Well (|x| < {W})"

        elif potential_mode == "Static Harmonic Oscillator":
            st.subheader("Harmonic Oscillator Settings")

            # k controls curvature / spacing between energy levels
            k = st.slider("Spring Constant (k)", 0.1, 5.0, 0.5, step=0.1)

            if st.button("Solve Harmonic"):
                x_full, dx, x_solver = make_grid(L, N_GRID)
                V_physics = 0.5 * k * x_solver**2
                V_physics = V_physics / np.max(V_physics) * 200  # scale for plot
                V_full_to_solve = np.pad(V_physics, (1,1), constant_values=1e10)
                run_analysis = True
                status_msg = f"Harmonic Oscillator (k = {k:.2f})"

# --- ANALYSIS AND PLOTTING ---
with col2:
    if run_analysis:
        
        # Determine which potential to solve
        if potential_mode == "Hand Gesture (Camera)":
            # Map the stored V_vision to physics grid
            V_raw = st.session_state.V_user_defined
            if V_raw is not None:
                x_full, dx, x_solver = make_grid(L, N_GRID)
                V_interpolated = np.interp(
                    np.linspace(0, 1, len(x_solver)), 
                    np.linspace(0, 1, len(V_raw)), 
                    V_raw
                )
                V_physics = V_interpolated * 200.0
                V_full_to_solve = np.pad(V_physics, (1, 1), 'constant', constant_values=1e10)
            
        
        if V_full_to_solve is not None:
            # 4. Compute
            with st.spinner(f"Solving Schrödinger Equation for {potential_mode}..."):
                T = kinetic_operator(N_GRID, dx)
                E, psi = solve(T, V_full_to_solve, dx)
                
            # 5. Plot
            st.write("### Quantum States")
            fig = plot_dead(E, psi, V_full_to_solve, x_full, nos=nos_user)
            st.pyplot(fig)

            # Reset state for next run
            if potential_mode == "Hand Gesture (Camera)":
                st.session_state.countdown_finished = False
                st.session_state.V_user_defined = None
        
        elif potential_mode == "Hand Gesture (Camera)":
            st.error("No potential data captured. Please try the snapshot sequence again.")
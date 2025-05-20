# Graph comparisons
import numpy as np
import matplotlib.pyplot as plt
import os # Added to check for file existence

data_filename = 'trajectory_output.dat'
ref_x_filename = 'X_coords_ref.npy' # Filename for reference X coordinates
ref_y_filename = 'Y_coords_ref.npy' # Filename for reference Y coordinates

params = {}
try:
    with open(data_filename, 'r') as f:
        for line in f:
            if not line.startswith('#'): break
            if line.startswith('# Trajectory Parameters:'): continue
            if line.startswith('# K') and '\t' in line: continue
            cleaned_line = line[1:].strip()
            parts = cleaned_line.split('=', 1)
            if len(parts) == 2:
                key = parts[0].strip()
                val_str = parts[1].strip()
                try:
                    # Try to convert to float, but keep as string if it fails
                    if val_str.lower() == 'nan': # Handle 'nan' string explicitly
                        params[key] = np.nan
                    else:
                        params[key] = float(val_str)
                except ValueError:
                    params[key] = val_str
except FileNotFoundError:
    print(f"Error: Data file '{data_filename}' not found.")
    exit()
except Exception as e: # Catch other potential errors during parameter parsing
    print(f"Error reading parameters from '{data_filename}': {e}")
    # Decide if you want to exit or continue; for now, let's try to continue

try:
    data = np.loadtxt(data_filename, comments='#')
    if data.ndim == 1 and data.shape[0] >= 5: data = data.reshape(1,data.shape[0])
    if data.shape[0] == 0 or data.shape[1] < 5: raise ValueError('No data or not enough columns in main data file.')
except Exception as e:
    print(f"Error loading data from '{data_filename}': {e}. Ensure it's not empty and has 5 columns.")
    exit()

K, r, phi, x, y = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]

# --- Load the additional X and Y coordinates for the reference trajectory ---
X_coords_ref = None
Y_coords_ref = None
plot_reference_trajectory = False

if os.path.exists(ref_x_filename) and os.path.exists(ref_y_filename):
    try:
        X_coords_ref = np.load(ref_x_filename)
        Y_coords_ref = np.load(ref_y_filename)
        if X_coords_ref.shape == Y_coords_ref.shape and X_coords_ref.ndim == 1 and X_coords_ref.size > 0:
            plot_reference_trajectory = True
            print(f"Successfully loaded reference trajectory data from '{ref_x_filename}' and '{ref_y_filename}'.")
        else:
            print(f"Warning: Reference trajectory files loaded, but shapes mismatch or not 1D/empty. X_shape: {X_coords_ref.shape}, Y_shape: {Y_coords_ref.shape}")
    except Exception as e:
        print(f"Error loading reference trajectory data: {e}")
else:
    print(f"Warning: Reference trajectory files ('{ref_x_filename}' or '{ref_y_filename}') not found. Plotting only primary trajectory.")


# --- Parameters for Title (handle potential NaN from params) ---
M_val = params.get('M', 1.0) # Default to 1.0 if not found
b_val_str = str(params.get('b_calculated', 'N/A')).lower() # Use 'N/A' for missing
b_val = np.nan if b_val_str == 'nan' or b_val_str == 'n/a' else float(b_val_str)

psi_val_rad_str = str(params.get('psi_rad_input', 'N/A')).lower()
psi_val_rad = np.nan if psi_val_rad_str == 'nan' or psi_val_rad_str == 'n/a' else float(psi_val_rad_str)
psi_val_deg = psi_val_rad * 180/np.pi if not np.isnan(psi_val_rad) else np.nan


fig, axs = plt.subplots(1, 2, figsize=(18, 7)) # Adjusted figure size

# --- x-y plot ---
axs[0].plot(x, y, label='Traj 1 (C code) x-y path', color='dodgerblue', lw=1.5) # Main trajectory

if plot_reference_trajectory:
    axs[0].plot(X_coords_ref, Y_coords_ref, label='Traj 2 (Reference) x-y path', color='orangered', linestyle='--', lw=1.5) # Additional trajectory

axs[0].plot(0,0,'ko', markersize=5, label='BH (r=0)')
event_horizon_radius_key = 'SchwarzschildBlackHoleRadius' # Key used in your plot_trajectory_data.c
event_horizon_radius = params.get(event_horizon_radius_key, 2 * M_val) # Use M_val if key not found
circle = plt.Circle((0, 0), event_horizon_radius, color='dimgray', alpha=0.4, fill=True, linestyle='--', label=f'Event Horizon (r={event_horizon_radius:.2f}M)')
axs[0].add_artist(circle)
axs[0].set_xlabel('x / M', fontsize=12)
axs[0].set_ylabel('y / M', fontsize=12)
axs[0].set_aspect('equal', 'box')
axs[0].legend(fontsize='small')
axs[0].grid(True, linestyle=':', alpha=0.7)
axs[0].tick_params(axis='both', which='major', labelsize=10)


# --- r(K) and phi(K) plot (for the primary trajectory from trajectory_output.dat) ---
axs[1].plot(K, r, label='r(K) (Traj 1)', color='blue', lw=1.5)
axs[1].set_xlabel('Affine Parameter K', fontsize=12)
axs[1].set_ylabel('r / M (Traj 1)', color='blue', fontsize=12)
axs[1].tick_params(axis='y', labelcolor='blue', labelsize=10)
axs[1].tick_params(axis='x', labelsize=10)

ax_phi = axs[1].twinx()
# Unwrap phi only if it's not empty or all NaNs
phi_to_unwrap = phi[~np.isnan(phi)] # Ensure we only unwrap valid numbers
if phi_to_unwrap.size > 0 :
    phi_unwrapped = np.unwrap(phi_to_unwrap)
    # We need to plot against K values corresponding to non-NaN phi
    K_for_phi = K[~np.isnan(phi)]
    if K_for_phi.size == phi_unwrapped.size : # Ensure sizes match after potential NaN removal
        ax_phi.plot(K_for_phi, phi_unwrapped, label='$\phi_{unwrapped}(K)$ (Traj 1)', color='green', linestyle=':', lw=1.5)
    else:
        print("Warning: Mismatch in K and unwrapped phi sizes, skipping phi plot for Traj 1.")
else:
    print("Warning: Phi data for Traj 1 is empty or all NaN, skipping phi plot.")

ax_phi.set_ylabel('Unwrapped $\phi$ (rad) (Traj 1)', color='green', fontsize=12)
ax_phi.tick_params(axis='y', labelcolor='green', labelsize=10)

h1, l1 = axs[1].get_legend_handles_labels()
h2, l2 = ax_phi.get_legend_handles_labels()
axs[1].legend(h1+h2, l1+l2, loc='best', fontsize='small')
axs[1].grid(True, linestyle=':', alpha=0.7)

# --- Super Title ---
psi_deg_str = f"{psi_val_deg:.1f}$^\\circ$" if not np.isnan(psi_val_deg) else "N/A"
b_str = f"{b_val:.3f}" if not np.isnan(b_val) else "N/A"
title_str = f'Photon Trajectory Comparison (Traj 1 from C): M={M_val:.2f}, b={b_str}, $\\psi$={psi_val_rad:.3f} rad ({psi_deg_str})'
fig.suptitle(title_str, fontsize=14)

plt.subplots_adjust(top=0.90, wspace=0.35) # Adjusted top for suptitle and wspace
output_filename = 'trajectory_comparison_plot.png'
plt.savefig(output_filename, dpi=150)
print(f"Python comparison plot saved to {output_filename}")
# plt.show() # Uncomment to display plot interactively
import numpy as np
import matplotlib.pyplot as plt

data_filename = 'trajectory_output.dat'

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
                    params[key] = float(val_str)
                except ValueError:
                    params[key] = val_str 
except FileNotFoundError:
    print(f"Error: Data file '{data_filename}' not found.")
    exit()

try:
    data = np.loadtxt(data_filename, comments='#')
    if data.ndim == 1 and data.shape[0] >= 5: data = data.reshape(1,data.shape[0]) 
    if data.shape[0] == 0 or data.shape[1] < 5: raise ValueError('No data or not enough columns.')
except Exception as e: 
    print(f"Error loading data from '{data_filename}': {e}. Ensure it's not empty and has 5 columns.")
    exit()

K, r, phi, x, y = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]

M_val = params.get('M', 1.000000000000000e+00)
b_val = params.get('b_calculated', np.nan if str(params.get('b_calculated', 'nan')).lower() == 'nan' else float(params.get('b_calculated', np.nan)))
psi_val_rad = params.get('psi_rad_input', -1.570796326794897e+00)

fig, axs = plt.subplots(1, 2, figsize=(15, 6))

# x-y plot
axs[0].plot(x, y, label='x-y path', color='dodgerblue', lw=1.5)
axs[0].plot(0,0,'ko', markersize=5, label='BH (r=0)') 
event_horizon_radius = params.get('SchwarzschildBlackHoleRadius', 2 * M_val)
circle = plt.Circle((0, 0), event_horizon_radius, color='dimgray', alpha=0.4, fill=True, linestyle='--', label=f'Event Horizon (r={event_horizon_radius:.2f}M)')
axs[0].add_artist(circle)
axs[0].set_xlabel('x / M')
axs[0].set_ylabel('y / M')
axs[0].set_aspect('equal', 'box')
axs[0].legend(fontsize='small')
axs[0].grid(True, linestyle=':', alpha=0.7)
axs[0].set_ylim(-5,5)

# r(K) and phi(K) plot
axs[1].plot(K, r, label='r(K)', color='blue', lw=1.5)
axs[1].set_xlabel('Affine Parameter K')
axs[1].set_ylabel('r / M', color='blue')
axs[1].tick_params(axis='y', labelcolor='blue')
ax_phi = axs[1].twinx()
phi_unwrapped = np.unwrap(phi) 
ax_phi.plot(K, phi_unwrapped, label='$\phi_{unwrapped}(K)$', color='green', linestyle=':', lw=1.5)
ax_phi.set_ylabel('Unwrapped $\phi$ (rad)', color='green')
ax_phi.tick_params(axis='y', labelcolor='green')
h1, l1 = axs[1].get_legend_handles_labels()
h2, l2 = ax_phi.get_legend_handles_labels()
axs[1].legend(h1+h2, l1+l2, loc='best', fontsize='small')
axs[1].grid(True, linestyle=':', alpha=0.7)

title_str = f'Photon Trajectory: M={M_val:.2f}, b={b_val:.3f}, $\\psi$={psi_val_rad:.3f} rad ({psi_val_rad*180/np.pi:.1f}$^\\circ$)'
fig.suptitle(title_str, fontsize=14)
plt.subplots_adjust(top=0.88, wspace=0.3)
plt.savefig('trajectory_plot_py.png', dpi=150)
print(f"Python plot saved to trajectory_plot_py.png")
plt.show() # Uncomment to display plot interactively
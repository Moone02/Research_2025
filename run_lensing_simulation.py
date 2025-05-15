# run_lensing_simulation.py

import numpy as np
# Import the functions from your saved file
from Python_Schwarzschild_Code_Commented import (
    Results_Radial,  # Or Results_Cartesian, Results_Radial_Light_Ring
    map_photons,
    # You might also import compute_trajectory or Image_map for direct testing
)

# --- Step 1: Generate Photon Mapping Data ---
print("--- Starting Photon Data Generation ---")
observer_x = 100.0
window_plane_x = 20.0
source_plane_x = -50.0
max_radial_on_window = 15.0
num_radial_steps = 50 # Keep this small for a quick test!
chunk_s = 10000       # Smaller chunk size for faster feedback during testing

# Call one of the Results_... functions
# This will create .npy files in the current directory
photon_data_files = Results_Radial(
    x_0=source_plane_x,
    x_1=window_plane_x,
    x_2=observer_x,
    R_max=max_radial_on_window,
    n=num_radial_steps,
    chunk_size=chunk_s
)
print(f"Photon data generation complete. Files: {photon_data_files}")
print("--- Photon Data Generation Finished ---")

# --- Step 2: Prepare Your Source Image ---
print("\n--- Preparing Source Image ---")
def create_simple_grid_image(height, width, grid_size=20):
    img = np.full((height, width, 3), 200, dtype=np.uint8) # Light gray background
    for r in range(0, height, grid_size):
        img[r, :, :] = 50 # Dark gray horizontal lines
    for c in range(0, width, grid_size):
        img[:, c, :] = 50 # Dark gray vertical lines
    return img

source_image_array = create_simple_grid_image(height=100, width=150)
# Or load from an actual image file:
# from PIL import Image
# source_image_path = "your_source_image.png"
# try:
#     img = Image.open(source_image_path).convert('RGB')
#     source_image_array = np.array(img)
# except FileNotFoundError:
#     print(f"ERROR: Source image {source_image_path} not found!")
#     exit()
print("Source image prepared.")
print("--- Source Image Preparation Finished ---")


# --- Step 3: Generate the Lensed Image ---
print("\n--- Starting Image Mapping ---")
output_image_filename = "lensed_grid_image.png"
output_pixel_shape = (300, 300) # (height, width) for the final lensed image

# Check if photon_data_files were actually created
if not photon_data_files or not all(isinstance(f, str) for f in photon_data_files):
    print("Error: No valid photon data files were generated. Cannot proceed with mapping.")
else:
    success = map_photons(
        image_source=source_image_array,
        photon_chunk_files=photon_data_files,
        save_path=output_image_filename,
        output_shape=output_pixel_shape, # Using fixed output shape mode
        default_color=(0, 0, 0),         # Black background
        flip_z_axis_render=True
    )

    if success:
        print(f"Lensed image successfully saved to {output_image_filename}")
    else:
        print(f"Image mapping failed for {output_image_filename}")
print("--- Image Mapping Finished ---")
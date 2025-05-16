//main_test.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> // For NAN if not in GSL headers directly for C89/C90
#include "schwarzschild_tracer.h"

// Helper to create a simple dummy PPM image for testing map_photons
PPMImage* create_dummy_ppm(int width, int height) {
    PPMImage* img = malloc(sizeof(PPMImage));
    if (!img) return NULL;
    img->width = width;
    img->height = height;
    img->channels = 3;
    img->data = malloc((size_t)width * height * 3);
    if (!img->data) {
        free(img);
        return NULL;
    }
    for (int r = 0; r < height; ++r) {
        for (int c = 0; c < width; ++c) {
            size_t idx = ((size_t)r * width + c) * 3;
            img->data[idx + 0] = (unsigned char)((c * 255) / width); // Red gradient
            img->data[idx + 1] = (unsigned char)((r * 255) / height); // Green gradient
            img->data[idx + 2] = 128;                         // Blue constant
        }
    }
    return img;
}


int main() {
    printf("--- Schwarzschild Tracer Test Suite ---\n\n");
    gsl_set_error_handler_off(); // Use return codes instead of aborting

    // --- Test 1: compute_trajectory ---
    printf("--- Test 1: compute_trajectory (scattering-like) ---\n");
    // Parameters: r0=20, phi0=0, M=1, psi=-1.25 (inward, should scatter), r_max=100
    // x_stop=Not active, x_targets=NULL, t_end=large, num_interp=10
    PhotonTrajectory *traj1 = compute_trajectory(20.0, 0.0, 1.0, -1.25, 100.0,
                                                NAN, false, NULL, DEFAULT_T_END, 10);
    if (traj1) {
        if (traj1->error_code == 0 && traj1->K && traj1->K->size > 0) {
            printf("  compute_trajectory SUCCESS. Trajectory has %zu points.\n", traj1->K->size);
            printf("  First point: K=%.2f, r=%.2f, phi=%.2f, x=%.2f, y=%.2f\n",
                   gsl_vector_get(traj1->K, 0), gsl_vector_get(traj1->r, 0),
                   gsl_vector_get(traj1->phi, 0), gsl_vector_get(traj1->x, 0),
                   gsl_vector_get(traj1->y, 0));
            size_t last_idx = traj1->K->size - 1;
            printf("  Last point:  K=%.2f, r=%.2f, phi=%.2f, x=%.2f, y=%.2f\n",
                   gsl_vector_get(traj1->K, last_idx), gsl_vector_get(traj1->r, last_idx),
                   gsl_vector_get(traj1->phi, last_idx), gsl_vector_get(traj1->x, last_idx),
                   gsl_vector_get(traj1->y, last_idx));
        } else {
            printf("  compute_trajectory FAILED or no points. Error code: %d\n", traj1->error_code);
        }
        free_photon_trajectory(traj1);
    } else {
        printf("  compute_trajectory FAILED (NULL result).\n");
    }
    printf("\n");

    // --- Test 2: compute_trajectory_crossings_only ---
    printf("--- Test 2: compute_trajectory_crossings_only (inward radial) ---\n");
    gsl_vector *x_targets_test2 = gsl_vector_alloc(3);
    gsl_vector_set(x_targets_test2, 0, 4.0);
    gsl_vector_set(x_targets_test2, 1, 3.0);
    gsl_vector_set(x_targets_test2, 2, 2.1); // Near horizon
    // Parameters: r0=5, phi0=0, M=1, psi=PI (radially inward), r_max=10
    // x_stop=active at x=0, t_end=large
    TrajectoryCrossings *crossings1 = compute_trajectory_crossings_only(
        5.0, 0.0, 1.0, M_PI, 10.0, 0.0, true, x_targets_test2, DEFAULT_T_END);

    if (crossings1) {
        if (crossings1->error_code == 0) {
            printf("  compute_trajectory_crossings_only SUCCESS.\n");
            for (size_t i = 0; i < crossings1->num_x_targets; ++i) {
                printf("    Target x=%.2f: %zu crossings found.\n",
                       gsl_vector_get(x_targets_test2, i),
                       crossings1->crossings_y_at_x_targets[i] ? crossings1->crossings_y_at_x_targets[i]->size : 0);
                if (crossings1->crossings_y_at_x_targets[i] && crossings1->crossings_y_at_x_targets[i]->size > 0) {
                    printf("      First y-crossing: %.4f\n", gsl_vector_get(crossings1->crossings_y_at_x_targets[i], 0));
                }
            }
        } else {
            printf("  compute_trajectory_crossings_only FAILED. Error code: %d\n", crossings1->error_code);
        }
        free_trajectory_crossings(crossings1);
    } else {
        printf("  compute_trajectory_crossings_only FAILED (NULL result).\n");
    }
    gsl_vector_free(x_targets_test2);
    printf("\n");

    // --- Test 3: image_map ---
    printf("--- Test 3: image_map (simple case) ---\n");
    // Observer at x=20, window at x=10, source at x=-5. Y=0.1, Z=0 on window.
    ImageMapResult im_res1 = image_map(0.1, 0.0, -5.0, 10.0, 20.0, 1.0);
    if (im_res1.error_code == 0) {
        if (!im_res1.miss) {
            printf("  image_map HIT: y_window=%.4f, y_image=%.4f\n", im_res1.y_window, im_res1.y_image);
        } else {
            printf("  image_map MISS.\n");
        }
    } else {
        printf("  image_map FAILED. Error code: %d\n", im_res1.error_code);
    }
    printf("\n");

    // --- Test 4: results_radial (minimal to test file creation) ---
    printf("--- Test 4: results_radial (minimal) ---\n");
    int num_radial_files = 0;
    // x0=-10, x1=5, x2=20, R_max=0.05, n=1 (few points), chunk_size=10
    char **radial_files = results_radial(-10.0, 5.0, 20.0, 0.05, 1, 10, &num_radial_files);
    if (radial_files) {
        printf("  results_radial potentially created %d files:\n", num_radial_files);
        for (int i = 0; i < num_radial_files; ++i) { // Iterate up to num_radial_files
            if (radial_files[i]) printf("    File: %s\n", radial_files[i]);
        }
        if (num_radial_files > 0 && radial_files[0]) {
            printf("  Attempting to load first chunk: %s\n", radial_files[0]);
            size_t points_loaded;
            PhotonMapDataPoint *loaded_chunk = load_photon_data_chunk(radial_files[0], &points_loaded);
            if (loaded_chunk) {
                printf("    Successfully loaded %zu points from chunk.\n", points_loaded);
                if (points_loaded > 0) {
                    printf("      First point: y_win_x=%.3f, y_win_y=%.3f, y_img_x=%.3f, y_img_y=%.3f\n",
                           loaded_chunk[0].y_window_cart_x, loaded_chunk[0].y_window_cart_y,
                           loaded_chunk[0].y_image_cart_x,  loaded_chunk[0].y_image_cart_y);
                }
                free(loaded_chunk);
            } else {
                printf("    Failed to load chunk or chunk empty.\n");
            }
        }
        // free_string_array will handle the NULL terminator if add_string_to_list added it
    } else {
        printf("  results_radial FAILED (NULL file list returned).\n");
    }
    printf("\n");


    // --- Test 5: map_photons (minimal) ---
    printf("--- Test 5: map_photons (minimal) ---\n");
    PPMImage *dummy_source_img = create_dummy_ppm(64, 64);
    if (dummy_source_img) {
        printf("  Created dummy source image (64x64).\n");
        // Use files from Test 4 if any were created
        if (radial_files && num_radial_files > 0) {
            // Define output shape for auto-bounds mode
            int output_dims[2] = {32, 32}; // height, width
            RGBColor bg_color = {10, 10, 10};

            printf("  Attempting to render map_photons with %d chunk file(s) to 'test_lensed_output.ppm'.\n", num_radial_files);
            int map_status = map_photons(
                NULL, dummy_source_img,     // Use raw data
                radial_files, num_radial_files, // Pass explicit count
                "test_lensed_output.ppm",
                NULL, 0.0,                  // Auto-bounds mode
                output_dims,
                bg_color, false);

            if (map_status == 0) {
                printf("  map_photons SUCCESS. Output image 'test_lensed_output.ppm' should be created.\n");
            } else {
                printf("  map_photons FAILED. Status: %d\n", map_status);
            }
        } else {
            printf("  Skipping map_photons test as no chunk files were generated by results_radial.\n");
        }
        free_ppm_image(dummy_source_img);
    } else {
        printf("  Failed to create dummy source image for map_photons test.\n");
    }

    // Cleanup files from Test 4
    if (radial_files) {
        free_string_array(radial_files, num_radial_files); // num_radial_files is the count of actual strings
    }


    printf("\n--- Test Suite Finished ---\n");
    return 0;
}
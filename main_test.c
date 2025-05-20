//main_test.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> // For NAN if not in GSL headers directly for C89/C90, and for fabs
#include "schwarzschild_tracer.h" // Your header file

// --- Helper function to create a dummy PPM image ---
PPMImage* create_dummy_ppm(int width, int height) {
    if (width < 0 || height < 0) {
        fprintf(stderr, "create_dummy_ppm: width and height cannot be negative.\n");
        return NULL;
    }
    PPMImage* img = malloc(sizeof(PPMImage));
    if (!img) {
        perror("create_dummy_ppm: malloc for PPMImage struct failed");
        return NULL;
    }
    img->width = width;
    img->height = height;
    img->channels = 3;
    size_t data_size = (size_t)width * height * img->channels;

    if (data_size == 0) {
        img->data = NULL;
    } else {
        img->data = malloc(data_size);
        if (!img->data) {
            perror("create_dummy_ppm: malloc for image data failed");
            free(img);
            return NULL;
        }
        for (int r_idx = 0; r_idx < height; ++r_idx) { // Renamed loop var to avoid conflict
            for (int c_idx = 0; c_idx < width; ++c_idx) { // Renamed loop var
                size_t idx = ((size_t)r_idx * width + c_idx) * img->channels;
                img->data[idx + 0] = (unsigned char)((width > 1) ? (c_idx * 255) / (width - 1) : (c_idx * 255));
                img->data[idx + 1] = (unsigned char)((height > 1) ? (r_idx * 255) / (height - 1) : (r_idx * 255));
                img->data[idx + 2] = 128;
            }
        }
    }
    return img;
}

// --- Unit Test Helper: Compare PhotonMapDataPoint ---
int compare_photon_map_data_points(const PhotonMapDataPoint* p1, const PhotonMapDataPoint* p2, double tol) {
    if (fabs(p1->y_window_cart_x - p2->y_window_cart_x) > tol ||
        fabs(p1->y_window_cart_y - p2->y_window_cart_y) > tol ||
        fabs(p1->y_image_cart_x - p2->y_image_cart_x) > tol ||
        fabs(p1->y_image_cart_y - p2->y_image_cart_y) > tol) {
        return 0; // Not equal
    }
    return 1; // Equal
}

// --- Unit Test for Photon Chunk I/O ---
void test_photon_chunk_io() {
    printf("--- Unit Test: Photon Chunk I/O ---\n");
    const char* test_chunk_filename = "test_photon_chunk.bin";
    size_t num_points_to_save = 3;
    PhotonMapDataPoint original_data[3];

    original_data[0] = (PhotonMapDataPoint){1.0, 2.0, 3.0, 4.0};
    original_data[1] = (PhotonMapDataPoint){-1.5, 0.5, 10.0, -20.0};
    original_data[2] = (PhotonMapDataPoint){0.0, 0.0, 0.0, 0.0};

    printf("  Saving %zu points to %s...\n", num_points_to_save, test_chunk_filename);
    if (save_photon_data_chunk(test_chunk_filename, original_data, num_points_to_save) == 0) {
        printf("    Save successful.\n");
    } else {
        printf("    Save FAILED.\n");
        remove(test_chunk_filename);
        return;
    }

    PhotonMapDataPoint* loaded_data = NULL;
    size_t num_points_loaded = 0;
    printf("  Loading points from %s...\n", test_chunk_filename);
    loaded_data = load_photon_data_chunk(test_chunk_filename, &num_points_loaded);

    if (loaded_data) {
        printf("    Load successful. Loaded %zu points.\n", num_points_loaded);
        if (num_points_loaded == num_points_to_save) {
            int all_match = 1;
            for (size_t i = 0; i < num_points_to_save; ++i) {
                if (!compare_photon_map_data_points(&original_data[i], &loaded_data[i], 1e-9)) {
                    printf("    Data mismatch at index %zu:\n", i);
                    printf("      Original: y_win_x=%.3f, y_win_y=%.3f, y_img_x=%.3f, y_img_y=%.3f\n",
                           original_data[i].y_window_cart_x, original_data[i].y_window_cart_y,
                           original_data[i].y_image_cart_x,  original_data[i].y_image_cart_y);
                    printf("      Loaded:   y_win_x=%.3f, y_win_y=%.3f, y_img_x=%.3f, y_img_y=%.3f\n",
                           loaded_data[i].y_window_cart_x, loaded_data[i].y_window_cart_y,
                           loaded_data[i].y_image_cart_x,  loaded_data[i].y_image_cart_y);
                    all_match = 0;
                }
            }
            if (all_match) {
                printf("    SUCCESS: Original and loaded data match.\n");
            } else {
                printf("    FAILURE: Data mismatch detected.\n");
            }
        } else {
            printf("    FAILURE: Number of points loaded (%zu) does not match saved (%zu).\n",
                   num_points_loaded, num_points_to_save);
        }
        free(loaded_data);
    } else {
        if (num_points_to_save > 0) {
            printf("    Load FAILED or no points loaded (expected %zu).\n", num_points_to_save);
        } else if (num_points_loaded == 0) { // This case is when num_points_to_save was 0
            printf("    Load of 0 points successful (got NULL and 0 points as expected).\n");
        } else { // num_points_to_save was 0, but load returned something unexpected
             printf("    Load of 0 points FAILED (num_points_loaded=%zu, pointer=%p).\n", num_points_loaded, (void*)loaded_data);
        }
    }

    const char* test_zero_chunk_filename = "test_zero_photon_chunk.bin";
    printf("  Saving 0 points to %s...\n", test_zero_chunk_filename);
    if (save_photon_data_chunk(test_zero_chunk_filename, NULL, 0) == 0) {
        printf("    Save of 0 points successful.\n");
        loaded_data = load_photon_data_chunk(test_zero_chunk_filename, &num_points_loaded);
        if (loaded_data == NULL && num_points_loaded == 0) {
            printf("    SUCCESS: Load of 0 points successful (got NULL and 0 points).\n");
        } else {
            printf("    FAILURE: Load of 0 points did not return NULL and 0 points (loaded %zu, ptr %p).\n", num_points_loaded, (void*)loaded_data);
            if(loaded_data) free(loaded_data);
        }
    } else {
        printf("    Save of 0 points FAILED.\n");
    }

    remove(test_chunk_filename);
    remove(test_zero_chunk_filename);
    printf("--- Photon Chunk I/O Test Finished ---\n\n");
}

// --- Unit Test Helper: Compare PPMImages ---
int compare_ppm_images(const PPMImage* img1, const PPMImage* img2, int samples_to_check) {
    if (!img1 && !img2) return 1;
    if (!img1 || !img2) {
        printf("    Image comparison error: one image is NULL, the other is not (img1: %p, img2: %p).\n", (void*)img1, (void*)img2);
        return 0;
    }
    if (img1->width != img2->width || img1->height != img2->height || img1->channels != img2->channels) {
        printf("    Image dimension mismatch:\n");
        printf("      Img1: %dx%d C%d\n", img1->width, img1->height, img1->channels);
        printf("      Img2: %dx%d C%d\n", img2->width, img2->height, img2->channels);
        return 0;
    }

    if (!img1->data && !img2->data) {
        if (img1->width == 0 && img1->height == 0) return 1;
        printf("    Warning: Both images have NULL data pointers but non-zero dimensions (%dx%d).\n", img1->width, img1->height);
        return 1;
    }
    if ((img1->data && !img2->data) || (!img1->data && img2->data)) {
        printf("    Image data pointer mismatch (one NULL, other not, img1->data: %p, img2->data: %p).\n", (void*)img1->data, (void*)img2->data);
        return 0;
    }

    size_t num_pixels = (size_t)img1->width * img1->height;
    if (num_pixels == 0) {
        return 1;
    }

    if (samples_to_check <= 0) samples_to_check = 1;
    if (samples_to_check > (int)num_pixels) samples_to_check = (int)num_pixels;

    for (int i = 0; i < samples_to_check; ++i) {
        size_t px_idx;
        if (samples_to_check == 1 || num_pixels == 1) {
            px_idx = 0;
        } else {
            px_idx = (i * (num_pixels - 1)) / (samples_to_check - 1);
        }
        
        size_t byte_idx = px_idx * img1->channels;

        if (img1->data[byte_idx + 0] != img2->data[byte_idx + 0] ||
            (img1->channels > 1 && img1->data[byte_idx + 1] != img2->data[byte_idx + 1]) ||
            (img1->channels > 2 && img1->data[byte_idx + 2] != img2->data[byte_idx + 2])) {
            printf("    Pixel data mismatch at approx pixel index %zu (byte %zu):\n", px_idx, byte_idx);
            printf("      Img1: (%u, %u, %u)\n", img1->data[byte_idx], (img1->channels > 1 ? img1->data[byte_idx+1] : 0), (img1->channels > 2 ? img1->data[byte_idx+2] : 0));
            printf("      Img2: (%u, %u, %u)\n", img2->data[byte_idx], (img2->channels > 1 ? img2->data[byte_idx+1] : 0), (img2->channels > 2 ? img2->data[byte_idx+2] : 0));
            return 0;
        }
    }
    return 1;
}

// --- Unit Test for PPM Image I/O ---
void test_ppm_image_io() {
    printf("--- Unit Test: PPM Image I/O ---\n");
    const char* test_ppm_filename = "test_image_output.ppm";
    int width = 10, height = 7;

    printf("  Creating dummy PPM image (%dx%d)...\n", width, height);
    PPMImage* original_image = create_dummy_ppm(width, height);
    if (!original_image) {
        printf("    FAILURE: Could not create dummy PPM image for test_ppm_image_io.\n");
        return;
    }
    printf("    Dummy image created.\n");

    printf("  Saving dummy image to %s...\n", test_ppm_filename);
    if (save_ppm_image(test_ppm_filename, original_image) == 0) {
        printf("    Save successful.\n");
    } else {
        printf("    Save FAILED.\n");
        free_ppm_image(original_image);
        remove(test_ppm_filename);
        return;
    }

    PPMImage* loaded_image = NULL;
    printf("  Loading image from %s...\n", test_ppm_filename);
    loaded_image = load_ppm_image(test_ppm_filename);

    if (loaded_image) {
        printf("    Load successful.\n");
        if (compare_ppm_images(original_image, loaded_image, 5)) {
            printf("    SUCCESS: Original and loaded PPM images match.\n");
        } else {
            printf("    FAILURE: PPM image data mismatch detected.\n");
        }
        free_ppm_image(loaded_image);
    } else {
        printf("    Load FAILED.\n");
    }

    free_ppm_image(original_image);
    remove(test_ppm_filename);
    printf("--- PPM Image I/O Test Finished ---\n\n");
}

int main() {
    printf("--- Schwarzschild Tracer Test Suite ---\n\n");
    gsl_set_error_handler_off();

    test_photon_chunk_io();
    test_ppm_image_io();

// --- Test 1: compute_trajectory ---
printf("--- Test 1: compute_trajectory (scattering-like with x_targets) ---\n");
gsl_vector *xt_test1 = gsl_vector_alloc(1);
gsl_vector_set(xt_test1, 0, 0.0); // Example: check crossing of y-axis (x=0)

PhotonTrajectory *traj1 = compute_trajectory(20.0, 0.0, 1.0, -1.25, 100.0,
                                            NAN, false, xt_test1, DEFAULT_T_END, 100); // Increased num_interp_points
if (traj1) {
    if (traj1->error_code == 0 && traj1->K && traj1->K->size > 0) {
        // ... (existing prints for first/last point) ...
        printf("  Checking x_target crossings for compute_trajectory:\n");
        if (traj1->crossings_y_at_x_targets && traj1->num_x_targets > 0) {
            for (size_t i = 0; i < traj1->num_x_targets; ++i) {
                printf("    Target x=%.2f: %zu crossings found.\n",
                       gsl_vector_get(xt_test1, i),
                       traj1->crossings_y_at_x_targets[i] ? traj1->crossings_y_at_x_targets[i]->size : 0);
                if (traj1->crossings_y_at_x_targets[i] && traj1->crossings_y_at_x_targets[i]->size > 0) {
                    printf("      First y-crossing: %.4f\n", gsl_vector_get(traj1->crossings_y_at_x_targets[i], 0));
                }
            }
        } else {
            printf("    No x_target crossings reported or num_x_targets is 0.\n");
        }
    } else {
        // ... (existing error prints) ...
    }
    free_photon_trajectory(traj1);
} else {
    // ... (existing error prints) ...
}
gsl_vector_free(xt_test1);
printf("\n");

    // --- Test 2: compute_trajectory_crossings_only (inward radial) ---
    printf("--- Test 2: compute_trajectory_crossings_only (inward radial) ---\n");
    gsl_vector *x_targets_test2 = gsl_vector_alloc(3);
    gsl_vector_set(x_targets_test2, 0, 4.0);
    gsl_vector_set(x_targets_test2, 1, 3.0);
    gsl_vector_set(x_targets_test2, 2, 2.1);
    TrajectoryCrossings *crossings1 = compute_trajectory_crossings_only(
        5.0, 0.0, 1.0, -M_PI/2, 10.0, 0.0, true, x_targets_test2, DEFAULT_T_END);

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

    // --- Test 3: image_map (using your new promising parameters) ---
    printf("--- Test 3: image_map (using promising parameters) ---\n");
    double test3_x0 = -10.0;
    double test3_x1 = 19.0;
    double test3_x2 = 20.0;
    double test3_Y = 0.5;
    double test3_Z = 0.0;
    double test3_a_scale = 1.0;

    printf("  Parameters: Y=%.2f, Z=%.2f, x0=%.2f, x1=%.2f, x2=%.2f, a_scale=%.2f\n",
           test3_Y, test3_Z, test3_x0, test3_x1, test3_x2, test3_a_scale);

    ImageMapResult im_res1 = image_map(test3_Y, test3_Z, test3_x0, test3_x1, test3_x2, test3_a_scale);
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

    // --- Test 4: results_radial (using your new promising parameters) ---
    printf("--- Test 4: results_radial (using promising parameters) ---\n");
    int num_radial_files = 0;
    char **radial_files = NULL;
    double test4_x0 = -10.0;
    double test4_x1 = 19.0;
    double test4_x2 = 20.0;
    double test4_R_max_sample = 1.0;
    int test4_n_param = 1000;
    size_t test4_chunk_size = 1e6;

    printf("  Parameters: x0=%.2f, x1=%.2f, x2=%.2f, R_max_sample=%.2f, n=%d, chunk_size=%zu\n",
           test4_x0, test4_x1, test4_x2, test4_R_max_sample, test4_n_param, test4_chunk_size);

    radial_files = results_radial(test4_x0, test4_x1, test4_x2,
                                  test4_R_max_sample, test4_n_param,
                                  test4_chunk_size, &num_radial_files);
    if (radial_files && num_radial_files > 0) {
        printf("  results_radial SUCCESS. Created %d file(s):\n", num_radial_files);
        for (int i = 0; i < num_radial_files; ++i) {
            if (radial_files[i]) printf("    File: %s\n", radial_files[i]);
        }
        if (radial_files[0]) {
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
    } else {
        printf("  results_radial FAILED or produced 0 files.\n");
    }
    printf("\n");

    // --- Test 5: map_photons (using loaded rainbow.ppm and new radial_files) ---
    printf("--- Test 5: map_photons (using loaded rainbow.ppm) ---\n");
    PPMImage *source_image_for_test5 = NULL;
    const char* source_image_filename_test5 = "rainbow.ppm";

    printf("  Loading source image '%s' for map_photons test...\n", source_image_filename_test5);
    source_image_for_test5 = load_ppm_image(source_image_filename_test5);

    if (source_image_for_test5) {
        printf("  Successfully loaded source image '%s' (%dx%d).\n",
            source_image_filename_test5, source_image_for_test5->width, source_image_for_test5->height);

        if (radial_files && num_radial_files > 0 && radial_files[0] != NULL) {
            int output_dims[2] = {256, 256};
            RGBColor bg_color = {10, 10, 10};
            const char* lensed_output_filename = "test_lensed_rainbow.ppm";

            printf("  Attempting to render map_photons with %d chunk file(s) from '%s' to '%s'.\n",
                   num_radial_files, source_image_filename_test5, lensed_output_filename);
            int map_status = map_photons(
                NULL, source_image_for_test5,
                radial_files, num_radial_files,
                lensed_output_filename,
                NULL, 0.0,
                output_dims,
                bg_color, false);

            if (map_status == 0) {
                printf("  map_photons SUCCESS. Output image '%s' should be created.\n", lensed_output_filename);
            } else {
                printf("  map_photons FAILED. Status: %d\n", map_status);
            }
        } else {
            printf("  Skipping map_photons test processing as no chunk files were generated by results_radial or list is empty.\n");
        }
        free_ppm_image(source_image_for_test5);
    } else {
        printf("  Failed to load source image '%s' for map_photons test. Ensure it's a valid P6 PPM file and in the correct directory.\n", source_image_filename_test5);
    }

    if (radial_files) {
        free_string_array(radial_files, num_radial_files);
    }

    printf("\n--- Test Suite Finished ---\n");


   // --- Test 6: results_cartesian and map_photons ---
printf("--- Test 6: results_cartesian & map_photons ---\n");
int num_cartesian_files = 0;
char **cartesian_files = NULL;

// Parameters for results_cartesian (adjust as needed)
double test6_x0 = -10.0;
double test6_x1 = 19.0;
double test6_x2 = 20.0;
double test6_a_param = 1.0; // Max Y/Z for sampling
int test6_n_param = 10;     // Grid density
size_t test6_chunk_size = DEFAULT_CHUNK_SIZE_PHOTONS; // From your header

printf("  Calling results_cartesian with: x0=%.2f, x1=%.2f, x2=%.2f, a=%.2f, n=%d, chunk_size=%zu\n",
       test6_x0, test6_x1, test6_x2, test6_a_param, test6_n_param, test6_chunk_size);

cartesian_files = results_cartesian(test6_x0, test6_x1, test6_x2,
                                   test6_a_param, test6_n_param,
                                   test6_chunk_size, &num_cartesian_files);

if (cartesian_files && num_cartesian_files > 0) {
    printf("  results_cartesian SUCCESS. Created %d file(s):\n", num_cartesian_files);
    for (int i = 0; i < num_cartesian_files; ++i) {
        if (cartesian_files[i]) printf("    File: %s\n", cartesian_files[i]);
    }

    // Now use these files with map_photons
    PPMImage *source_image_test6 = load_ppm_image("rainbow.ppm");
    if (source_image_test6) {
        printf("  Successfully loaded source image 'rainbow.ppm' for Cartesian map_photons test.\n");
        int output_dims_test6[2] = {256, 256};
        RGBColor bg_color_test6 = {15, 0, 15}; // Different background for easy identification
        const char* lensed_cartesian_filename = "test_lensed_cartesian.ppm";

        printf("  Attempting map_photons with Cartesian chunks to '%s'.\n", lensed_cartesian_filename);
        int map_status_test6 = map_photons(
            NULL, source_image_test6,
            cartesian_files, num_cartesian_files,
            lensed_cartesian_filename,
            NULL, 0.0, // Using auto-bounds mode
            output_dims_test6,
            bg_color_test6, false);

        if (map_status_test6 == 0) {
            printf("  map_photons SUCCESS with Cartesian data. Output image '%s' should be created.\n", lensed_cartesian_filename);
        } else {
            printf("  map_photons FAILED with Cartesian data. Status: %d\n", map_status_test6);
        }
        free_ppm_image(source_image_test6);
    } else {
        printf("  Failed to load 'rainbow.ppm' for Cartesian map_photons test.\n");
    }
    free_string_array(cartesian_files, num_cartesian_files);
} else {
    printf("  results_cartesian FAILED or produced 0 files.\n");
}
printf("\n");


    return 0;

}
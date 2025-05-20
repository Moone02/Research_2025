// main_test.c
// Top
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> // For NAN, fabs, ceil, round, M_PI
#include <float.h> // For DBL_MAX
#include "schwarzschild_tracer.h" // Your main library header

// --- Helper function to create a dummy PPM image (copied from earlier version for completeness if not in a shared test_util) ---
// If you create a test_utils.c/h, this could live there.
PPMImage* main_test_create_dummy_ppm(int width, int height) { // Renamed to avoid conflict if linked with unit tests
    if (width < 0 || height < 0) {
        fprintf(stderr, "main_test_create_dummy_ppm: width and height cannot be negative.\n");
        return NULL;
    }
    PPMImage* img = malloc(sizeof(PPMImage));
    if (!img) {
        perror("main_test_create_dummy_ppm: malloc for PPMImage struct failed");
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
            perror("main_test_create_dummy_ppm: malloc for image data failed");
            free(img);
            return NULL;
        }
        for (int r_idx = 0; r_idx < height; ++r_idx) {
            for (int c_idx = 0; c_idx < width; ++c_idx) {
                size_t idx = ((size_t)r_idx * width + c_idx) * img->channels;
                img->data[idx + 0] = (unsigned char)((width > 1 && width - 1 != 0) ? (c_idx * 255) / (width - 1) : (c_idx * 255));
                img->data[idx + 1] = (unsigned char)((height > 1 && height -1 != 0) ? (r_idx * 255) / (height - 1) : (r_idx * 255));
                img->data[idx + 2] = 128;
            }
        }
    }
    return img;
}


int main() {
    printf("--- Schwarzschild Tracer Integration Test Suite ---\n\n");
    gsl_set_error_handler_off(); 

    // --- Test I/O Utilities (Optional, if not covered by individual unit tests) ---
    // For brevity, assuming test_photon_chunk_io() and test_ppm_image_io() 
    // are run via 'make run_unit_tests' targeting their respective .c files.
    // If you want them here too, uncomment and ensure their definitions are available or linked.
    // test_photon_chunk_io();
    // test_ppm_image_io();

    char **radial_files = NULL;
    int num_radial_files = 0;
    char **cartesian_files = NULL;
    int num_cartesian_files = 0;

    // Parameters for results_... functions
    double x0_param = -10.0;
    double x1_param = 19.0;
    double x2_param = 20.0;

    // --- Test: results_radial (generates data for map_photons) ---
    printf("--- Integration Test: results_radial ---\n");
    double rr_R_max_sample = 1.0; 
    int    rr_n_param = 100;     // Adjusted for quicker integration testing
    size_t rr_chunk_size = DEFAULT_CHUNK_SIZE_PHOTONS / 10; 

    printf("  Parameters: x0=%.2f, x1=%.2f, x2=%.2f, R_max_sample=%.2f, n=%d, chunk_size=%zu\n",
           x0_param, x1_param, x2_param, rr_R_max_sample, rr_n_param, rr_chunk_size);

    radial_files = results_radial(x0_param, x1_param, x2_param,
                                  rr_R_max_sample, rr_n_param,
                                  rr_chunk_size, &num_radial_files);
    if (radial_files && num_radial_files > 0) {
        printf("  results_radial SUCCESS. Created %d file(s). First file: %s\n", num_radial_files, radial_files[0] ? radial_files[0] : "NULL");
    } else {
        printf("  results_radial FAILED or produced 0 files.\n");
    }
    printf("\n");


    // --- Test: map_photons with data from results_radial ---
    printf("--- Integration Test: map_photons with results_radial data ---\n");
    PPMImage *source_img_r = load_ppm_image("rainbow.ppm"); 
    if (source_img_r) {
        printf("  Loaded source image 'rainbow.ppm' (%dx%d) for map_photons radial test.\n", source_img_r->width, source_img_r->height);
        if (radial_files && num_radial_files > 0) {
            RGBColor bg_color_r = {0, 0, 0}; // Black background
            
            // Test Mode A (Full View of Source Hits)
            double source_log_width_rA = 40.0; 
            int target_output_w_rA = 300; 
            const char* lensed_file_rA = "main_test_lensed_radial_ModeA.ppm"; // Distinct filename
            printf("  map_photons Mode A (Full View): src_width=%.1f, target_out_w=%d to '%s'\n", 
                   source_log_width_rA, target_output_w_rA, lensed_file_rA);
            int status_rA = map_photons(NULL, source_img_r, radial_files, num_radial_files, lensed_file_rA,
                                        source_log_width_rA, target_output_w_rA, 
                                        NULL, // observer_window_bounds is NULL for Mode A
                                        bg_color_r);
            if (status_rA == 0) printf("    SUCCESS: %s created.\n", lensed_file_rA);
            else printf("    FAILURE: map_photons Mode A status %d.\n", status_rA);

            // Test Mode B (Windowed View of Source Hits)
            double source_log_width_rB = 25.0; 
            int target_output_w_rB = 200;
            double obs_window_rB[4] = {-0.5, 0.5, -0.25, 0.25}; // y0min, y0max, z0min, z0max
            const char* lensed_file_rB = "main_test_lensed_radial_ModeB.ppm"; // Distinct filename
            printf("  map_photons Mode B (Windowed): src_width=%.1f, target_out_w=%d, window=[%.2f,%.2f,%.2f,%.2f] to '%s'\n", 
                   source_log_width_rB, target_output_w_rB, obs_window_rB[0], obs_window_rB[1], obs_window_rB[2], obs_window_rB[3], lensed_file_rB);
            int status_rB = map_photons(NULL, source_img_r, radial_files, num_radial_files, lensed_file_rB,
                                        source_log_width_rB, target_output_w_rB, 
                                        obs_window_rB, 
                                        bg_color_r);
            if (status_rB == 0) printf("    SUCCESS: %s created.\n", lensed_file_rB);
            else printf("    FAILURE: map_photons Mode B status %d.\n", status_rB);

        } else { printf("  Skipping map_photons with radial_files as none were generated.\n"); }
        free_ppm_image(source_img_r);
        source_img_r = NULL;
    } else { printf("  Failed to load rainbow.ppm for map_photons radial test.\n");}
    
    if (radial_files) { 
        free_string_array(radial_files, num_radial_files); 
        radial_files = NULL; 
        num_radial_files = 0;
    }
    printf("\n");


    // --- Test: results_cartesian (generates data for map_photons) ---
    printf("--- Integration Test: results_cartesian ---\n");
    double rc_a_param = 1.0; 
    int    rc_n_param = 200;  // Using your increased value
    size_t rc_chunk_size = DEFAULT_CHUNK_SIZE_PHOTONS; 

    printf("  Parameters: x0=%.2f, x1=%.2f, x2=%.2f, a=%.2f, n=%d, chunk_size=%zu\n",
           x0_param, x1_param, x2_param, rc_a_param, rc_n_param, rc_chunk_size);

    cartesian_files = results_cartesian(x0_param, x1_param, x2_param,
                                       rc_a_param, rc_n_param,
                                       rc_chunk_size, &num_cartesian_files);
    if (cartesian_files && num_cartesian_files > 0) {
        printf("  results_cartesian SUCCESS. Created %d file(s). First file: %s\n", num_cartesian_files, cartesian_files[0] ? cartesian_files[0] : "NULL");
    } else {
        printf("  results_cartesian FAILED or produced 0 files.\n");
    }
    printf("\n");

    // --- Test: map_photons with data from results_cartesian ---
    printf("--- Integration Test: map_photons with results_cartesian data ---\n");
    PPMImage *source_img_c = load_ppm_image("rainbow.ppm");
    if (source_img_c) {
        printf("  Loaded source image 'rainbow.ppm' (%dx%d) for map_photons cartesian test.\n", source_img_c->width, source_img_c->height);
        if (cartesian_files && num_cartesian_files > 0) {
            RGBColor bg_color_c = {0, 0, 0}; // Black background
            
            // Test Mode B (Windowed View) with Cartesian data
            double source_log_width_cB = 30.0; 
            int target_output_w_cB = 280;
            double obs_window_cB[4] = {-1.0, 1.0, -1.0, 1.0}; // Window from your last successful run's Mode 1
            const char* lensed_file_cB = "main_test_lensed_cartesian_ModeB.ppm"; // Distinct filename
            printf("  map_photons Mode B (Windowed): src_width=%.1f, target_out_w=%d, window=[%.1f,%.1f,%.1f,%.1f] to '%s'\n", 
                   source_log_width_cB, target_output_w_cB, obs_window_cB[0], obs_window_cB[1], obs_window_cB[2], obs_window_cB[3], lensed_file_cB);
            int status_cB = map_photons(NULL, source_img_c, cartesian_files, num_cartesian_files, lensed_file_cB,
                                        source_log_width_cB, target_output_w_cB, 
                                        obs_window_cB, 
                                        bg_color_c);
            if (status_cB == 0) printf("    SUCCESS: %s created.\n", lensed_file_cB);
            else printf("    FAILURE: map_photons Mode B status %d.\n", status_cB);

        } else { printf("  Skipping map_photons with cartesian_files as none were generated.\n"); }
        free_ppm_image(source_img_c);
        source_img_c = NULL;
    } else { printf("  Failed to load rainbow.ppm for map_photons cartesian test.\n");}
    
    if (cartesian_files) { 
        free_string_array(cartesian_files, num_cartesian_files); 
        cartesian_files = NULL;
        num_cartesian_files = 0;
    }
    printf("\n");

    printf("--- Integration Test Suite Finished ---\n");
    return 0;
}

// End of main_test.c
// Bottom
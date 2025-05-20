// photon_rendering.c
// Top

#define _POSIX_C_SOURCE 200809L 
#include "schwarzschild_tracer.h"
// Includes like stdio.h, stdlib.h, string.h, math.h, float.h, gsl/gsl_sys.h
// should be covered by schwarzschild_tracer.h

// Note: load_ppm_image, save_ppm_image, free_ppm_image,
// load_photon_data_chunk are declared in schwarzschild_tracer.h
// and will be defined in their respective .c files (ppm_image_io.c, photon_chunk_io.c).

int map_photons(
    const char *image_source_path, PPMImage *image_source_data_param,
    char **photon_chunk_files, int num_chunk_files_param,
    const char *save_path,
    double source_logical_width_input,
    int target_output_pixel_width,
    const double *observer_window_bounds, // If NULL, Mode A, else Mode B
    RGBColor default_color_param) {

    int mapped_h = 1, mapped_w = 1;
    double dest_map_bound_min_y = 0.0, dest_map_bound_max_y = 0.0;
    double dest_map_bound_min_z = 0.0, dest_map_bound_max_z = 0.0;

    PPMImage *source_image_loaded = NULL;
    const PPMImage *source_image_to_use = NULL;
    unsigned char *mapped_image_pixels = NULL;
    double *sum_r_array = NULL, *sum_g_array = NULL, *sum_b_array = NULL;
    int64_t *count_array_pixels = NULL;
    int final_status = -1;

    // --- 1. Validate Inputs & Load Source Image ---
    if (target_output_pixel_width <= 0) {
        fprintf(stderr, "map_photons: target_output_pixel_width must be positive.\n");
        goto cleanup_map_photons_early;
    }
    if (source_logical_width_input <= EPSILON_GENERAL) {
        fprintf(stderr, "map_photons: source_logical_width_input must be positive.\n");
        goto cleanup_map_photons_early;
    }

    if (image_source_data_param) {
        source_image_to_use = image_source_data_param;
    } else if (image_source_path) {
        source_image_loaded = load_ppm_image(image_source_path);
        if (!source_image_loaded) {
            fprintf(stderr, "map_photons: Failed to load source image from '%s'.\n", image_source_path);
            goto cleanup_map_photons_early;
        }
        source_image_to_use = source_image_loaded;
    } else {
        fprintf(stderr, "map_photons: No source image provided (path or data).\n"); goto cleanup_map_photons_early;
    }
    if (!source_image_to_use || source_image_to_use->width <= 0 || source_image_to_use->height <= 0 || source_image_to_use->channels != 3) {
        fprintf(stderr, "map_photons: Invalid source image properties.\n"); goto cleanup_map_photons_early;
    }
    int orig_h = source_image_to_use->height; int orig_w = source_image_to_use->width;
    double source_pixel_aspect = (orig_h > 0) ? ((double)orig_w / orig_h) : 1.0;
    printf("map_photons: Source image '%s' loaded (%dx%d, aspect: %.3f).\n",
           image_source_path ? image_source_path : "from data", orig_w, orig_h, source_pixel_aspect);

    // --- Determine num_chunks_to_process EARLY ---
    int num_chunks_to_process = 0;
    if (photon_chunk_files) {
        if (num_chunk_files_param >= 0) {
            num_chunks_to_process = num_chunk_files_param;
        } else {
            for (num_chunks_to_process = 0; photon_chunk_files[num_chunks_to_process] != NULL; ++num_chunks_to_process);
        }
    }
    printf("map_photons: Number of photon chunk files to process: %d\n", num_chunks_to_process);

    // --- 2. Define Source Logical Area to Sample From (common to both modes) ---
    double source_logical_height_calculated;
    if (source_pixel_aspect > EPSILON_GENERAL) {
        source_logical_height_calculated = source_logical_width_input / source_pixel_aspect;
    } else {
        fprintf(stderr, "map_photons: Warning - source image aspect ratio invalid. Setting source logical height equal to width.\n");
        source_logical_height_calculated = source_logical_width_input;
    }
    source_logical_height_calculated = fmax(source_logical_height_calculated, EPSILON_GENERAL);

    double source_map_bound_y1_extent = source_logical_width_input / 2.0; 
    double source_map_bound_z1_extent = source_logical_height_calculated / 2.0;
    double source_denom_y1 = fmax(source_logical_width_input, EPSILON_GENERAL);    
    double source_denom_z1 = fmax(source_logical_height_calculated, EPSILON_GENERAL); 
    printf("Source logical sampling area: y1 in [%.3f, %.3f], z1 in [%.3f, %.3f]\n",
           -source_map_bound_y1_extent, source_map_bound_y1_extent,
           -source_map_bound_z1_extent, source_map_bound_z1_extent);

    // --- 3. Determine Mode and Output Pixel/Logical Dimensions ---
    mapped_w = target_output_pixel_width; 

    if (observer_window_bounds == NULL) { // MODE A: "Full View of Source Hits"
        printf("Operating in Mode A: Full View of Source Hits.\n");
        printf("  Z-axis rendering: Logical MIN_z maps to image top, logical MAX_z to image bottom.\n");

        if (source_pixel_aspect > EPSILON_GENERAL) {
            mapped_h = (int)round((double)mapped_w / source_pixel_aspect);
        } else {
            mapped_h = mapped_w; 
        }
        if (mapped_h <= 0) mapped_h = 1;
        printf("  Output pixel dimensions (WxH): %d x %d (matches source aspect).\n", mapped_w, mapped_h);

        double min_y0_scan = DBL_MAX, max_y0_scan = -DBL_MAX;
        double min_z0_scan = DBL_MAX, max_z0_scan = -DBL_MAX;
        bool found_photons_for_mode_a = false;

        if (num_chunks_to_process > 0) {
            printf("  Scanning photon chunks for Mode A destination bounds...\n");
            for (int ci_scan = 0; ci_scan < num_chunks_to_process; ++ci_scan) {
                size_t num_pts_scan;
                PhotonMapDataPoint* chunk_scan = load_photon_data_chunk(photon_chunk_files[ci_scan], &num_pts_scan);
                if (chunk_scan && num_pts_scan > 0) {
                    for (size_t k = 0; k < num_pts_scan; ++k) {
                        if (fabs(chunk_scan[k].y_image_cart_x) > source_map_bound_y1_extent + EPSILON_GENERAL ||
                            fabs(chunk_scan[k].y_image_cart_y) > source_map_bound_z1_extent + EPSILON_GENERAL) {
                            continue; 
                        }
                        if (gsl_isnan(chunk_scan[k].y_window_cart_x) || gsl_isnan(chunk_scan[k].y_window_cart_y)) continue;

                        found_photons_for_mode_a = true;
                        if (chunk_scan[k].y_window_cart_x < min_y0_scan) min_y0_scan = chunk_scan[k].y_window_cart_x;
                        if (chunk_scan[k].y_window_cart_x > max_y0_scan) max_y0_scan = chunk_scan[k].y_window_cart_x;
                        if (chunk_scan[k].y_window_cart_y < min_z0_scan) min_z0_scan = chunk_scan[k].y_window_cart_y;
                        if (chunk_scan[k].y_window_cart_y > max_z0_scan) max_z0_scan = chunk_scan[k].y_window_cart_y;
                    }
                }
                free(chunk_scan);
            }
        }

        if (!found_photons_for_mode_a) {
            printf("  Warning: No photons found hitting the source area for Mode A. Output will be background color.\n");
            dest_map_bound_min_y = -1.0; dest_map_bound_max_y = 1.0;
            dest_map_bound_min_z = -1.0; dest_map_bound_max_z = 1.0;
        } else {
            double scanned_logical_width = fmax(max_y0_scan - min_y0_scan, EPSILON_GENERAL);
            double scanned_logical_height = fmax(max_z0_scan - min_z0_scan, EPSILON_GENERAL);
            double scanned_aspect = (scanned_logical_height > EPSILON_GENERAL) ? (scanned_logical_width / scanned_logical_height) : 1.0;
            printf("  Mode A: Scanned photon window extents: y0=[%.3f, %.3f], z0=[%.3f, %.3f] (Aspect: %.3f)\n",
                   min_y0_scan, max_y0_scan, min_z0_scan, max_z0_scan, scanned_aspect);

            if (scanned_aspect > source_pixel_aspect) { 
                double scanned_z_center = (min_z0_scan + max_z0_scan) / 2.0;
                double new_logical_height = scanned_logical_width / source_pixel_aspect;
                dest_map_bound_min_y = min_y0_scan;
                dest_map_bound_max_y = max_y0_scan;
                dest_map_bound_min_z = scanned_z_center - new_logical_height / 2.0;
                dest_map_bound_max_z = scanned_z_center + new_logical_height / 2.0;
            } else { 
                double scanned_y_center = (min_y0_scan + max_y0_scan) / 2.0;
                double new_logical_width = scanned_logical_height * source_pixel_aspect;
                dest_map_bound_min_z = min_z0_scan;
                dest_map_bound_max_z = max_z0_scan;
                dest_map_bound_min_y = scanned_y_center - new_logical_width / 2.0;
                dest_map_bound_max_y = scanned_y_center + new_logical_width / 2.0;
            }
        }
        printf("  Mode A: Final logical window to render: y0=[%.3f, %.3f], z0=[%.3f, %.3f] (Aspect: %.3f)\n",
               dest_map_bound_min_y, dest_map_bound_max_y, dest_map_bound_min_z, dest_map_bound_max_z,
               ( (dest_map_bound_max_z - dest_map_bound_min_z) > EPSILON_GENERAL ? 
                 (dest_map_bound_max_y - dest_map_bound_min_y) / (dest_map_bound_max_z - dest_map_bound_min_z) : 0.0 ) );

    } else { // MODE B: "Windowed View of Source Hits"
        printf("Operating in Mode B: Windowed View of Source Hits.\n");
        printf("  Z-axis rendering: Logical MIN_z maps to image top, logical MAX_z to image bottom.\n");

        dest_map_bound_min_y = observer_window_bounds[0];
        dest_map_bound_max_y = observer_window_bounds[1];
        dest_map_bound_min_z = observer_window_bounds[2];
        dest_map_bound_max_z = observer_window_bounds[3];

        if (dest_map_bound_min_y >= dest_map_bound_max_y || dest_map_bound_min_z >= dest_map_bound_max_z) {
            fprintf(stderr, "map_photons: Invalid observer_window_bounds (min >= max).\n");
            goto cleanup_map_photons_early; 
        }

        double window_logical_width = dest_map_bound_max_y - dest_map_bound_min_y;
        double window_logical_height = dest_map_bound_max_z - dest_map_bound_min_z;
        double observer_window_aspect = (window_logical_height > EPSILON_GENERAL) ? (window_logical_width / window_logical_height) : 1.0;

        if (observer_window_aspect > EPSILON_GENERAL) {
            mapped_h = (int)round((double)mapped_w / observer_window_aspect);
        } else {
            mapped_h = mapped_w; 
        }
        if (mapped_h <= 0) mapped_h = 1;
        printf("  Output pixel dimensions (WxH): %d x %d (matches observer window aspect %.3f).\n", mapped_w, mapped_h, observer_window_aspect);
        printf("  User-defined observer window (logical): y0 in [%.3f, %.3f], z0 in [%.3f, %.3f]\n",
               dest_map_bound_min_y, dest_map_bound_max_y, dest_map_bound_min_z, dest_map_bound_max_z);
    }

    // --- 4. Allocate Buffers & Handle Empty Chunk List ---
    size_t num_output_pixels = (size_t)mapped_h * mapped_w;
    mapped_image_pixels = malloc(num_output_pixels * 3 * sizeof(unsigned char));
    sum_r_array = calloc(num_output_pixels, sizeof(double));
    sum_g_array = calloc(num_output_pixels, sizeof(double));
    sum_b_array = calloc(num_output_pixels, sizeof(double));
    count_array_pixels = calloc(num_output_pixels, sizeof(int64_t));
    if (!mapped_image_pixels || !sum_r_array || !sum_g_array || !sum_b_array || !count_array_pixels) {
        fprintf(stderr, "map_photons: Failed to allocate image buffers.\n"); goto cleanup_map_photons;
    }
    for (size_t i = 0; i < num_output_pixels; ++i) {
        mapped_image_pixels[i*3 + 0] = default_color_param.r;
        mapped_image_pixels[i*3 + 1] = default_color_param.g;
        mapped_image_pixels[i*3 + 2] = default_color_param.b;
    }

    if (num_chunks_to_process == 0 && save_path && strlen(save_path) > 0) {
        printf("map_photons: No photon chunk files provided. Saving default background image.\n");
        PPMImage temp_out_img = {mapped_image_pixels, mapped_w, mapped_h, 3};
        if (save_ppm_image(save_path, &temp_out_img) == 0) final_status = 0;
        goto cleanup_map_photons; 
    }
    
    double dest_denom_y0 = fmax(dest_map_bound_max_y - dest_map_bound_min_y, EPSILON_GENERAL);
    double dest_denom_z0 = fmax(dest_map_bound_max_z - dest_map_bound_min_z, EPSILON_GENERAL);

    // --- 5. Loop Through Chunks for Actual Pixel Processing & Accumulation ---
    printf("Processing photon data for rendering...\n");
    for (int ci = 0; ci < num_chunks_to_process; ++ci) {
        if ((ci + 1) % 10 == 0 || ci == 0 || ci == num_chunks_to_process - 1 || num_chunks_to_process < 10) {
             printf("\r  Processing chunk %d/%d: %s...", ci + 1, num_chunks_to_process, photon_chunk_files[ci]); fflush(stdout);
        }
        size_t num_pts_in_chunk;
        PhotonMapDataPoint* chunk_data = load_photon_data_chunk(photon_chunk_files[ci], &num_pts_in_chunk);
        if (!chunk_data || num_pts_in_chunk == 0) { free(chunk_data); continue; }

        for (size_t k_photon = 0; k_photon < num_pts_in_chunk; ++k_photon) {
            double y0_f = chunk_data[k_photon].y_window_cart_x; double z0_f = chunk_data[k_photon].y_window_cart_y;
            double y1_f = chunk_data[k_photon].y_image_cart_x;  double z1_f = chunk_data[k_photon].y_image_cart_y;

            if (gsl_isnan(y0_f) || gsl_isnan(z0_f) || gsl_isnan(y1_f) || gsl_isnan(z1_f)) { continue; }

            if (fabs(y1_f) > source_map_bound_y1_extent + EPSILON_GENERAL ||
                fabs(z1_f) > source_map_bound_z1_extent + EPSILON_GENERAL) {
                continue;
            }
            if (y0_f < dest_map_bound_min_y - EPSILON_GENERAL || y0_f > dest_map_bound_max_y + EPSILON_GENERAL ||
                z0_f < dest_map_bound_min_z - EPSILON_GENERAL || z0_f > dest_map_bound_max_z + EPSILON_GENERAL) {
                continue;
            }

            int src_col_idx = (int)round(((y1_f + source_map_bound_y1_extent) / source_denom_y1) * (orig_w - 1));
            int src_row_idx = (int)round(((z1_f + source_map_bound_z1_extent) / source_denom_z1) * (orig_h - 1));
            src_col_idx = (src_col_idx < 0) ? 0 : (src_col_idx >= orig_w ? orig_w - 1 : src_col_idx);
            src_row_idx = (src_row_idx < 0) ? 0 : (src_row_idx >= orig_h ? orig_h - 1 : src_row_idx);

            int dest_col_idx = (int)round(((y0_f - dest_map_bound_min_y) / dest_denom_y0) * (mapped_w - 1));
            int dest_row_idx = (int)round(((z0_f - dest_map_bound_min_z) / dest_denom_z0) * (mapped_h - 1)); 

            if (dest_col_idx < 0 || dest_col_idx >= mapped_w || dest_row_idx < 0 || dest_row_idx >= mapped_h) { continue; }

            size_t dest_pixel_offset = (size_t)dest_row_idx * mapped_w + dest_col_idx;
            size_t src_pixel_offset  = (size_t)src_row_idx * orig_w * 3 + src_col_idx * 3;
            
            sum_r_array[dest_pixel_offset] += source_image_to_use->data[src_pixel_offset + 0];
            sum_g_array[dest_pixel_offset] += source_image_to_use->data[src_pixel_offset + 1];
            sum_b_array[dest_pixel_offset] += source_image_to_use->data[src_pixel_offset + 2];
            count_array_pixels[dest_pixel_offset]++;
        }
        free(chunk_data);
    }
    printf("\nPhoton processing complete.\n");

    // --- 6. Finalize Image by Averaging ---
    for (size_t i_pix = 0; i_pix < num_output_pixels; ++i_pix) {
        if (count_array_pixels[i_pix] > 0) {
            mapped_image_pixels[i_pix*3 + 0] = (unsigned char)fmax(0, fmin(255, round(sum_r_array[i_pix] / count_array_pixels[i_pix])));
            mapped_image_pixels[i_pix*3 + 1] = (unsigned char)fmax(0, fmin(255, round(sum_g_array[i_pix] / count_array_pixels[i_pix])));
            mapped_image_pixels[i_pix*3 + 2] = (unsigned char)fmax(0, fmin(255, round(sum_b_array[i_pix] / count_array_pixels[i_pix])));
        }
    }

    // --- 7. Save and Cleanup ---
    PPMImage final_img = {mapped_image_pixels, mapped_w, mapped_h, 3};
    if (save_path && strlen(save_path) > 0) { 
        if (save_ppm_image(save_path, &final_img) == 0) {
            final_status = 0;
            printf("Successfully saved lensed image to '%s'.\n", save_path);
        } else {
            fprintf(stderr, "map_photons: Failed to save final image to '%s'.\n", save_path);
            // final_status remains -1
        }
    } else {
        printf("map_photons: Warning - no save_path provided or path is empty. Image not saved.\n");
        final_status = 0; // Processing was successful, just no save attempted
    }

cleanup_map_photons: 
    if(mapped_image_pixels) { free(mapped_image_pixels); mapped_image_pixels = NULL;}
    if(sum_r_array) { free(sum_r_array); sum_r_array = NULL;}
    if(sum_g_array) { free(sum_g_array); sum_g_array = NULL;}
    if(sum_b_array) { free(sum_b_array); sum_b_array = NULL;}
    if(count_array_pixels) { free(count_array_pixels); count_array_pixels = NULL;}
cleanup_map_photons_early: 
    if(source_image_loaded) { free_ppm_image(source_image_loaded); source_image_loaded = NULL;}
    return final_status;
}


#ifdef UNIT_TEST_PHOTON_RENDERING
// For a proper unit test of map_photons, you'd need:
// 1. Mocked/dummy photon chunk files with known data.
// 2. A known small source image.
// 3. Pre-calculated expected output image for a given set of parameters.
// This is more of an integration test.
// The main_test.c provides better integration testing.
// Here, we can do a very basic smoke test if chunk files are created by other tests.

int main() {
    printf("--- Unit Test for photon_rendering (map_photons smoke test) ---\n");
    gsl_set_error_handler_off();

    // Create a dummy source image for the test
    PPMImage* dummy_source = create_dummy_ppm(64, 64); // Using helper from main_test.c (needs to be available or reimplemented)
                                                      // For a true unit test, this helper might be part of a test utils file.
                                                      // For now, this will only work if linked with main_test.o or if create_dummy_ppm is also in this file.
                                                      // Let's assume for now this file would be compiled standalone with its own helpers if needed.
                                                      // Or, better yet, use a tiny actual PPM file for testing.
    if (!dummy_source) {
        printf("Failed to create dummy source for map_photons test.\n");
        return 1;
    }
    // Create a dummy photon chunk file
    const char* dummy_chunk_file = "test_map_photons_chunk.bin";
    PhotonMapDataPoint pdata[2];
    pdata[0] = (PhotonMapDataPoint){0.1, 0.1, 0.5, 0.5}; // y0x, y0y, y1x, y1y
    pdata[1] = (PhotonMapDataPoint){-0.1, -0.1, -0.5, -0.5};
    if (save_photon_data_chunk(dummy_chunk_file, pdata, 2) != 0) {
        printf("Failed to save dummy chunk for map_photons test.\n");
        free_ppm_image(dummy_source);
        return 1;
    }
    char* chunk_list[] = {(char*)dummy_chunk_file, NULL};

    RGBColor bg = {10,20,30};
    const char* out_file_A = "test_map_photons_ModeA.ppm";
    const char* out_file_B = "test_map_photons_ModeB.ppm";
    double obs_window_B[4] = {-0.2, 0.2, -0.2, 0.2};

    printf("Testing map_photons Mode A...\n");
    int statusA = map_photons(NULL, dummy_source, chunk_list, 1, out_file_A,
                              1.0, // source_logical_width_input
                              100, // target_output_pixel_width
                              NULL, // observer_window_bounds (Mode A)
                              bg);
    if (statusA == 0) printf("Mode A smoke test PASSED (file %s potentially created).\n", out_file_A);
    else printf("Mode A smoke test FAILED (status %d).\n", statusA);

    printf("Testing map_photons Mode B...\n");
    int statusB = map_photons(NULL, dummy_source, chunk_list, 1, out_file_B,
                              1.0, // source_logical_width_input
                              100, // target_output_pixel_width
                              obs_window_B, // observer_window_bounds (Mode B)
                              bg);
    if (statusB == 0) printf("Mode B smoke test PASSED (file %s potentially created).\n", out_file_B);
    else printf("Mode B smoke test FAILED (status %d).\n", statusB);

    free_ppm_image(dummy_source);
    remove(dummy_chunk_file);
    remove(out_file_A);
    remove(out_file_B);

    printf("--- Unit Test for photon_rendering Finished ---\n");
    return 0;
}
// Need helper function create_dummy_ppm if not linked with main_test.o
// For simplicity, assume it would be defined here or in a common test util for standalone unit test.
#ifndef UNIT_TEST_MAIN_FILE // Avoid redefinition if main_test.c includes this as a source
PPMImage* create_dummy_ppm(int width, int height) {
    // ... (Simplified dummy implementation for standalone test if needed) ...
    // ... (This is not ideal, better to link with a common test utils or use actual small files) ...
    // For this exercise, assume this block would be filled or linked if UNIT_TEST_PHOTON_RENDERING is defined
    // and this file is compiled as a standalone test. The version in main_test.c is more complete.
    PPMImage* img = malloc(sizeof(PPMImage));
    if (!img) return NULL;
    img->width = width; img->height = height; img->channels = 3;
    if (width > 0 && height > 0) {
        img->data = calloc((size_t)width * height * 3, 1);
        if (!img->data) { free(img); return NULL;}
    } else {
        img->data = NULL;
    }
    return img;
}
#endif // UNIT_TEST_MAIN_FILE

#endif // UNIT_TEST_PHOTON_RENDERING


// photon_rendering.c
// Bottom
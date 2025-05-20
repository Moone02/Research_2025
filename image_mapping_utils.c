// image_mapping_utils.c
// Top

#include "schwarzschild_tracer.h" // For ImageMapResult, TrajectoryCrossings, gsl_vector, NAN, EPSILON_GENERAL, DEFAULT_T_END, M_PI
// #include <math.h> // For sqrt, fabs, acos, fmax, fmin - already in schwarzschild_tracer.h
// #include <gsl/gsl_vector.h> // Already in schwarzschild_tracer.h

// --- Image Mapping Functions ---

ImageMapResult image_map(
    double Y, double Z, double x_0_plane, double x_1_window, double x_2_observer, double a_scale) {
    ImageMapResult im_result = { .miss = 1, .y_window = NAN, .y_image = NAN, .error_code = 0 };

    // Ensure x_2_observer > x_1_window > x_0_plane for typical reverse ray tracing setup
    // and M=1 context as per Python reference for these mapping functions.
    // The 'a_scale' in Python's Image_map is related to the window size for r_max calculation.
    // Python: a_prime = (x_2 - x_0) / (x_2 - x_1) * a if (x_2 - x_1) != 0 else a
    // Python: r_max = np.sqrt(2 * a_prime**2 + x_0**2)
    // Python also adds a base amount to r_max: if r_max < 2*M + 100: r_max = 2*M + 100
    // For M=1, this is r_max < 102, r_max = 102.

    double M_mapping = 1.0; // Fixed Mass for these specific mapping functions as per Python.

    double observer_window_diff = x_2_observer - x_1_window;
    double a_prime;
    if (fabs(observer_window_diff) < EPSILON_GENERAL) {
        // Avoid division by zero if observer is at the window plane
        a_prime = a_scale;
    } else {
        a_prime = ((x_2_observer - x_0_plane) / observer_window_diff) * a_scale;
    }

    double r_max_integration = sqrt(2.0 * a_prime * a_prime + x_0_plane * x_0_plane);
    // Ensure r_max is sufficiently large, especially if geometry is compact
    if (r_max_integration < (2.0 * M_mapping + 100.0)) {
        r_max_integration = (2.0 * M_mapping + 100.0);
    }
    // Also ensure r_max is greater than observer position x_2_observer if x_2_observer is large
    if (r_max_integration < x_2_observer + 10.0) { // Add some buffer
         r_max_integration = x_2_observer + 10.0;
    }


    double rho_sq = Y * Y + Z * Z;
    double observer_to_window_dist_sq_on_plane = observer_window_diff * observer_window_diff + rho_sq;

    if (observer_to_window_dist_sq_on_plane <= EPSILON_GENERAL) { // Avoid sqrt of zero/negative or division by zero
        im_result.error_code = 1; // Geometrically problematic
        return im_result;
    }

    double arg_acos = sqrt(rho_sq) / sqrt(observer_to_window_dist_sq_on_plane);
    // Clamp arg_acos to [-1, 1] to prevent domain errors with acos due to precision
    arg_acos = fmax(-1.0, fmin(1.0, arg_acos));
    double psi_calc = -acos(arg_acos); // Matches Python's psi calculation for Image_map

    gsl_vector *x_targets_map = gsl_vector_alloc(2);
    if (!x_targets_map) {
        im_result.error_code = -1; // Allocation error
        return im_result;
    }
    gsl_vector_set(x_targets_map, 0, x_0_plane);  // Image/Source plane
    gsl_vector_set(x_targets_map, 1, x_1_window); // Window plane

    // x_stop value should prevent trajectory from going too far "behind" the source plane
    double x_stop_map = x_0_plane - 1.0; // Or some other reasonable offset

    // Call compute_trajectory_crossings_only
    // r_0 is observer's position, phi_0 is 0 by convention for these 2D plane calculations
    TrajectoryCrossings *crossings = compute_trajectory_crossings_only(
        x_2_observer, 0.0, M_mapping, psi_calc, r_max_integration,
        x_stop_map, true, // x_stop active
        x_targets_map, DEFAULT_T_END);

    gsl_vector_free(x_targets_map); // Free after use

    if (!crossings) {
        im_result.error_code = -1; // Crossing computation failed to allocate struct
        return im_result;
    }
    if (crossings->error_code != 0) {
        im_result.error_code = crossings->error_code; // Propagate error from trajectory computation
        free_trajectory_crossings(crossings);
        return im_result;
    }

    // Check if we have crossings for both targets
    // y_list[0] in Python corresponds to x_targets_map[0] (x_0_plane)
    // y_list[1] in Python corresponds to x_targets_map[1] (x_1_window)
    if (crossings->num_x_targets == 2 &&
        crossings->crossings_y_at_x_targets &&
        crossings->crossings_y_at_x_targets[0] && crossings->crossings_y_at_x_targets[0]->size > 0 &&
        crossings->crossings_y_at_x_targets[1] && crossings->crossings_y_at_x_targets[1]->size > 0) {
        
        im_result.y_image = gsl_vector_get(crossings->crossings_y_at_x_targets[0], 0); // First crossing of source plane
        im_result.y_window = gsl_vector_get(crossings->crossings_y_at_x_targets[1], 0); // First crossing of window plane
        im_result.miss = 0; // It's a hit
    } else {
        im_result.miss = 1; // Miss if any target wasn't crossed or data structure is unexpected
    }

    free_trajectory_crossings(crossings);
    return im_result;
}

ImageMapResult image_map_radial(
    double r_rho, double x_0_plane, double x_1_window, double x_2_observer, double a_scale_factor) {
    // This function simply calls the Cartesian version with Z=0.0 and Y=r_rho
    return image_map(r_rho, 0.0, x_0_plane, x_1_window, x_2_observer, a_scale_factor);
}


#ifdef UNIT_TEST_IMAGE_MAPPING_UTILS
// Requires trajectory_module.c to be linked for compute_trajectory_crossings_only
// and schwarzschild_tracer.h for definitions.

int main() {
    printf("--- Unit Test for image_mapping_utils ---\n");
    gsl_set_error_handler_off();

    // Test Case 1: Radial shot (Y=0, Z=0 or r_rho=0)
    // Expect y_window and y_image to be ~0.0
    double x0_test1 = 3.0;
    double x1_test1 = 10.0;
    double x2_test1 = 20.0;
    double a_scale_test1 = 1.0;
    printf("Test 1: image_map_radial (r_rho=0)...\n");
    ImageMapResult res_radial = image_map_radial(0.0, x0_test1, x1_test1, x2_test1, a_scale_test1);
    if (res_radial.error_code == 0 && !res_radial.miss) {
        printf("  HIT: y_window=%.6e, y_image=%.6e\n", res_radial.y_window, res_radial.y_image);
        if (fabs(res_radial.y_window) < 1e-6 && fabs(res_radial.y_image) < 1e-6) {
            printf("  Test 1 PASSED.\n");
        } else {
            printf("  Test 1 FAILED: Expected y_window and y_image near 0 for radial shot.\n");
        }
    } else {
        printf("  Test 1 FAILED: Miss or error_code %d.\n", res_radial.error_code);
    }

    // Test Case 2: Cartesian shot from your main_test.c
    double x0_test2 = -10.0;
    double x1_test2 = 19.0;
    double x2_test2 = 20.0;
    double Y_test2 = 0.5;
    double Z_test2 = 0.0; // Keep it in the plane for simplicity like radial
    double a_scale_test2 = 1.0;
    printf("Test 2: image_map (Y=0.5, Z=0)...\n");
    ImageMapResult res_cart = image_map(Y_test2, Z_test2, x0_test2, x1_test2, x2_test2, a_scale_test2);
     if (res_cart.error_code == 0 && !res_cart.miss) {
        printf("  HIT: y_window=%.4f, y_image=%.4f\n", res_cart.y_window, res_cart.y_image);
        // Expected values from your previous successful main_test.c output:
        // y_window=0.5264, y_image=7.8063
        if (fabs(res_cart.y_window - 0.5264) < 1e-3 && fabs(res_cart.y_image - 7.8063) < 1e-3) {
             printf("  Test 2 PASSED (matches expected values).\n");
        } else {
             printf("  Test 2 FAILED: Values deviate from expected (y_w_exp=0.5264, y_i_exp=7.8063).\n");
        }
    } else {
        printf("  Test 2 FAILED: Miss or error_code %d.\n", res_cart.error_code);
    }

    // Test Case 3: Observer inside event horizon (should ideally be caught by compute_trajectory_crossings_only or its core)
    // compute_trajectory_crossings_only itself checks if r_0 (which is x_2_observer here) <= 2M
    // M_mapping is 1.0 in image_map.
    printf("Test 3: image_map (observer at x=1.0, M=1.0)...\n");
    ImageMapResult res_inside = image_map(0.0, 0.0, -10.0, 5.0, 1.0, 1.0);
    if (res_inside.miss == 1 && res_inside.error_code == 0) { // Expecting a miss, but no allocation error from image_map itself
        // The underlying compute_trajectory_crossings_only should return empty crossings or an error.
        // If it returns empty crossings, image_map sets miss=1.
        printf("  Test 3 PASSED (correctly reported as MISS or underlying error caught by trajectory func).\n");
    } else {
        printf("  Test 3 FAILED: Expected miss or error for observer inside horizon (miss=%d, error_code=%d).\n", res_inside.miss, res_inside.error_code);
    }


    printf("--- Unit Test for image_mapping_utils Finished ---\n");
    return 0;
}
#endif // UNIT_TEST_IMAGE_MAPPING_UTILS

// image_mapping_utils.c
// Bottom
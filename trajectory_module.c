// trajectory_module.c
// Top

#define _POSIX_C_SOURCE 200809L 
#include "schwarzschild_tracer.h" 
#include <float.h> // For DBL_MAX, DBL_EPSILON (ensure this is here)
// Other includes like math.h, string.h, stdio.h, stdlib.h, gsl/* should be in schwarzschild_tracer.h

// --- Helper Macros ---
#define CHECK_GSL_ALLOC_VEC_CORE(ptr, func_name, item_name_str, err_var, err_val, cleanup_label) \
    if ((ptr) == NULL) { fprintf(stderr, "%s: GSL vector allocation failed for %s\n", func_name, item_name_str); err_var = err_val; goto cleanup_label; }
#define CHECK_ALLOC_GEN_CORE(ptr, func_name, type_name_str, err_var, err_val, cleanup_label) \
    if ((ptr) == NULL) { fprintf(stderr, "%s: Standard allocation failed for %s\n", func_name, type_name_str); err_var = err_val; goto cleanup_label; }

// --- Internal Struct Definitions for this Compilation Unit ---
typedef struct {
    EventFunctionParams *event_params_root; 
    double K_prev_step, K_curr_step;
    double y_state_prev_step[2];
    double y_state_curr_step[2];
    double (*get_event_val_func_root)(double, const double[], EventFunctionParams*);
} GslRootEventParams;

typedef struct {
    gsl_vector *K_pts;
    gsl_vector *r_pts;
    gsl_vector *phi_pts;
} TrajectorySegmentDataInternal;


// --- Static Function Definitions (internal helpers) ---

static int schwarzschild_geodesic_eqs(double K, const double y_state[], double f_derivs[], void *params_ptr) {
    (void)K; 
    ODEParams *params = (ODEParams *)params_ptr;
    double M = params->M; double b = params->b; int sign_val = *(params->sign_dr_dk);
    double r = y_state[0];

    if (r <= (2.0 * M + EVENT_DETECTION_TOLERANCE * 0.01)) { 
        f_derivs[0] = 0.0; f_derivs[1] = 0.0; return GSL_SUCCESS;
    }
    double metric_term = (1.0 - 2.0 * M / r);
    if (metric_term < 0 && r > 2.0 * M) { metric_term = 0.0; }

    double term_b_r_sq = (b / r) * (b / r) ;
    double fr_val = 1.0 - metric_term * term_b_r_sq;
    if (fr_val < 0) { fr_val = 0.0; }

    f_derivs[0] = sign_val * sqrt(fr_val);
    if (r < EPSILON_GENERAL) { f_derivs[1] = 0.0; }
    else { f_derivs[1] = b / (r * r); }
    return GSL_SUCCESS;
}

static double get_event_val_fr_zero(double K, const double y_state[], EventFunctionParams *params) {
    (void)K; double r = y_state[0]; double M = params->ode_p_event->M; double b = params->ode_p_event->b;
    if (r <= (2.0 * M + EVENT_DETECTION_TOLERANCE * 0.1)) { return 0.0; }
    double metric_term = (1.0 - 2.0 * M / r); double term_b_r_sq = (b / r) * (b / r);
    return 1.0 - metric_term * term_b_r_sq;
}
static double get_event_val_r_leq_2M(double K, const double y_state[], EventFunctionParams *params) {
    (void)K; return y_state[0] - (2.0 * params->M + EVENT_DETECTION_TOLERANCE * 0.1); }
static double get_event_val_r_max(double K, const double y_state[], EventFunctionParams *params) {
    (void)K; return y_state[0] - params->r_max_event; }
static double get_event_val_x_stop(double K, const double y_state[], EventFunctionParams *params) {
    (void)K; double r = y_state[0]; double phi = y_state[1];
    if (r < EPSILON_GENERAL) { return 1.0; }
    return r * cos(phi) - params->x_stop_event; }
static double get_event_val_x_target(double K, const double y_state[], EventFunctionParams *params) {
    (void)K; double r = y_state[0]; double phi = y_state[1];
    if (r < EPSILON_GENERAL) { return 1.0; }
    return r * cos(phi) - params->x_target_event; }

static double gsl_root_event_wrapper(double K_test, void *params_ptr) {
    GslRootEventParams *root_p = (GslRootEventParams *)params_ptr;
    double y_interp[2];
    if (fabs(root_p->K_curr_step - root_p->K_prev_step) < DBL_EPSILON) {
        y_interp[0] = root_p->y_state_curr_step[0];
        y_interp[1] = root_p->y_state_curr_step[1];
    } else {
        double fraction = (K_test - root_p->K_prev_step) / (root_p->K_curr_step - root_p->K_prev_step);
        fraction = fmax(0.0, fmin(1.0, fraction));
        y_interp[0] = root_p->y_state_prev_step[0] + fraction * (root_p->y_state_curr_step[0] - root_p->y_state_prev_step[0]);
        y_interp[1] = root_p->y_state_prev_step[1] + fraction * (root_p->y_state_curr_step[1] - root_p->y_state_prev_step[1]);
    }
    return root_p->get_event_val_func_root(K_test, y_interp, root_p->event_params_root);
}

static int find_event_time_gsl(
    double K_start_interval, double y_start_interval[],
    double K_end_interval, double y_end_interval[],
    double (*get_event_val_func_gsl)(double, const double[], EventFunctionParams*),
    EventFunctionParams *event_p_gsl, int crossing_direction,
    double *K_event_found_gsl, double y_event_found_gsl[]) {

    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T_solver = gsl_root_fsolver_brent;
    gsl_root_fsolver *s_solver = gsl_root_fsolver_alloc(T_solver);
    if (!s_solver) { /*fprintf(stderr, "find_event_time_gsl: Failed to allocate GSL root solver.\n");*/ return 0; }

    gsl_function F_gsl_root;
    GslRootEventParams root_params_instance;
    root_params_instance.get_event_val_func_root = get_event_val_func_gsl;
    root_params_instance.event_params_root = event_p_gsl;
    root_params_instance.K_prev_step = K_start_interval;
    root_params_instance.K_curr_step = K_end_interval;
    memcpy(root_params_instance.y_state_prev_step, y_start_interval, 2 * sizeof(double));
    memcpy(root_params_instance.y_state_curr_step, y_end_interval, 2 * sizeof(double));

    F_gsl_root.function = &gsl_root_event_wrapper;
    F_gsl_root.params = &root_params_instance;

    double val_start = get_event_val_func_gsl(K_start_interval, y_start_interval, event_p_gsl);
    double val_end = get_event_val_func_gsl(K_end_interval, y_end_interval, event_p_gsl);

    bool bracketed = false;
    if (crossing_direction == 0) {
        if (val_start * val_end <= EVENT_DETECTION_TOLERANCE * EVENT_DETECTION_TOLERANCE) { bracketed = true; }
    } else if (crossing_direction < 0) {
        if (val_start > -EVENT_DETECTION_TOLERANCE && val_end <= EVENT_DETECTION_TOLERANCE) { bracketed = true; }
    } else {
        if (val_start < EVENT_DETECTION_TOLERANCE && val_end >= -EVENT_DETECTION_TOLERANCE) { bracketed = true; }
    }
    if (fabs(val_start) < EVENT_DETECTION_TOLERANCE) {
        *K_event_found_gsl = K_start_interval; memcpy(y_event_found_gsl, y_start_interval, 2 * sizeof(double));
        gsl_root_fsolver_free(s_solver); return 1;
    }
    if (fabs(val_end) < EVENT_DETECTION_TOLERANCE && K_end_interval > K_start_interval) {
        *K_event_found_gsl = K_end_interval; memcpy(y_event_found_gsl, y_end_interval, 2 * sizeof(double));
        gsl_root_fsolver_free(s_solver); return 1;
    }
    if (!bracketed) { gsl_root_fsolver_free(s_solver); return 0; }

    int status_set = gsl_root_fsolver_set(s_solver, &F_gsl_root, K_start_interval, K_end_interval);
    if (status_set != GSL_SUCCESS) { gsl_root_fsolver_free(s_solver); return 0; }

    int status_iterate;
    do {
        iter++;
        status_iterate = gsl_root_fsolver_iterate(s_solver);
        if (status_iterate != GSL_SUCCESS && status_iterate != GSL_CONTINUE) {
             gsl_root_fsolver_free(s_solver); return 0;
        }
        *K_event_found_gsl = gsl_root_fsolver_root(s_solver);
        double K_low = gsl_root_fsolver_x_lower(s_solver);
        double K_high = gsl_root_fsolver_x_upper(s_solver);
        status_iterate = gsl_root_test_interval(K_low, K_high, 0, EVENT_DETECTION_TOLERANCE);
    } while (status_iterate == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s_solver);
    if (status_iterate != GSL_SUCCESS) { return 0; }

    if (fabs(K_end_interval - K_start_interval) < DBL_EPSILON) {
         memcpy(y_event_found_gsl, y_end_interval, 2 * sizeof(double));
    } else {
        double fraction = (*K_event_found_gsl - K_start_interval) / (K_end_interval - K_start_interval);
        fraction = fmax(0.0, fmin(1.0, fraction));
        y_event_found_gsl[0] = y_start_interval[0] + fraction * (y_end_interval[0] - y_start_interval[0]);
        y_event_found_gsl[1] = y_start_interval[1] + fraction * (y_end_interval[1] - y_start_interval[1]);
    }
    return 1;
}

static int reallocate_gsl_vector_if_needed(
    gsl_vector **vec_ptr,                   
    size_t current_logical_element_count,   
    size_t *current_physical_capacity_ptr,  
    size_t initial_capacity_val) {          

    if (!vec_ptr || !current_physical_capacity_ptr) {
        return -1; 
    }
    if (current_logical_element_count >= *current_physical_capacity_ptr) {
        size_t new_physical_capacity;
        if (*current_physical_capacity_ptr > 0) {
            new_physical_capacity = *current_physical_capacity_ptr * 2; 
        } else {
            new_physical_capacity = initial_capacity_val; 
        }
        if (new_physical_capacity <= current_logical_element_count) {
            new_physical_capacity = current_logical_element_count + initial_capacity_val;
            if (new_physical_capacity <= current_logical_element_count) { 
                new_physical_capacity = current_logical_element_count + 1; 
            }
        }
        gsl_vector *new_v = gsl_vector_alloc(new_physical_capacity);
        if (!new_v) {
            return -1; 
        }
        gsl_vector *old_v = *vec_ptr; 
        if (old_v) { 
            size_t num_elements_to_copy = current_logical_element_count;
            if (num_elements_to_copy > old_v->size) {
                num_elements_to_copy = old_v->size; 
            }
            if (num_elements_to_copy > 0) {
                 gsl_vector_const_view old_subvector_to_copy = gsl_vector_const_subvector(old_v, 0, num_elements_to_copy);
                 gsl_vector_view new_subvector_to_paste_into = gsl_vector_subvector(new_v, 0, num_elements_to_copy);
                 gsl_vector_memcpy(&new_subvector_to_paste_into.vector, &old_subvector_to_copy.vector);
            }
            gsl_vector_free(old_v); 
        }
        *vec_ptr = new_v; 
        *current_physical_capacity_ptr = new_physical_capacity; 
    }
    return 0; 
}

int gsl_vector_dynamic_append(gsl_vector **vec_ptr, double value_param, size_t *current_capacity_ptr) {
    if (!vec_ptr || !current_capacity_ptr) { return -1; }

    if (!(*vec_ptr)) { 
        *current_capacity_ptr = INITIAL_CROSSING_POINTS_CAPACITY; 
        *vec_ptr = gsl_vector_alloc(*current_capacity_ptr); 
        if (!(*vec_ptr)) { 
            fprintf(stderr, "gsl_vector_dynamic_append: Initial allocation failed.\n");
            return -1; 
        }
        (*vec_ptr)->size = 0; 
    }
    
    // The third argument to reallocate_gsl_vector_if_needed IS current_capacity_ptr itself
    // because reallocate_gsl_vector_if_needed expects a (size_t *)
    if (reallocate_gsl_vector_if_needed(vec_ptr, (*vec_ptr)->size, current_capacity_ptr, INITIAL_CROSSING_POINTS_CAPACITY) != 0) {
        return -1; 
    }
    
    size_t index_to_set = (*vec_ptr)->size;

    // CORRECTED LINE: Use the function's parameter name 'current_capacity_ptr'
    if ((*vec_ptr)->data && index_to_set < *current_capacity_ptr) { 
        (*vec_ptr)->data[index_to_set] = value_param; 
    } else {
        // Also use 'current_capacity_ptr' in the error message
        fprintf(stderr, "DYNAMIC_APPEND_ERROR: Invalid state for direct write. Index: %zu, Tracked Capacity: %zu, Data Ptr: %p\n",
                index_to_set, *current_capacity_ptr, (void*)((*vec_ptr) ? (*vec_ptr)->data : NULL) ); 
        return -2; 
    }

    (*vec_ptr)->size++; 
    return 0;
}

double normalize_phi(double phi) {
    phi = fmod(phi, 2.0 * M_PI);
    if (phi < 0.0) {
        phi += 2.0 * M_PI;
    }
    return phi;
}

static void unwrap_phi_values(const gsl_vector* K_raw, const gsl_vector* phi_raw, gsl_vector* phi_unwrapped) {
    if (!K_raw || !phi_raw || !phi_unwrapped || K_raw->size != phi_raw->size || K_raw->size != phi_unwrapped->size || K_raw->size == 0) {
        return;
    }
    gsl_vector_set(phi_unwrapped, 0, gsl_vector_get(phi_raw, 0));
    double prev_unwrapped_val = gsl_vector_get(phi_unwrapped, 0);

    for (size_t i = 1; i < K_raw->size; ++i) {
        double current_raw_val = gsl_vector_get(phi_raw, i);
        double diff_from_prev_unwrapped = current_raw_val - prev_unwrapped_val;
        
        double adjusted_delta = fmod(diff_from_prev_unwrapped + M_PI, 2.0 * M_PI) - M_PI;
        if (adjusted_delta == -M_PI && diff_from_prev_unwrapped > 0) { 
             adjusted_delta = M_PI;
        }

        double new_unwrapped_val = prev_unwrapped_val + adjusted_delta;
        gsl_vector_set(phi_unwrapped, i, new_unwrapped_val);
        prev_unwrapped_val = new_unwrapped_val;
    }
}

// Static variable to count calls to the core function
static int integrate_photon_trajectory_core_call_count = 0; // File-scope static

// Core integration logic - kept static as it's an internal implementation detail
static int integrate_photon_trajectory_core(
    double r_0, double phi_0, double M_val, double psi, double r_max_val,
    double x_stop_val, bool x_stop_active_flag,
    const gsl_vector *x_targets_vec,
    double t_end_max_overall, 
    double rtol, double atol,
    PhotonTrajectory *full_traj_output,   
    TrajectoryCrossings *crossings_output, 
    int num_interp_points_for_full_traj
) {
    // *** PASTED THE FULL BODY OF integrate_photon_trajectory_core HERE ***
    // *** FROM schwarzschild_tracer_c.txt (THE VERSION WE DEBUGGED EXTENSIVELY) ***
    integrate_photon_trajectory_core_call_count++;
    int current_call_instance = integrate_photon_trajectory_core_call_count;

    int core_error_code = 0;
    double K_final_reached_integration = 0.0;

    gsl_odeiv2_driver *driver = NULL;
    TrajectorySegmentDataInternal *segments_collected_list = NULL;
    size_t num_segments_collected = 0;
    size_t segments_collected_capacity = 0;
    gsl_vector *current_segment_K_temp = NULL, *current_segment_r_temp = NULL, *current_segment_phi_temp = NULL;
    size_t current_segment_point_count = 0; 
    size_t current_segment_K_capacity = 0;    
    size_t current_segment_r_capacity = 0;    
    size_t current_segment_phi_capacity = 0;  
    gsl_vector **crossings_collector_ptr = NULL;
    size_t *crossings_collector_capacities = NULL;
    size_t num_x_targets_val = 0;
    gsl_vector *K_all_raw = NULL, *r_all_raw = NULL, *phi_all_raw_temp = NULL, *phi_unwrapped_all = NULL;
    gsl_spline *spline_r_interp = NULL, *spline_phi_interp = NULL;
    gsl_interp_accel *acc_r_interp = NULL, *acc_phi_interp = NULL;

    if (r_0 <= (2.0 * M_val + EPSILON_GENERAL) || r_0 >= r_max_val) {
        K_final_reached_integration = 0.0;
        if (full_traj_output) {
            int n_pts_trivial = (num_interp_points_for_full_traj > 0) ? num_interp_points_for_full_traj : 1;
            if (n_pts_trivial <=0) n_pts_trivial = 1; // Ensure at least one point
            full_traj_output->K = gsl_vector_calloc(n_pts_trivial); CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->K, "core (trivial K)", "K", core_error_code, -1, cleanup_core);
            full_traj_output->r = gsl_vector_calloc(n_pts_trivial); CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->r, "core (trivial r)", "r", core_error_code, -1, cleanup_core);
            full_traj_output->phi = gsl_vector_calloc(n_pts_trivial); CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->phi, "core (trivial phi)", "phi", core_error_code, -1, cleanup_core);
            full_traj_output->x = gsl_vector_calloc(n_pts_trivial); CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->x, "core (trivial x)", "x", core_error_code, -1, cleanup_core);
            full_traj_output->y = gsl_vector_calloc(n_pts_trivial); CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->y, "core (trivial y)", "y", core_error_code, -1, cleanup_core);

            double phi_norm_init = normalize_phi(phi_0);
            for (int i = 0; i < n_pts_trivial; ++i) {
                gsl_vector_set(full_traj_output->K, i, 0.0);
                gsl_vector_set(full_traj_output->r, i, r_0);
                gsl_vector_set(full_traj_output->phi, i, phi_norm_init);
                gsl_vector_set(full_traj_output->x, i, r_0 * cos(phi_norm_init));
                gsl_vector_set(full_traj_output->y, i, r_0 * sin(phi_norm_init));
            }
            full_traj_output->num_x_targets = x_targets_vec ? x_targets_vec->size : 0;
            if (full_traj_output->num_x_targets > 0) {
                full_traj_output->crossings_y_at_x_targets = calloc(full_traj_output->num_x_targets, sizeof(gsl_vector*));
                CHECK_ALLOC_GEN_CORE(full_traj_output->crossings_y_at_x_targets, "core (trivial cr_y)", "gsl_vec**", core_error_code, -1, cleanup_core);
                for (size_t i = 0; i < full_traj_output->num_x_targets; ++i) {
                    full_traj_output->crossings_y_at_x_targets[i] = gsl_vector_alloc(0); 
                    CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->crossings_y_at_x_targets[i], "core (trivial cr[i])", "cr_y[i]", core_error_code, -1, cleanup_core);
                }
            } else { full_traj_output->crossings_y_at_x_targets = NULL; }
             full_traj_output->error_code = 0; 
        }
        if (crossings_output) {
            crossings_output->num_x_targets = x_targets_vec ? x_targets_vec->size : 0;
            if (crossings_output->num_x_targets > 0) {
                crossings_output->crossings_y_at_x_targets = calloc(crossings_output->num_x_targets, sizeof(gsl_vector*));
                CHECK_ALLOC_GEN_CORE(crossings_output->crossings_y_at_x_targets, "core (trivial cr_out_arr)", "gsl_vec**", core_error_code, -1, cleanup_core);
                for (size_t i = 0; i < crossings_output->num_x_targets; ++i) {
                    crossings_output->crossings_y_at_x_targets[i] = gsl_vector_alloc(0);
                    CHECK_GSL_ALLOC_VEC_CORE(crossings_output->crossings_y_at_x_targets[i], "core (trivial cr_out[i])", "cr_out[i]", core_error_code, -1, cleanup_core);
                }
            } else { crossings_output->crossings_y_at_x_targets = NULL; }
            crossings_output->error_code = 0; 
        }
        goto cleanup_core; 
    }

    double metric_r0_term = (1.0 - 2.0 * M_val / r_0);
    if (metric_r0_term <= EPSILON_GENERAL) {K_final_reached_integration = 0.0; goto cleanup_core;} // Should be caught by r_0 <= 2M

    double b_val = cos(psi) * r_0 / sqrt(metric_r0_term);
    int current_sign_dr_dk;
    // ... (Logic for current_sign_dr_dk as in your file) ...
    if (fabs(b_val) < 1e-9) { 
        double sin_psi_val = sin(psi);
        if (sin_psi_val > EPSILON_GENERAL) { current_sign_dr_dk = 1; }
        else if (sin_psi_val < -EPSILON_GENERAL) { current_sign_dr_dk = -1; }
        else { current_sign_dr_dk = 1; } // Default for ambiguous purely radial
    } else { 
        if (sin(psi) > EPSILON_GENERAL) { current_sign_dr_dk = 1; }
        else if (sin(psi) < -EPSILON_GENERAL) { current_sign_dr_dk = -1; }
        else { current_sign_dr_dk = 1; } // Default for purely tangential
    }

    ODEParams ode_params_instance;
    ode_params_instance.M = M_val;
    ode_params_instance.b = b_val;
    ode_params_instance.sign_dr_dk = &current_sign_dr_dk;

    gsl_odeiv2_system sys = {schwarzschild_geodesic_eqs, NULL, 2, &ode_params_instance};
    double initial_driver_step = (rtol < 1e-9 || atol < 1e-9) ? 1e-5 : 1e-6; // Adjusted initial step
    driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, initial_driver_step, rtol, atol);
    CHECK_ALLOC_GEN_CORE(driver, "core (driver)", "gsl_odeiv2_driver", core_error_code, -2, cleanup_core);

    double K_loop_variable = 0.0; 
    double y_current_state[2] = {r_0, phi_0};

    if (full_traj_output) {
        segments_collected_capacity = INITIAL_SEGMENTS_CAPACITY;
        segments_collected_list = malloc(segments_collected_capacity * sizeof(TrajectorySegmentDataInternal));
        CHECK_ALLOC_GEN_CORE(segments_collected_list, "core (segments_list)", "SegList", core_error_code, -1, cleanup_core);
        for(size_t i=0; i<segments_collected_capacity; ++i) { // Initialize pointers
            segments_collected_list[i].K_pts = NULL; segments_collected_list[i].r_pts = NULL; segments_collected_list[i].phi_pts = NULL;
        }
        current_segment_K_capacity = INITIAL_RAW_POINTS_CAPACITY;
        current_segment_r_capacity = INITIAL_RAW_POINTS_CAPACITY;
        current_segment_phi_capacity = INITIAL_RAW_POINTS_CAPACITY;
        current_segment_K_temp = gsl_vector_alloc(current_segment_K_capacity); 
        CHECK_GSL_ALLOC_VEC_CORE(current_segment_K_temp, "core (K_temp init)", "K_temp", core_error_code, -1, cleanup_core);
        current_segment_r_temp = gsl_vector_alloc(current_segment_r_capacity); 
        CHECK_GSL_ALLOC_VEC_CORE(current_segment_r_temp, "core (r_temp init)", "r_temp", core_error_code, -1, cleanup_core);
        current_segment_phi_temp = gsl_vector_alloc(current_segment_phi_capacity); 
        CHECK_GSL_ALLOC_VEC_CORE(current_segment_phi_temp, "core (phi_temp init)", "phi_temp", core_error_code, -1, cleanup_core);
        
        gsl_vector_set(current_segment_K_temp, 0, K_loop_variable); 
        gsl_vector_set(current_segment_r_temp, 0, y_current_state[0]);
        gsl_vector_set(current_segment_phi_temp, 0, y_current_state[1]);
        current_segment_point_count = 1;
    }

    if (x_targets_vec && (full_traj_output || crossings_output)) {
        num_x_targets_val = x_targets_vec->size;
        if (num_x_targets_val > 0) {
            crossings_collector_ptr = calloc(num_x_targets_val, sizeof(gsl_vector*)); 
            CHECK_ALLOC_GEN_CORE(crossings_collector_ptr, "core (cr_coll_ptr)","gsl_vec**", core_error_code, -1, cleanup_core);
            crossings_collector_capacities = calloc(num_x_targets_val, sizeof(size_t)); 
            CHECK_ALLOC_GEN_CORE(crossings_collector_capacities, "core (cr_coll_cap)","size_t*", core_error_code, -1, cleanup_core);
            for (size_t i = 0; i < num_x_targets_val; ++i) {
                crossings_collector_capacities[i] = INITIAL_CROSSING_POINTS_CAPACITY; // Init capacity for each
                crossings_collector_ptr[i] = gsl_vector_alloc(crossings_collector_capacities[i]); 
                CHECK_GSL_ALLOC_VEC_CORE(crossings_collector_ptr[i], "core (cr_coll[i])","coll[i]", core_error_code, -1, cleanup_core);
                crossings_collector_ptr[i]->size = 0; // Important: GSL vector size tracks used elements
            }
        }
    }
    
    EventFunctionParams event_params_instances[MAX_EVENT_TYPES_CORE]; // MAX_EVENT_TYPES_CORE from .h
    double (*event_get_val_func_ptrs[MAX_EVENT_TYPES_CORE])(double, const double[], EventFunctionParams*);
    bool event_is_terminal[MAX_EVENT_TYPES_CORE];
    int event_crossing_dir[MAX_EVENT_TYPES_CORE];
    int event_is_x_target_idx_map[MAX_EVENT_TYPES_CORE]; 
    size_t num_active_event_funcs = 0;

    // ... (Event setup logic - copied verbatim from your file) ...
    event_params_instances[num_active_event_funcs].ode_p_event = &ode_params_instance;
    event_get_val_func_ptrs[num_active_event_funcs] = get_event_val_fr_zero;
    event_is_terminal[num_active_event_funcs] = true; event_crossing_dir[num_active_event_funcs] = 0;
    event_is_x_target_idx_map[num_active_event_funcs] = -1; num_active_event_funcs++;

    event_params_instances[num_active_event_funcs].M = M_val; 
    event_get_val_func_ptrs[num_active_event_funcs] = get_event_val_r_leq_2M;
    event_is_terminal[num_active_event_funcs] = true; event_crossing_dir[num_active_event_funcs] = -1;
    event_is_x_target_idx_map[num_active_event_funcs] = -1; num_active_event_funcs++;

    event_params_instances[num_active_event_funcs].r_max_event = r_max_val;
    event_get_val_func_ptrs[num_active_event_funcs] = get_event_val_r_max;
    event_is_terminal[num_active_event_funcs] = true; event_crossing_dir[num_active_event_funcs] = 1; 
    event_is_x_target_idx_map[num_active_event_funcs] = -1; num_active_event_funcs++;
    
    if (x_stop_active_flag) {
        event_params_instances[num_active_event_funcs].x_stop_event = x_stop_val;
        event_get_val_func_ptrs[num_active_event_funcs] = get_event_val_x_stop;
        event_is_terminal[num_active_event_funcs] = true; event_crossing_dir[num_active_event_funcs] = -1; 
        event_is_x_target_idx_map[num_active_event_funcs] = -1; num_active_event_funcs++;
    }

    if (x_targets_vec && num_x_targets_val > 0) {
        for (size_t i = 0; i < num_x_targets_val; ++i) {
            if (num_active_event_funcs >= MAX_EVENT_TYPES_CORE) { 
                fprintf(stderr, "integrate_photon_trajectory_core: Too many event types defined.\n"); 
                core_error_code = -2; goto cleanup_core; // Use main cleanup
            }
            event_params_instances[num_active_event_funcs].x_target_event = gsl_vector_get(x_targets_vec, i);
            event_get_val_func_ptrs[num_active_event_funcs] = get_event_val_x_target;
            event_is_terminal[num_active_event_funcs] = false; 
            event_crossing_dir[num_active_event_funcs] = -1;   
            event_is_x_target_idx_map[num_active_event_funcs] = (int)i; 
            num_active_event_funcs++;
        }
    }

    int integration_stop_code = 0; 
    int safety_break_counter = 0;
    const int MAX_SAFETY_BREAK_CORE = 200000; 
    const double MIN_K_PROGRESS_PER_ITER = 1e-12; 

        while (K_loop_variable < t_end_max_overall && integration_stop_code == 0 && safety_break_counter < MAX_SAFETY_BREAK_CORE) {
        safety_break_counter++;
        double K_iter_start_val = K_loop_variable; 
        double y_iter_start_state[2]; 
        memcpy(y_iter_start_state, y_current_state, 2 * sizeof(double));

        double K_target_for_this_gsl_step = K_iter_start_val + MAX_K_PER_GSL_SUPER_STEP;
        if (K_target_for_this_gsl_step > t_end_max_overall) {
            K_target_for_this_gsl_step = t_end_max_overall;
        }
        
        if (K_target_for_this_gsl_step <= K_iter_start_val + MIN_K_PROGRESS_PER_ITER && K_iter_start_val >= t_end_max_overall - EPSILON_GENERAL) {
            integration_stop_code = 1; 
            K_final_reached_integration = t_end_max_overall;
            K_loop_variable = t_end_max_overall;
            break;
        }
        if (K_target_for_this_gsl_step <= K_iter_start_val + MIN_K_PROGRESS_PER_ITER) {
            // Try a much smaller fraction of the super step if the full super step is too small
            K_target_for_this_gsl_step = K_iter_start_val + MAX_K_PER_GSL_SUPER_STEP * 0.01; 
             if (K_target_for_this_gsl_step > t_end_max_overall) K_target_for_this_gsl_step = t_end_max_overall;
             if (K_target_for_this_gsl_step <= K_iter_start_val + MIN_K_PROGRESS_PER_ITER){
                  integration_stop_code = 1;
                  K_final_reached_integration = t_end_max_overall;
                  K_loop_variable = t_end_max_overall;
                  break;
             }
        }
        
        
        double K_gsl_inout = K_iter_start_val; 
        int status = gsl_odeiv2_driver_apply(driver, &K_gsl_inout, K_target_for_this_gsl_step, y_current_state);
        double K_gsl_step_end = K_gsl_inout; 


        if (status != GSL_SUCCESS) {
            fprintf(stderr, "GSL ODE solver error: %s (K_start_iter=%.2e, r=%.2e, phi=%.2e)\n", gsl_strerror(status), K_iter_start_val, y_iter_start_state[0], y_iter_start_state[1]);
            integration_stop_code = 3; K_final_reached_integration = K_iter_start_val; break;
        }

        if (fabs(K_gsl_step_end - K_iter_start_val) < DBL_EPSILON * fabs(K_iter_start_val) + 10.0 * DBL_MIN) {
            if (safety_break_counter > 50) { 
                fprintf(stderr, "Warning: GSL ODE solver stalled at K=%.15e (K_start_iter=%.15e, h=%.3e). Terminating integration.\n", 
                        K_gsl_step_end, K_iter_start_val, driver->h);
                integration_stop_code = 4; 
                K_final_reached_integration = K_gsl_step_end; 
                break;
            }
        }
        K_final_reached_integration = K_gsl_step_end; 

        if (full_traj_output) {
            // Pass the correct individual capacity variable for each vector
            if (reallocate_gsl_vector_if_needed(&current_segment_K_temp, current_segment_point_count, &current_segment_K_capacity, INITIAL_RAW_POINTS_CAPACITY) != 0 ||
                reallocate_gsl_vector_if_needed(&current_segment_r_temp, current_segment_point_count, &current_segment_r_capacity, INITIAL_RAW_POINTS_CAPACITY) != 0 ||
                reallocate_gsl_vector_if_needed(&current_segment_phi_temp, current_segment_point_count, &current_segment_phi_capacity, INITIAL_RAW_POINTS_CAPACITY) != 0) {
                integration_stop_code = 3; K_final_reached_integration = K_iter_start_val; break; 
            }
            gsl_vector_set(current_segment_K_temp, current_segment_point_count, K_gsl_step_end);
            gsl_vector_set(current_segment_r_temp, current_segment_point_count, y_current_state[0]); 
            gsl_vector_set(current_segment_phi_temp, current_segment_point_count, y_current_state[1]);

            
            current_segment_point_count++;
        }

        double earliest_K_event_this_step = K_gsl_step_end + 1.0; 
        int triggered_event_func_original_idx = -1;
        double y_at_earliest_event_this_step[2];
        y_at_earliest_event_this_step[0] = NAN; y_at_earliest_event_this_step[1] = NAN;

        for (size_t ev_idx = 0; ev_idx < num_active_event_funcs; ++ev_idx) {
            double K_event_found_type; double y_event_found_type[2];
            int event_found_flag = find_event_time_gsl(
                K_iter_start_val, y_iter_start_state,    
                K_gsl_step_end, y_current_state,         
                event_get_val_func_ptrs[ev_idx], &event_params_instances[ev_idx],
                event_crossing_dir[ev_idx], &K_event_found_type, y_event_found_type
            );

            if (event_found_flag) {
                bool accept_this_specific_event = false;
                if (event_is_terminal[ev_idx]) {
                    if (K_event_found_type >= K_iter_start_val - EVENT_DETECTION_TOLERANCE &&
                        K_event_found_type <= K_gsl_step_end + EVENT_DETECTION_TOLERANCE) {
                        accept_this_specific_event = true;
                    }
                } else { 
                    if (K_event_found_type > K_iter_start_val + (EVENT_DETECTION_TOLERANCE * 0.1) && 
                        K_event_found_type <= K_gsl_step_end + EVENT_DETECTION_TOLERANCE) {
                        accept_this_specific_event = true;
                    }
                }

                if (accept_this_specific_event) {
                    if (K_event_found_type < earliest_K_event_this_step) {
                        earliest_K_event_this_step = K_event_found_type;
                        triggered_event_func_original_idx = ev_idx;
                        memcpy(y_at_earliest_event_this_step, y_event_found_type, 2 * sizeof(double));
                    } else if (fabs(K_event_found_type - earliest_K_event_this_step) < EVENT_DETECTION_TOLERANCE * 0.01) {
                        if (event_get_val_func_ptrs[ev_idx] == get_event_val_fr_zero) { 
                            earliest_K_event_this_step = K_event_found_type; triggered_event_func_original_idx = ev_idx; memcpy(y_at_earliest_event_this_step, y_event_found_type, 2 * sizeof(double));
                        } else if (event_get_val_func_ptrs[ev_idx] == get_event_val_r_leq_2M && event_get_val_func_ptrs[triggered_event_func_original_idx] != get_event_val_fr_zero) {
                            earliest_K_event_this_step = K_event_found_type; triggered_event_func_original_idx = ev_idx; memcpy(y_at_earliest_event_this_step, y_event_found_type, 2 * sizeof(double));
                        }
                    }
                }
            }
        } 

        if (triggered_event_func_original_idx != -1) { 

            K_loop_variable = earliest_K_event_this_step; 
            memcpy(y_current_state, y_at_earliest_event_this_step, 2 * sizeof(double)); 
            K_final_reached_integration = K_loop_variable; 

            if (full_traj_output && current_segment_point_count > 0) {
                size_t last_idx_in_temp = current_segment_point_count - 1; 
                gsl_vector_set(current_segment_K_temp, last_idx_in_temp, K_loop_variable); 
                gsl_vector_set(current_segment_r_temp, last_idx_in_temp, y_current_state[0]); 
                gsl_vector_set(current_segment_phi_temp, last_idx_in_temp, y_current_state[1]); 
            }

            int x_target_list_idx = event_is_x_target_idx_map[triggered_event_func_original_idx];
            if (x_target_list_idx != -1) { 
                if (crossings_collector_ptr && x_target_list_idx < (int)num_x_targets_val) {
                    double r_event_val = y_at_earliest_event_this_step[0];
                    if (r_event_val >= (2 * M_val + EPSILON_GENERAL) && r_event_val >= 1e-10) { 
                        double y_val_to_add_crossing = y_at_earliest_event_this_step[0] * sin(y_at_earliest_event_this_step[1]);

                        if (gsl_vector_dynamic_append(&crossings_collector_ptr[x_target_list_idx], y_val_to_add_crossing, &crossings_collector_capacities[x_target_list_idx]) != 0) {
                            integration_stop_code = 3; K_final_reached_integration = K_loop_variable; break;
                        }
                    }
                }
                K_loop_variable += EVENT_DETECTION_TOLERANCE * 10.0; 
                K_final_reached_integration = K_loop_variable; 
            } else if (event_is_terminal[triggered_event_func_original_idx]) { 
                if (event_get_val_func_ptrs[triggered_event_func_original_idx] == get_event_val_fr_zero) { 
                    current_sign_dr_dk *= -1;
                    //y_current_state[0] += (EVENT_DETECTION_TOLERANCE * 10.0) * current_sign_dr_dk; 
                    //K_loop_variable += EVENT_DETECTION_TOLERANCE * 10.0;
                    y_current_state[0] += (1e-7) * current_sign_dr_dk; // Match Python's r nudge
                    K_loop_variable += 1e-7;  
                    K_final_reached_integration = K_loop_variable;

                    if (y_current_state[0] <= (2*M_val + EPSILON_GENERAL) || y_current_state[0] >= r_max_val ||
                        (x_stop_active_flag && (y_current_state[0]*cos(y_current_state[1]) <= x_stop_val + EPSILON_GENERAL))) {


                        integration_stop_code = 2; 
                    }

                    if (full_traj_output && integration_stop_code != 2) { 
                        if (num_segments_collected >= segments_collected_capacity) {
                            segments_collected_capacity = (segments_collected_capacity == 0) ? INITIAL_SEGMENTS_CAPACITY : segments_collected_capacity * 2;
                            TrajectorySegmentDataInternal *temp_realloc = realloc(segments_collected_list, segments_collected_capacity * sizeof(TrajectorySegmentDataInternal));
                            CHECK_ALLOC_GEN_CORE(temp_realloc, "core (realloc seg_list)", "SegList", integration_stop_code, 3, cleanup_core); // Use main cleanup on alloc fail here
                            segments_collected_list = temp_realloc;
                        }
                        
                        
                        segments_collected_list[num_segments_collected].K_pts = gsl_vector_alloc(current_segment_point_count); CHECK_GSL_ALLOC_VEC_CORE(segments_collected_list[num_segments_collected].K_pts, "core (seg K)", "K_seg", integration_stop_code, 3, cleanup_core);
                        segments_collected_list[num_segments_collected].r_pts = gsl_vector_alloc(current_segment_point_count); CHECK_GSL_ALLOC_VEC_CORE(segments_collected_list[num_segments_collected].r_pts, "core (seg r)", "r_seg", integration_stop_code, 3, cleanup_core);
                        segments_collected_list[num_segments_collected].phi_pts = gsl_vector_alloc(current_segment_point_count); CHECK_GSL_ALLOC_VEC_CORE(segments_collected_list[num_segments_collected].phi_pts, "core (seg phi)", "phi_seg", integration_stop_code, 3, cleanup_core);
                        
                        for(size_t k=0; k<current_segment_point_count; ++k) {
                            gsl_vector_set(segments_collected_list[num_segments_collected].K_pts, k, gsl_vector_get(current_segment_K_temp, k));
                            gsl_vector_set(segments_collected_list[num_segments_collected].r_pts, k, gsl_vector_get(current_segment_r_temp, k));
                            gsl_vector_set(segments_collected_list[num_segments_collected].phi_pts, k, gsl_vector_get(current_segment_phi_temp, k));
                        }
                        num_segments_collected++;
                        current_segment_point_count = 0; 
                        // Start new segment with the nudged state
                        if (reallocate_gsl_vector_if_needed(&current_segment_K_temp, 0, &current_segment_K_capacity, INITIAL_RAW_POINTS_CAPACITY) != 0 ||
                        reallocate_gsl_vector_if_needed(&current_segment_r_temp, 0, &current_segment_r_capacity, INITIAL_RAW_POINTS_CAPACITY) != 0 ||
                        reallocate_gsl_vector_if_needed(&current_segment_phi_temp, 0, &current_segment_phi_capacity, INITIAL_RAW_POINTS_CAPACITY) != 0) {
                        integration_stop_code = 3; K_final_reached_integration = K_loop_variable; break;
                    }
                    // Since current_segment_point_count is 0, these set the first element (index 0) of the new segment
                    gsl_vector_set(current_segment_K_temp, 0, K_loop_variable);
                    gsl_vector_set(current_segment_r_temp, 0, y_current_state[0]);
                    gsl_vector_set(current_segment_phi_temp, 0, y_current_state[1]);
                    current_segment_point_count = 1; 
                

                        if (current_call_instance == 1 && num_segments_collected == 1) { // We are building the second segment
                        if (current_segment_point_count == 255 || current_segment_point_count == 256 || current_segment_point_count == 257) {
                        printf("DEBUG_SEG1_BUILD (Iter %d, Seg1_LocalIdx %zu): K_gsl_step_end=%.15e, r=%.3e, phi=%.3e. Stored K=%.15e\n",
                       safety_break_counter,
                       current_segment_point_count, // This is the index where it's about to be written
                       K_gsl_step_end,
                       y_current_state[0],
                       y_current_state[1],
                       gsl_vector_get(current_segment_K_temp, current_segment_point_count) // Read back what was just written
                       );
                        fflush(stdout);
                        }
                        

                        }


                        current_segment_point_count = 1;
                     
                    }
                    } else {  integration_stop_code = 2; // This 'else' handles OTHER terminal events (NOT fr_zero)
                    char event_name[50] = "Unknown";
                    if (event_get_val_func_ptrs[triggered_event_func_original_idx] == get_event_val_r_leq_2M) {
                        strcpy(event_name, "r_leq_2M");
                    } else if (event_get_val_func_ptrs[triggered_event_func_original_idx] == get_event_val_r_max) {
                        strcpy(event_name, "r_max");
                    } else if (event_get_val_func_ptrs[triggered_event_func_original_idx] == get_event_val_x_stop) {
                        strcpy(event_name, "x_stop");
                    }
                    // Add more 'else if' for any other distinct terminal event types

                    integration_stop_code = 2; // Terminate integration due to this event
                }


                

                
            }

        // Near the end of the if (triggered_event_func_original_idx != -1) block
        } else { 
            K_loop_variable = K_gsl_step_end; 
        }

        if (integration_stop_code != 0) { break; } 

        if (K_loop_variable >= t_end_max_overall - EPSILON_GENERAL) { 
            integration_stop_code = 1; 
            K_final_reached_integration = t_end_max_overall; 
            K_loop_variable = t_end_max_overall; 
        }
        if (fabs(K_loop_variable - K_iter_start_val) < MIN_K_PROGRESS_PER_ITER && safety_break_counter > 100 && integration_stop_code == 0) {
            fprintf(stderr, "Warning: Minimal K progress in iteration %d (K_start_iter=%.15e, K_loop_var_end_iter=%.15e). Terminating.\n",
                    safety_break_counter, K_iter_start_val, K_loop_variable);
            integration_stop_code = 5; 
            K_final_reached_integration = K_loop_variable;
            break;
        }
    } 

    if (safety_break_counter >= MAX_SAFETY_BREAK_CORE) {
        fprintf(stderr, "Warning: Max safety break counter (%d) reached. K_final=%.15e\n", MAX_SAFETY_BREAK_CORE, K_final_reached_integration);
        if (integration_stop_code == 0) { integration_stop_code = 2; } 
    }
    // ... (other termination messages integration_stop_code 4, 5) ...

    // --- Final Segment Storage (if any points remain in current_segment_...) ---
    if (full_traj_output && current_segment_point_count > 0 && integration_stop_code != 3) { 
        if (num_segments_collected >= segments_collected_capacity) {
            segments_collected_capacity = (segments_collected_capacity == 0) ? INITIAL_SEGMENTS_CAPACITY : segments_collected_capacity * 2;
            TrajectorySegmentDataInternal *temp_realloc = realloc(segments_collected_list, segments_collected_capacity * sizeof(TrajectorySegmentDataInternal));
            CHECK_ALLOC_GEN_CORE(temp_realloc, "core (realloc final seg_list)", "SegList", integration_stop_code, 3, cleanup_core_post_loop); 
            segments_collected_list = temp_realloc;
        }
        // This segment K_pts etc should be allocated with current_segment_point_count
        segments_collected_list[num_segments_collected].K_pts = gsl_vector_alloc(current_segment_point_count); CHECK_GSL_ALLOC_VEC_CORE(segments_collected_list[num_segments_collected].K_pts, "core (final K_seg)", "K_seg", integration_stop_code, 3, cleanup_core_post_loop);
        segments_collected_list[num_segments_collected].r_pts = gsl_vector_alloc(current_segment_point_count); CHECK_GSL_ALLOC_VEC_CORE(segments_collected_list[num_segments_collected].r_pts, "core (final r_seg)", "r_seg", integration_stop_code, 3, cleanup_core_post_loop);
        segments_collected_list[num_segments_collected].phi_pts = gsl_vector_alloc(current_segment_point_count); CHECK_GSL_ALLOC_VEC_CORE(segments_collected_list[num_segments_collected].phi_pts, "core (final phi_seg)", "phi_seg", integration_stop_code, 3, cleanup_core_post_loop);
        
        for(size_t k=0; k<current_segment_point_count; ++k) {
            gsl_vector_set(segments_collected_list[num_segments_collected].K_pts, k, gsl_vector_get(current_segment_K_temp, k));
            gsl_vector_set(segments_collected_list[num_segments_collected].r_pts, k, gsl_vector_get(current_segment_r_temp, k));
            gsl_vector_set(segments_collected_list[num_segments_collected].phi_pts, k, gsl_vector_get(current_segment_phi_temp, k));
        }

        
        num_segments_collected++;
    }
cleanup_core_post_loop: {} 



    // --- Interpolation Logic ---
    // Check if full trajectory output is requested and no critical error occurred during integration
    if (full_traj_output && integration_stop_code != 3 && integration_stop_code != 4 && integration_stop_code != 5) {
        
        // Scenario 1: Trivial case where output might have been pre-filled by the function's initial r0 checks.
        if (num_segments_collected == 0 && 
            fabs(K_final_reached_integration - 0.0) < DBL_EPSILON &&
            (full_traj_output->K && full_traj_output->K->size >= 1 && fabs(gsl_vector_get(full_traj_output->K,0) - 0.0) < DBL_EPSILON)) {
            // Data was likely set by the initial trivial case handler. No further interpolation needed.

        } 
        // Scenario 2: Segments were collected, and integration ended reasonably. Proceed to interpolate.
        else if (num_segments_collected > 0 && K_final_reached_integration >= -EPSILON_GENERAL) {
            size_t total_raw_points = 0;
            for (size_t s = 0; s < num_segments_collected; ++s) {
                if (segments_collected_list[s].K_pts) { 
                     total_raw_points += segments_collected_list[s].K_pts->size;
                } else {
                    fprintf(stderr, "Warning C%d: Null K_pts in segment %zu during interpolation prep.\n", current_call_instance, s);
                }
            }


            if (total_raw_points >= 1) { 
                int n_pts_interp_actual = (num_interp_points_for_full_traj > 0) ? num_interp_points_for_full_traj : 1;
                
                if (full_traj_output->K == NULL || full_traj_output->K->size != (size_t)n_pts_interp_actual) {
                    gsl_vector_free(full_traj_output->K); full_traj_output->K = NULL;
                    gsl_vector_free(full_traj_output->r); full_traj_output->r = NULL;
                    gsl_vector_free(full_traj_output->phi); full_traj_output->phi = NULL;
                    gsl_vector_free(full_traj_output->x); full_traj_output->x = NULL;
                    gsl_vector_free(full_traj_output->y); full_traj_output->y = NULL;
                    
                    full_traj_output->K = gsl_vector_calloc(n_pts_interp_actual); // Use calloc for safety
                    CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->K, "core (interp K)", "K", integration_stop_code, 3, cleanup_core_interp);
                    full_traj_output->r = gsl_vector_calloc(n_pts_interp_actual); 
                    CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->r, "core (interp r)", "r", integration_stop_code, 3, cleanup_core_interp);
                    full_traj_output->phi = gsl_vector_calloc(n_pts_interp_actual); 
                    CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->phi, "core (interp phi)", "phi", integration_stop_code, 3, cleanup_core_interp);
                    full_traj_output->x = gsl_vector_calloc(n_pts_interp_actual); 
                    CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->x, "core (interp x)", "x", integration_stop_code, 3, cleanup_core_interp);
                    full_traj_output->y = gsl_vector_calloc(n_pts_interp_actual); 
                    CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->y, "core (interp y)", "y", integration_stop_code, 3, cleanup_core_interp);
                }

                K_all_raw = gsl_vector_alloc(total_raw_points); 
                CHECK_GSL_ALLOC_VEC_CORE(K_all_raw, "core (K_all_raw)", "K_all_raw", integration_stop_code, 3, cleanup_core_interp_arrays);
                r_all_raw = gsl_vector_alloc(total_raw_points); 
                CHECK_GSL_ALLOC_VEC_CORE(r_all_raw, "core (r_all_raw)", "r_all_raw", integration_stop_code, 3, cleanup_core_interp_arrays);
                phi_all_raw_temp = gsl_vector_alloc(total_raw_points); 
                CHECK_GSL_ALLOC_VEC_CORE(phi_all_raw_temp, "core (phi_all_raw_temp)", "phi_all_raw_temp", integration_stop_code, 3, cleanup_core_interp_arrays);

                size_t current_concat_idx = 0;
                for (size_t s = 0; s < num_segments_collected; ++s) {
                     if (!segments_collected_list[s].K_pts || segments_collected_list[s].K_pts->size == 0) continue; 
                    for (size_t i = 0; i < segments_collected_list[s].K_pts->size; ++i) {
                        if (current_concat_idx < total_raw_points) { 
                            gsl_vector_set(K_all_raw, current_concat_idx, gsl_vector_get(segments_collected_list[s].K_pts, i));
                            gsl_vector_set(r_all_raw, current_concat_idx, gsl_vector_get(segments_collected_list[s].r_pts, i));
                            gsl_vector_set(phi_all_raw_temp, current_concat_idx, gsl_vector_get(segments_collected_list[s].phi_pts, i));
                            current_concat_idx++;
                        } else {
                            fprintf(stderr, "Warning C%d: Concatenation index %zu exceeded total_raw_points %zu.\n", current_call_instance, current_concat_idx, total_raw_points);
                            goto end_concat_loop; 
                        }
                    }
                }
                end_concat_loop:;
                
                if (current_concat_idx < total_raw_points) {
                    fprintf(stderr, "Warning C%d: Final concatenated points (%zu) < initially calculated total_raw_points (%zu). Using actual count.\n", current_call_instance, current_concat_idx, total_raw_points);
                    total_raw_points = current_concat_idx;
                     if (total_raw_points == 0) {
                        fprintf(stderr, "Error C%d: No valid raw points to interpolate after concat adjustments.\n", current_call_instance);
                        if (full_traj_output->K) { gsl_vector_free(full_traj_output->K); full_traj_output->K = NULL;}
                        // (and r,phi,x,y)
                        goto cleanup_core_interp_arrays; 
                     }
                }
                

                phi_unwrapped_all = gsl_vector_alloc(total_raw_points); 
                CHECK_GSL_ALLOC_VEC_CORE(phi_unwrapped_all, "core (phi_unwrapped)", "phi_unwrapped", integration_stop_code, 3, cleanup_core_interp_arrays);
                
                gsl_vector_const_view K_all_raw_v = gsl_vector_const_subvector(K_all_raw, 0, total_raw_points);
                gsl_vector_const_view phi_all_raw_temp_v = gsl_vector_const_subvector(phi_all_raw_temp, 0, total_raw_points);
                unwrap_phi_values(&K_all_raw_v.vector, &phi_all_raw_temp_v.vector, phi_unwrapped_all);
                
                


                if (total_raw_points == 1) { 
                    for (size_t i = 0; i < full_traj_output->K->size; ++i) { 
                        gsl_vector_set(full_traj_output->K, i, gsl_vector_get(K_all_raw, 0));
                        double r_val = gsl_vector_get(r_all_raw, 0);
                        double phi_val_norm = normalize_phi(gsl_vector_get(phi_unwrapped_all, 0));
                        gsl_vector_set(full_traj_output->r, i, r_val);
                        gsl_vector_set(full_traj_output->phi, i, phi_val_norm);
                        gsl_vector_set(full_traj_output->x, i, r_val * cos(phi_val_norm));
                        gsl_vector_set(full_traj_output->y, i, r_val * sin(phi_val_norm));
                    }
                } else { // total_raw_points > 1

                    const gsl_interp_type *spline_type_cspline = gsl_interp_cspline;
                    const gsl_interp_type *spline_type_linear = gsl_interp_linear;
                    const gsl_interp_type *spline_type = (total_raw_points >= gsl_interp_type_min_size(spline_type_cspline)) ? spline_type_cspline : spline_type_linear;

                    spline_r_interp = gsl_spline_alloc(spline_type, total_raw_points); 
                    CHECK_ALLOC_GEN_CORE(spline_r_interp, "core (spline_r)", "gsl_spline", integration_stop_code, 3, cleanup_core_interp_splines);
                    spline_phi_interp = gsl_spline_alloc(spline_type, total_raw_points); 
                    CHECK_ALLOC_GEN_CORE(spline_phi_interp, "core (spline_phi)", "gsl_spline", integration_stop_code, 3, cleanup_core_interp_splines);
                    acc_r_interp = gsl_interp_accel_alloc(); 
                    CHECK_ALLOC_GEN_CORE(acc_r_interp, "core (acc_r)", "gsl_interp_accel", integration_stop_code, 3, cleanup_core_interp_splines);
                    acc_phi_interp = gsl_interp_accel_alloc(); 
                    CHECK_ALLOC_GEN_CORE(acc_phi_interp, "core (acc_phi)", "gsl_interp_accel", integration_stop_code, 3, cleanup_core_interp_splines);

                    gsl_vector_const_view r_all_raw_v = gsl_vector_const_subvector(r_all_raw, 0, total_raw_points);
                    
                    gsl_spline_init(spline_r_interp, K_all_raw_v.vector.data, r_all_raw_v.vector.data, total_raw_points);
                    gsl_spline_init(spline_phi_interp, K_all_raw_v.vector.data, phi_unwrapped_all->data, total_raw_points);

                    double K_interp_val_start = gsl_vector_get(&K_all_raw_v.vector, 0);
                    double K_interp_val_end = K_final_reached_integration; 
                                        
                    if (K_interp_val_end < K_interp_val_start - EPSILON_GENERAL) { 
                         fprintf(stderr, "Warning C%d: K_interp_val_end (%.15e) < K_interp_val_start (%.15e) for spline. Using K_interp_val_start as end.\n", current_call_instance, K_interp_val_end, K_interp_val_start);
                         K_interp_val_end = K_interp_val_start;
                    }

                    for (size_t i = 0; i < full_traj_output->K->size; ++i) {
                        double K_val_to_eval = (full_traj_output->K->size > 1) ?
                            (K_interp_val_start + (double)i * (K_interp_val_end - K_interp_val_start) / (double)(full_traj_output->K->size - 1)) : K_interp_val_end;
                        
                        double K_spline_data_start = gsl_vector_get(K_all_raw, 0); 
                        double K_spline_data_end = gsl_vector_get(K_all_raw, total_raw_points - 1);
                        K_val_to_eval = fmax(K_spline_data_start, fmin(K_spline_data_end, K_val_to_eval)); 

                        double r_val = gsl_spline_eval(spline_r_interp, K_val_to_eval, acc_r_interp);
                        double phi_unwrapped_val = gsl_spline_eval(spline_phi_interp, K_val_to_eval, acc_phi_interp);
                        double phi_val_norm = normalize_phi(phi_unwrapped_val);

                        gsl_vector_set(full_traj_output->K, i, K_val_to_eval);
                        gsl_vector_set(full_traj_output->r, i, r_val);
                        gsl_vector_set(full_traj_output->phi, i, phi_val_norm);
                        gsl_vector_set(full_traj_output->x, i, r_val * cos(phi_val_norm));
                        gsl_vector_set(full_traj_output->y, i, r_val * sin(phi_val_norm));
                    }
                cleanup_core_interp_splines: 
                    if(spline_r_interp) { gsl_spline_free(spline_r_interp); spline_r_interp = NULL; }
                    if(spline_phi_interp) { gsl_spline_free(spline_phi_interp); spline_phi_interp = NULL; }
                    if(acc_r_interp) { gsl_interp_accel_free(acc_r_interp); acc_r_interp = NULL; }
                    if(acc_phi_interp) { gsl_interp_accel_free(acc_phi_interp); acc_phi_interp = NULL; }
                } 
            cleanup_core_interp_arrays: 
                if(K_all_raw) { gsl_vector_free(K_all_raw); K_all_raw = NULL; }
                if(r_all_raw) { gsl_vector_free(r_all_raw); r_all_raw = NULL; }
                if(phi_all_raw_temp) { gsl_vector_free(phi_all_raw_temp); phi_all_raw_temp = NULL; }
                if(phi_unwrapped_all) { gsl_vector_free(phi_unwrapped_all); phi_unwrapped_all = NULL; }
                 // If integration_stop_code was set to 3 by a CHECK_ macro in this block,
                 // full_traj_output->K might be non-NULL but other vectors might be NULL.
                 // The error_code setting later will handle this.
            }  else { // total_raw_points < 1
                fprintf(stderr, "Warning C%d: No raw points available for interpolation. Output K will be NULL.\n", current_call_instance);
                if (full_traj_output->K) { gsl_vector_free(full_traj_output->K); full_traj_output->K = NULL;}
                if (full_traj_output->r) { gsl_vector_free(full_traj_output->r); full_traj_output->r = NULL;}
                if (full_traj_output->phi) { gsl_vector_free(full_traj_output->phi); full_traj_output->phi = NULL;}
                if (full_traj_output->x) { gsl_vector_free(full_traj_output->x); full_traj_output->x = NULL;}
                if (full_traj_output->y) { gsl_vector_free(full_traj_output->y); full_traj_output->y = NULL;}
            }
        cleanup_core_interp:; 
        } else if (integration_stop_code != 3 && integration_stop_code != 4 && integration_stop_code != 5) {
            // This branch means: (num_segments_collected == 0 AND not the initial trivial case) OR K_final_reached_integration was negative.
            // Implies something went wrong very early, or no valid segments.
            fprintf(stderr, "Warning C%d: No segments collected or K_final invalid for interpolation. K_final=%.3e. Output K will be NULL.\n", current_call_instance, K_final_reached_integration);
            if (full_traj_output && full_traj_output->K) { gsl_vector_free(full_traj_output->K); full_traj_output->K = NULL;}
            if (full_traj_output && full_traj_output->r) { gsl_vector_free(full_traj_output->r); full_traj_output->r = NULL;}
            if (full_traj_output && full_traj_output->phi) { gsl_vector_free(full_traj_output->phi); full_traj_output->phi = NULL;}
            if (full_traj_output && full_traj_output->x) { gsl_vector_free(full_traj_output->x); full_traj_output->x = NULL;}
            if (full_traj_output && full_traj_output->y) { gsl_vector_free(full_traj_output->y); full_traj_output->y = NULL;}
        }
        
        if (full_traj_output) { 
            if (full_traj_output->crossings_y_at_x_targets) { 
                for (size_t i = 0; i < full_traj_output->num_x_targets; ++i) {
                    if(full_traj_output->crossings_y_at_x_targets[i]) {
                        gsl_vector_free(full_traj_output->crossings_y_at_x_targets[i]);
                    }
                }
                free(full_traj_output->crossings_y_at_x_targets);
            }
            full_traj_output->crossings_y_at_x_targets = crossings_collector_ptr;
            full_traj_output->num_x_targets = num_x_targets_val;
            crossings_collector_ptr = NULL; 
            if(crossings_collector_capacities) { free(crossings_collector_capacities); crossings_collector_capacities = NULL; }
        }
        
        if (full_traj_output && full_traj_output->error_code == 0 && (integration_stop_code >= 3 && integration_stop_code <= 5) ) { 
            full_traj_output->error_code = -2; 
        }
    } // end if (full_traj_output && integration_stop_code not critical error)

    if (crossings_output) {
        if (crossings_collector_ptr) { 
            if (crossings_output->crossings_y_at_x_targets) { 
                for (size_t i = 0; i < crossings_output->num_x_targets; ++i) {
                    if(crossings_output->crossings_y_at_x_targets[i]) {
                        gsl_vector_free(crossings_output->crossings_y_at_x_targets[i]);
                    }
                }
                free(crossings_output->crossings_y_at_x_targets);
            }
            crossings_output->crossings_y_at_x_targets = crossings_collector_ptr;
            crossings_output->num_x_targets = num_x_targets_val;
            crossings_collector_ptr = NULL; 
            if(crossings_collector_capacities) { free(crossings_collector_capacities); crossings_collector_capacities = NULL; }
        }
        crossings_output->error_code = (integration_stop_code >= 3 && integration_stop_code <=5) ? -2 : 0;
    }

cleanup_core:
    if(driver) gsl_odeiv2_driver_free(driver);
    if(current_segment_K_temp) gsl_vector_free(current_segment_K_temp);
    if(current_segment_r_temp) gsl_vector_free(current_segment_r_temp);
    if(current_segment_phi_temp) gsl_vector_free(current_segment_phi_temp);
    if(segments_collected_list){
        for(size_t i=0; i<num_segments_collected; ++i){
            if(segments_collected_list[i].K_pts) gsl_vector_free(segments_collected_list[i].K_pts);
            if(segments_collected_list[i].r_pts) gsl_vector_free(segments_collected_list[i].r_pts);
            if(segments_collected_list[i].phi_pts) gsl_vector_free(segments_collected_list[i].phi_pts);
        }
        free(segments_collected_list);
    }
    if(crossings_collector_ptr){
        for(size_t i=0; i<num_x_targets_val; ++i) {
            if(crossings_collector_ptr[i]) gsl_vector_free(crossings_collector_ptr[i]);
        }
        free(crossings_collector_ptr);
    }
    if(crossings_collector_capacities) free(crossings_collector_capacities);
    if(K_all_raw) gsl_vector_free(K_all_raw);
    if(r_all_raw) gsl_vector_free(r_all_raw);
    if(phi_all_raw_temp) gsl_vector_free(phi_all_raw_temp);
    if(phi_unwrapped_all) gsl_vector_free(phi_unwrapped_all);
    if(spline_r_interp) gsl_spline_free(spline_r_interp);
    if(spline_phi_interp) gsl_spline_free(spline_phi_interp);
    if(acc_r_interp) gsl_interp_accel_free(acc_r_interp);
    if(acc_phi_interp) gsl_interp_accel_free(acc_phi_interp);

    // Propagate errors
    if (full_traj_output && full_traj_output->error_code == 0 && core_error_code != 0) full_traj_output->error_code = core_error_code;
    if (crossings_output && crossings_output->error_code == 0 && core_error_code != 0) crossings_output->error_code = core_error_code;
    
    // If an integration error happened (not alloc error), it's in integration_stop_code
    if (full_traj_output && full_traj_output->error_code == 0 && (integration_stop_code >= 3 && integration_stop_code <= 5)) full_traj_output->error_code = -2; // General integration error
    if (crossings_output && crossings_output->error_code == 0 && (integration_stop_code >= 3 && integration_stop_code <= 5)) crossings_output->error_code = -2;

    if (core_error_code !=0) return core_error_code; // Prioritize allocation errors
    return (integration_stop_code >= 3 && integration_stop_code <= 5) ? -1 : 0; // Return -1 for integration failures
}


// --- Public API Functions ---

PhotonTrajectory* compute_trajectory(
    double r_0, double phi_0, double M, double psi, double r_max,
    double x_stop_val, bool x_stop_active,
    const gsl_vector *x_targets,
    double t_end, int num_interp_points
) {
    PhotonTrajectory *result = calloc(1, sizeof(PhotonTrajectory));
    if (!result) { 
        perror("compute_trajectory: calloc PhotonTrajectory failed"); 
        return NULL; 
    }
    result->error_code = 0;
    result->K = NULL; result->r = NULL; result->phi = NULL; 
    result->x = NULL; result->y = NULL;
    result->crossings_y_at_x_targets = NULL; 
    result->num_x_targets = 0;

    int n_interp = (num_interp_points > 0) ? num_interp_points : DEFAULT_NUM_INTERP_POINTS;
    if (num_interp_points == 0) n_interp = 1; 

    int core_status = integrate_photon_trajectory_core(
        r_0, phi_0, M, psi, r_max, x_stop_val, x_stop_active, x_targets, t_end,
        1e-10, 1e-10, result, NULL, n_interp);

    if (core_status != 0 && result->error_code == 0) {
        result->error_code = core_status; 
    }
    if (result->K == NULL && result->error_code == 0 && n_interp > 0) {
        fprintf(stderr, "compute_trajectory: Warning - K is NULL post-integration despite no error_code set by core.\n");
        result->error_code = -10; 
    }
    return result;
}

void free_photon_trajectory(PhotonTrajectory *traj) {
    if (!traj) { return; }
    if(traj->K) { gsl_vector_free(traj->K); traj->K = NULL; }
    if(traj->r) { gsl_vector_free(traj->r); traj->r = NULL; }
    if(traj->phi) { gsl_vector_free(traj->phi); traj->phi = NULL; }
    if(traj->x) { gsl_vector_free(traj->x); traj->x = NULL; }
    if(traj->y) { gsl_vector_free(traj->y); traj->y = NULL; }
    if (traj->crossings_y_at_x_targets) {
        for (size_t i = 0; i < traj->num_x_targets; ++i) {
            if(traj->crossings_y_at_x_targets[i]) { 
                gsl_vector_free(traj->crossings_y_at_x_targets[i]);
            }
        }
        free(traj->crossings_y_at_x_targets);
        traj->crossings_y_at_x_targets = NULL;
    }
    free(traj);
}

TrajectoryCrossings* compute_trajectory_crossings_only(
    double r_0, double phi_0, double M, double psi, double r_max,
    double x_stop_val, bool x_stop_active,
    const gsl_vector *x_targets, double t_end
) {
    if (!x_targets) { 
        fprintf(stderr, "compute_trajectory_crossings_only: x_targets cannot be NULL.\n");
        return NULL;
    }
    TrajectoryCrossings *result = calloc(1, sizeof(TrajectoryCrossings));
    if (!result) { 
        perror("compute_trajectory_crossings_only: calloc TrajectoryCrossings failed"); 
        return NULL; 
    }
    result->error_code = 0;
    result->crossings_y_at_x_targets = NULL; 
    result->num_x_targets = 0; 

    int core_status = integrate_photon_trajectory_core(
        r_0, phi_0, M, psi, r_max, x_stop_val, x_stop_active, x_targets, t_end,
        1e-8, 1e-10, NULL, result, 0 ); 

    if (core_status != 0 && result->error_code == 0) {
        result->error_code = core_status;
    }
     if (result->crossings_y_at_x_targets == NULL && result->error_code == 0 && x_targets->size > 0) {
        fprintf(stderr, "compute_trajectory_crossings_only: Warning - crossings_y_at_x_targets is NULL post-integration.\n");
        result->error_code = -11; 
    }
    return result;
}

void free_trajectory_crossings(TrajectoryCrossings *crossings) {
    if (!crossings) { return; }
    if (crossings->crossings_y_at_x_targets) {
        for (size_t i = 0; i < crossings->num_x_targets; ++i) {
            if(crossings->crossings_y_at_x_targets[i]) { 
                gsl_vector_free(crossings->crossings_y_at_x_targets[i]);
            }
        }
        free(crossings->crossings_y_at_x_targets);
        crossings->crossings_y_at_x_targets = NULL;
    }
    free(crossings);
}


#ifdef UNIT_TEST_TRAJECTORY_MODULE
// Example unit tests for compute_trajectory and compute_trajectory_crossings_only
// Ensure M_PI is available (from math.h, typically included via schwarzschild_tracer.h)
int main() {
    printf("--- Unit Test for Trajectory Module (compute_trajectory*, free_*) ---\n");
    gsl_set_error_handler_off();

    // Test 1: compute_trajectory - scattering
    printf("Test 1: compute_trajectory (scattering)...\n");
    PhotonTrajectory *traj1 = compute_trajectory(20.0, 0.0, 1.0, -1.25, 100.0, 
                                                 gsl_nan(), false, // x_stop_val, x_stop_active
                                                 NULL, // x_targets
                                                 1000.0, 10); // t_end, num_interp_points
    if (traj1 && traj1->error_code == 0 && traj1->K && traj1->K->size > 0 && traj1->r && traj1->r->size == traj1->K->size) {
        printf("  Scatter test K size: %zu, final r: %g\n", traj1->K->size, gsl_vector_get(traj1->r, traj1->r->size-1));
        if (gsl_vector_get(traj1->r, traj1->r->size-1) > 19.0) { 
            printf("  Scatter test PASSED.\n");
        } else {
            printf("  Scatter test FAILED (final r %.2f too small).\n", gsl_vector_get(traj1->r, traj1->r->size-1));
        }
    } else {
        printf("  Scatter test FAILED (error_code %d or NULL/mismatched K/r or K empty).\n", traj1 ? traj1->error_code : -999);
    }
    free_photon_trajectory(traj1);

    // Test 2: compute_trajectory_crossings_only - inward radial
    printf("Test 2: compute_trajectory_crossings_only (inward radial)...\n");
    gsl_vector *xt_test2 = gsl_vector_alloc(2);
    gsl_vector_set(xt_test2, 0, 4.0);
    gsl_vector_set(xt_test2, 1, 3.0);
    TrajectoryCrossings *cross2 = compute_trajectory_crossings_only(5.0, 0.0, 1.0, -M_PI/2, 100.0, 
                                                                    0.0, true, // x_stop_val, x_stop_active
                                                                    xt_test2, 1000.0); 
    if (cross2 && cross2->error_code == 0 && cross2->num_x_targets == 2 && 
        cross2->crossings_y_at_x_targets && 
        cross2->crossings_y_at_x_targets[0] && cross2->crossings_y_at_x_targets[0]->size > 0 &&
        cross2->crossings_y_at_x_targets[1] && cross2->crossings_y_at_x_targets[1]->size > 0) {
        printf("  Crossings test: Target x=%.1f, y_cross=%.4f; Target x=%.1f, y_cross=%.4f\n",
                gsl_vector_get(xt_test2, 0), gsl_vector_get(cross2->crossings_y_at_x_targets[0],0),
                gsl_vector_get(xt_test2, 1), gsl_vector_get(cross2->crossings_y_at_x_targets[1],0));
        if (fabs(gsl_vector_get(cross2->crossings_y_at_x_targets[0],0)) < 1e-6 &&
            fabs(gsl_vector_get(cross2->crossings_y_at_x_targets[1],0)) < 1e-6) {
             printf("  Crossings test PASSED (y ~ 0 for radial).\n");
        } else {
             printf("  Crossings test FAILED (y not close to 0).\n");
        }
    } else {
        printf("  Crossings test FAILED (error_code %d or unexpected structure/counts).\n", cross2 ? cross2->error_code : -999);
         if (cross2) {
            printf("    Details: num_x_targets=%zu, crossings_ptr=%p\n", cross2->num_x_targets, (void*)cross2->crossings_y_at_x_targets);
            if (cross2->crossings_y_at_x_targets && cross2->num_x_targets >0) {
                 printf("      Target 0: ptr=%p, size=%zu\n", (void*)cross2->crossings_y_at_x_targets[0], cross2->crossings_y_at_x_targets[0] ? cross2->crossings_y_at_x_targets[0]->size : 0);
                 if (cross2->num_x_targets >1) printf("      Target 1: ptr=%p, size=%zu\n", (void*)cross2->crossings_y_at_x_targets[1], cross2->crossings_y_at_x_targets[1] ? cross2->crossings_y_at_x_targets[1]->size : 0);
            }
        }
    }
    free_trajectory_crossings(cross2);
    gsl_vector_free(xt_test2);

    printf("--- Unit Test for Trajectory Module Finished ---\n");
    return 0;
}
#endif // UNIT_TEST_TRAJECTORY_MODULE

// trajectory_module.c
// Bottom
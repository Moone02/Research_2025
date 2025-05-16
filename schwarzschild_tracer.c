// schwarzschild_tracer.c
#define _POSIX_C_SOURCE 200809L // For strdup

#include "schwarzschild_tracer.h"

// --- Helper Macros ---
#define CHECK_GSL_ALLOC_VEC_CORE(ptr, func_name, item_name_str, err_var, err_val, cleanup_label) \
    if ((ptr) == NULL) { fprintf(stderr, "%s: GSL vector allocation failed for %s\n", func_name, item_name_str); err_var = err_val; goto cleanup_label; }
#define CHECK_ALLOC_GEN_CORE(ptr, func_name, type_name_str, err_var, err_val, cleanup_label) \
    if ((ptr) == NULL) { fprintf(stderr, "%s: Standard allocation failed for %s\n", func_name, type_name_str); err_var = err_val; goto cleanup_label; }

// --- Struct Definition for GSL Root Finder Wrapper (Moved Up) ---
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
    if (!s_solver) { fprintf(stderr, "find_event_time_gsl: Failed to allocate GSL root solver.\n"); return 0; }

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

static int reallocate_gsl_vector_if_needed(gsl_vector **vec_ptr, size_t current_count, size_t *current_capacity_ptr, size_t initial_capacity_val) {
    if (!vec_ptr || !current_capacity_ptr) { return -1; }

    if (current_count >= *current_capacity_ptr) {
        size_t new_capacity = (*current_capacity_ptr > 0) ? (*current_capacity_ptr * 2) : initial_capacity_val;
        if (new_capacity <= current_count) { 
            new_capacity = current_count + initial_capacity_val; 
            if (new_capacity <= current_count) { 
                new_capacity = current_count + 1;
            }
        }

        gsl_vector *new_v = gsl_vector_alloc(new_capacity);
        if (!new_v) {
            fprintf(stderr, "reallocate_gsl_vector_if_needed: Failed to allocate new vector of size %zu\n", new_capacity);
            return -1; 
        }

        gsl_vector *old_v = *vec_ptr;
        if (old_v && old_v->size > 0) { 
            size_t copy_size = old_v->size; 
            if (copy_size > 0) {
                 gsl_vector_const_view old_sub = gsl_vector_const_subvector(old_v, 0, copy_size);
                 gsl_vector_view new_sub = gsl_vector_subvector(new_v, 0, copy_size);
                 gsl_vector_memcpy(&new_sub.vector, &old_sub.vector);
            }
            new_v->size = old_v->size; 
        } else {
            new_v->size = 0; 
        }

        gsl_vector_free(old_v);
        *vec_ptr = new_v;
        *current_capacity_ptr = new_capacity;
    }
    return 0;
}


static int gsl_vector_dynamic_append(gsl_vector **vec_ptr, double value, size_t *current_capacity_ptr) {
    if (!vec_ptr || !current_capacity_ptr) { return -1; }
    gsl_vector *vec = *vec_ptr;

    if (!vec) { 
        *current_capacity_ptr = INITIAL_CROSSING_POINTS_CAPACITY; 
        vec = gsl_vector_alloc(*current_capacity_ptr);
        if (!vec) { fprintf(stderr, "gsl_vector_dynamic_append: Initial allocation failed.\n"); return -1; }
        vec->size = 0; 
    }
    
    if (reallocate_gsl_vector_if_needed(&vec, vec->size, current_capacity_ptr, INITIAL_CROSSING_POINTS_CAPACITY) != 0) {
        return -1; 
    }
    
    gsl_vector_set(vec, vec->size, value);
    vec->size++;
    *vec_ptr = vec;
    return 0;
}

static void format_coord_str_for_filename(char *buffer, size_t buffer_size, double val) {
    char temp[128];
    snprintf(temp, sizeof(temp), "%.6g", val);
    size_t j = 0;
    for (size_t i = 0; temp[i] != '\0' && j < buffer_size - 1; ++i) {
        if (temp[i] == '.') {
            if (j < buffer_size - 1) { buffer[j++] = 'p'; }
        } else if (temp[i] == '-') {
            if (j + 3 < buffer_size - 1) {
                buffer[j++] = 'n'; buffer[j++] = 'e'; buffer[j++] = 'g';
            }
        } else {
            if (j < buffer_size - 1) { buffer[j++] = temp[i]; }
        }
    }
    buffer[j] = '\0';
}

static double normalize_phi(double phi) {
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

static int add_string_to_list(char*** list_ptr, int* count_ptr, int* capacity_ptr, const char* str_to_add) {
    if (!list_ptr || !count_ptr || !capacity_ptr) { return -1; }

    if (*count_ptr + 1 >= *capacity_ptr ) { // Need space for current string + NULL terminator
        int new_capacity = (*capacity_ptr == 0) ? INITIAL_SAVED_FILES_CAPACITY : *capacity_ptr * 2;
        if (new_capacity <= *count_ptr + 1) { new_capacity = *count_ptr + 1 + INITIAL_SAVED_FILES_CAPACITY; }


        char** temp_list = realloc(*list_ptr, new_capacity * sizeof(char*));
        if (!temp_list) {
            fprintf(stderr, "add_string_to_list: Failed to realloc file list to capacity %d.\n", new_capacity);
            return -1;
        }
        *list_ptr = temp_list;
        *capacity_ptr = new_capacity;
    }

    if (str_to_add != NULL) {
        (*list_ptr)[*count_ptr] = strdup(str_to_add);
        if (!(*list_ptr)[*count_ptr]) {
            fprintf(stderr, "add_string_to_list: Failed to strdup filepath.\n");
            return -1;
        }
        (*count_ptr)++; 
        (*list_ptr)[*count_ptr] = NULL; 
    } else { 
        (*list_ptr)[*count_ptr] = NULL;
    }
    return 0;
}


static void skip_ppm_comments_and_whitespace(FILE *fp) {
    int ch;
    while ((ch = fgetc(fp)) != EOF && isspace(ch)) {}
    if (ch == '#') {
        while ((ch = fgetc(fp)) != EOF && ch != '\n' && ch != '\r') {}
        skip_ppm_comments_and_whitespace(fp);
    } else if (ch != EOF) {
        ungetc(ch, fp);
    }
}

// --- Core Trajectory Computation (integrate_photon_trajectory_core) ---
static int integrate_photon_trajectory_core(
    double r_0, double phi_0, double M_val, double psi, double r_max_val,
    double x_stop_val, bool x_stop_active_flag,
    const gsl_vector *x_targets_vec,
    double t_end_max,
    double rtol, double atol,
    PhotonTrajectory *full_traj_output, 
    TrajectoryCrossings *crossings_output,
    int num_interp_points_for_full_traj
) {
    int core_error_code = 0;
    double K_final_reached_integration = 0.0;

    gsl_odeiv2_driver *driver = NULL;
    TrajectorySegmentDataInternal *segments_collected_list = NULL;
    size_t num_segments_collected = 0;
    size_t segments_collected_capacity = 0;
    gsl_vector *current_segment_K_temp = NULL, *current_segment_r_temp = NULL, *current_segment_phi_temp = NULL;
    size_t current_segment_point_count = 0;
    size_t current_segment_point_capacity = 0;
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
                CHECK_ALLOC_GEN_CORE(full_traj_output->crossings_y_at_x_targets, "core (trivial crossings_y_at_x_targets)", "gsl_vector**", core_error_code, -1, cleanup_core);
                for (size_t i = 0; i < full_traj_output->num_x_targets; ++i) {
                    full_traj_output->crossings_y_at_x_targets[i] = gsl_vector_alloc(0);
                    CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->crossings_y_at_x_targets[i], "core (trivial crossings[i])", "crossings_y_at_x_targets[i]", core_error_code, -1, cleanup_core);
                }
            } else { full_traj_output->crossings_y_at_x_targets = NULL; }
        }
        if (crossings_output) {
            crossings_output->num_x_targets = x_targets_vec ? x_targets_vec->size : 0;
            if (crossings_output->num_x_targets > 0) {
                crossings_output->crossings_y_at_x_targets = calloc(crossings_output->num_x_targets, sizeof(gsl_vector*));
                CHECK_ALLOC_GEN_CORE(crossings_output->crossings_y_at_x_targets, "core (trivial crossings_output array)", "gsl_vector**", core_error_code, -1, cleanup_core);
                for (size_t i = 0; i < crossings_output->num_x_targets; ++i) {
                    crossings_output->crossings_y_at_x_targets[i] = gsl_vector_alloc(0);
                    CHECK_GSL_ALLOC_VEC_CORE(crossings_output->crossings_y_at_x_targets[i], "core (trivial crossings_output[i])", "crossings_y_at_x_targets[i]", core_error_code, -1, cleanup_core);
                }
            } else { crossings_output->crossings_y_at_x_targets = NULL; }
        }
        goto cleanup_core;
    }

    double metric_r0_term = (1.0 - 2.0 * M_val / r_0);
    if (metric_r0_term <= EPSILON_GENERAL) { K_final_reached_integration = 0.0; goto cleanup_core; }

    double b_val = cos(psi) * r_0 / sqrt(metric_r0_term);
    int current_sign_dr_dk = 1;
    double psi_norm = normalize_phi(psi);
    if (psi_norm > M_PI && psi_norm < 2.0 * M_PI) { current_sign_dr_dk = -1; }

    ODEParams ode_params_instance;
    ode_params_instance.M = M_val;
    ode_params_instance.b = b_val;
    ode_params_instance.sign_dr_dk = &current_sign_dr_dk;

    gsl_odeiv2_system sys = {schwarzschild_geodesic_eqs, NULL, 2, &ode_params_instance};
    
    driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-6, rtol, atol);
    CHECK_ALLOC_GEN_CORE(driver, "core (driver)", "gsl_odeiv2_driver", core_error_code, -2, cleanup_core);

    double K_current = 0.0;
    double y_current_state[2] = {r_0, phi_0};

    if (full_traj_output) {
        segments_collected_capacity = INITIAL_SEGMENTS_CAPACITY;
        segments_collected_list = malloc(segments_collected_capacity * sizeof(TrajectorySegmentDataInternal));
        CHECK_ALLOC_GEN_CORE(segments_collected_list, "core (segments_collected_list)", "TrajectorySegmentDataInternal array", core_error_code, -1, cleanup_core);
        for(size_t i=0; i<segments_collected_capacity; ++i) {
            segments_collected_list[i].K_pts = NULL; segments_collected_list[i].r_pts = NULL; segments_collected_list[i].phi_pts = NULL;
        }
        current_segment_point_capacity = INITIAL_RAW_POINTS_CAPACITY;
        current_segment_K_temp = gsl_vector_alloc(current_segment_point_capacity); CHECK_GSL_ALLOC_VEC_CORE(current_segment_K_temp, "core (K_temp)", "current_segment_K_temp", core_error_code, -1, cleanup_core);
        current_segment_r_temp = gsl_vector_alloc(current_segment_point_capacity); CHECK_GSL_ALLOC_VEC_CORE(current_segment_r_temp, "core (r_temp)", "current_segment_r_temp", core_error_code, -1, cleanup_core);
        current_segment_phi_temp = gsl_vector_alloc(current_segment_point_capacity); CHECK_GSL_ALLOC_VEC_CORE(current_segment_phi_temp, "core (phi_temp)", "current_segment_phi_temp", core_error_code, -1, cleanup_core);
        
        gsl_vector_set(current_segment_K_temp, 0, K_current);
        gsl_vector_set(current_segment_r_temp, 0, y_current_state[0]);
        gsl_vector_set(current_segment_phi_temp, 0, y_current_state[1]);
        current_segment_point_count = 1;
    }

    if (x_targets_vec && (full_traj_output || crossings_output)) {
        num_x_targets_val = x_targets_vec->size;
        if (num_x_targets_val > 0) {
            crossings_collector_ptr = calloc(num_x_targets_val, sizeof(gsl_vector*)); CHECK_ALLOC_GEN_CORE(crossings_collector_ptr, "core (crossings_collector_ptr)", "gsl_vector**", core_error_code, -1, cleanup_core);
            crossings_collector_capacities = calloc(num_x_targets_val, sizeof(size_t)); CHECK_ALLOC_GEN_CORE(crossings_collector_capacities, "core (crossings_collector_capacities)", "size_t*", core_error_code, -1, cleanup_core);
            for (size_t i = 0; i < num_x_targets_val; ++i) {
                crossings_collector_ptr[i] = gsl_vector_alloc(INITIAL_CROSSING_POINTS_CAPACITY); CHECK_GSL_ALLOC_VEC_CORE(crossings_collector_ptr[i], "core (crossings_collector_ptr[i])", "crossings_collector_ptr[i]", core_error_code, -1, cleanup_core);
                gsl_vector_set_zero(crossings_collector_ptr[i]);
                crossings_collector_ptr[i]->size = 0;
                crossings_collector_capacities[i] = INITIAL_CROSSING_POINTS_CAPACITY;
            }
        }
    }

    #define MAX_EVENT_TYPES_CORE 100
    EventFunctionParams event_params_instances[MAX_EVENT_TYPES_CORE];
    double (*event_get_val_func_ptrs[MAX_EVENT_TYPES_CORE])(double, const double[], EventFunctionParams*);
    bool event_is_terminal[MAX_EVENT_TYPES_CORE];
    int event_crossing_dir[MAX_EVENT_TYPES_CORE];
    int event_is_x_target_idx_map[MAX_EVENT_TYPES_CORE];
    size_t num_active_event_funcs = 0;

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
            if (num_active_event_funcs >= MAX_EVENT_TYPES_CORE) { core_error_code = -2; goto cleanup_core; }
            event_params_instances[num_active_event_funcs].x_target_event = gsl_vector_get(x_targets_vec, i);
            event_get_val_func_ptrs[num_active_event_funcs] = get_event_val_x_target;
            event_is_terminal[num_active_event_funcs] = false;
            event_crossing_dir[num_active_event_funcs] = -1;
            event_is_x_target_idx_map[num_active_event_funcs] = i; num_active_event_funcs++;
        }
    }

    int integration_stop_code = 0;
    int safety_break_counter = 0;
    const int MAX_SAFETY_BREAK_CORE = 200000;

    while (K_current < t_end_max && integration_stop_code == 0 && safety_break_counter < MAX_SAFETY_BREAK_CORE) {
        safety_break_counter++;
        double K_before_this_step = K_current;
        double y_before_this_step[2]; memcpy(y_before_this_step, y_current_state, 2 * sizeof(double));

        int status = gsl_odeiv2_driver_apply(driver, &K_current, t_end_max, y_current_state);

        if (status != GSL_SUCCESS) {
            fprintf(stderr, "GSL ODE solver error: %s (K=%.2e, r=%.2e, phi=%.2e)\n", gsl_strerror(status), K_before_this_step, y_before_this_step[0], y_before_this_step[1]);
            integration_stop_code = 3; K_final_reached_integration = K_before_this_step; break;
        }
        K_final_reached_integration = K_current;

        if (full_traj_output) {
            if (reallocate_gsl_vector_if_needed(&current_segment_K_temp, current_segment_point_count, &current_segment_point_capacity, INITIAL_RAW_POINTS_CAPACITY) != 0 ||
                reallocate_gsl_vector_if_needed(&current_segment_r_temp, current_segment_point_count, &current_segment_point_capacity, INITIAL_RAW_POINTS_CAPACITY) != 0 ||
                reallocate_gsl_vector_if_needed(&current_segment_phi_temp, current_segment_point_count, &current_segment_point_capacity, INITIAL_RAW_POINTS_CAPACITY) != 0) {
                integration_stop_code = 3; K_final_reached_integration = K_before_this_step; break;
            }
            gsl_vector_set(current_segment_K_temp, current_segment_point_count, K_current);
            gsl_vector_set(current_segment_r_temp, current_segment_point_count, y_current_state[0]);
            gsl_vector_set(current_segment_phi_temp, current_segment_point_count, y_current_state[1]);
            current_segment_point_count++;
        }

        double earliest_K_event_this_step = K_current + 1.0;
        int triggered_event_func_original_idx = -1;
        double y_at_earliest_event_this_step[2];
        y_at_earliest_event_this_step[0] = NAN; y_at_earliest_event_this_step[1] = NAN; // Initialize

        for (size_t ev_idx = 0; ev_idx < num_active_event_funcs; ++ev_idx) {
            double K_event_found_type; double y_event_found_type[2];
            int event_found_flag = find_event_time_gsl(
                K_before_this_step, y_before_this_step, K_current, y_current_state,
                event_get_val_func_ptrs[ev_idx], &event_params_instances[ev_idx],
                event_crossing_dir[ev_idx], &K_event_found_type, y_event_found_type
            );

            if (event_found_flag && K_event_found_type >= K_before_this_step - EVENT_DETECTION_TOLERANCE && K_event_found_type <= K_current + EVENT_DETECTION_TOLERANCE) {
                if (K_event_found_type < earliest_K_event_this_step) {
                    earliest_K_event_this_step = K_event_found_type;
                    triggered_event_func_original_idx = ev_idx;
                    memcpy(y_at_earliest_event_this_step, y_event_found_type, 2 * sizeof(double));
                } else if (fabs(K_event_found_type - earliest_K_event_this_step) < EVENT_DETECTION_TOLERANCE * 0.1) {
                    if (event_get_val_func_ptrs[ev_idx] == get_event_val_fr_zero) {
                        earliest_K_event_this_step = K_event_found_type; triggered_event_func_original_idx = ev_idx; memcpy(y_at_earliest_event_this_step, y_event_found_type, 2 * sizeof(double));
                    } else if (event_get_val_func_ptrs[ev_idx] == get_event_val_r_leq_2M && event_get_val_func_ptrs[triggered_event_func_original_idx] != get_event_val_fr_zero) {
                        earliest_K_event_this_step = K_event_found_type; triggered_event_func_original_idx = ev_idx; memcpy(y_at_earliest_event_this_step, y_event_found_type, 2 * sizeof(double));
                    }
                }
            }
        }

        if (triggered_event_func_original_idx != -1) {
            K_final_reached_integration = earliest_K_event_this_step;
            K_current = earliest_K_event_this_step;
            memcpy(y_current_state, y_at_earliest_event_this_step, 2 * sizeof(double));

            if (full_traj_output && current_segment_point_count > 0) {
                size_t last_idx_in_temp = current_segment_point_count - 1;
                gsl_vector_set(current_segment_K_temp, last_idx_in_temp, K_current);
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
                            integration_stop_code = 3; break;
                        }
                    }
                }
            } else if (event_is_terminal[triggered_event_func_original_idx]) {
                if (event_get_val_func_ptrs[triggered_event_func_original_idx] == get_event_val_fr_zero) {
                    current_sign_dr_dk *= -1;
                    y_current_state[0] += 1e-7 * current_sign_dr_dk;
                    K_current += 1e-7;
                    K_final_reached_integration = K_current;

                    if (y_current_state[0] <= (2*M_val + EPSILON_GENERAL) || y_current_state[0] >= r_max_val ||
                        (x_stop_active_flag && (y_current_state[0]*cos(y_current_state[1]) <= x_stop_val + EPSILON_GENERAL))) {
                        integration_stop_code = 2; break;
                    }

                    if (full_traj_output) {
                        if (num_segments_collected >= segments_collected_capacity) {
                            segments_collected_capacity = (segments_collected_capacity == 0) ? INITIAL_SEGMENTS_CAPACITY : segments_collected_capacity * 2;
                            TrajectorySegmentDataInternal *temp_realloc = realloc(segments_collected_list, segments_collected_capacity * sizeof(TrajectorySegmentDataInternal));
                            CHECK_ALLOC_GEN_CORE(temp_realloc, "core (realloc segments_list)", "TrajectorySegmentDataInternal array", integration_stop_code, 3, cleanup_core_loop);
                            segments_collected_list = temp_realloc;
                        }
                        segments_collected_list[num_segments_collected].K_pts = gsl_vector_alloc(current_segment_point_count); CHECK_GSL_ALLOC_VEC_CORE(segments_collected_list[num_segments_collected].K_pts, "core (seg K)", "segment K_pts", integration_stop_code, 3, cleanup_core_loop);
                        segments_collected_list[num_segments_collected].r_pts = gsl_vector_alloc(current_segment_point_count); CHECK_GSL_ALLOC_VEC_CORE(segments_collected_list[num_segments_collected].r_pts, "core (seg r)", "segment r_pts", integration_stop_code, 3, cleanup_core_loop);
                        segments_collected_list[num_segments_collected].phi_pts = gsl_vector_alloc(current_segment_point_count); CHECK_GSL_ALLOC_VEC_CORE(segments_collected_list[num_segments_collected].phi_pts, "core (seg phi)", "segment phi_pts", integration_stop_code, 3, cleanup_core_loop);
                        
                        for(size_t k=0; k<current_segment_point_count; ++k) {
                            gsl_vector_set(segments_collected_list[num_segments_collected].K_pts, k, gsl_vector_get(current_segment_K_temp, k));
                            gsl_vector_set(segments_collected_list[num_segments_collected].r_pts, k, gsl_vector_get(current_segment_r_temp, k));
                            gsl_vector_set(segments_collected_list[num_segments_collected].phi_pts, k, gsl_vector_get(current_segment_phi_temp, k));
                        }
                        num_segments_collected++;
                        current_segment_point_count = 0; 
                        if (reallocate_gsl_vector_if_needed(&current_segment_K_temp, current_segment_point_count, &current_segment_point_capacity, INITIAL_RAW_POINTS_CAPACITY) != 0 ||
                            reallocate_gsl_vector_if_needed(&current_segment_r_temp, current_segment_point_count, &current_segment_point_capacity, INITIAL_RAW_POINTS_CAPACITY) != 0 ||
                            reallocate_gsl_vector_if_needed(&current_segment_phi_temp, current_segment_point_count, &current_segment_point_capacity, INITIAL_RAW_POINTS_CAPACITY) != 0) {
                            integration_stop_code = 3; K_final_reached_integration = K_before_this_step; break;
                        }
                        gsl_vector_set(current_segment_K_temp, current_segment_point_count, K_current);
                        gsl_vector_set(current_segment_r_temp, current_segment_point_count, y_current_state[0]);
                        gsl_vector_set(current_segment_phi_temp, current_segment_point_count, y_current_state[1]);
                        current_segment_point_count = 1;
                    }
                } else {
                    integration_stop_code = 2; break;
                }
            }
        }
    cleanup_core_loop: {} 
        if (integration_stop_code == 3) { break; }

        if (K_current >= t_end_max - EPSILON_GENERAL) { integration_stop_code = 1; K_final_reached_integration = t_end_max; }
    }

    if (safety_break_counter >= MAX_SAFETY_BREAK_CORE) {
        fprintf(stderr, "Warning: Max safety break counter reached in integration loop.\n");
        if (integration_stop_code == 0) { integration_stop_code = 2; }
    }

    if (full_traj_output && current_segment_point_count > 0 && integration_stop_code != 3) {
        if (num_segments_collected >= segments_collected_capacity) {
            segments_collected_capacity = (segments_collected_capacity == 0) ? INITIAL_SEGMENTS_CAPACITY : segments_collected_capacity * 2;
            TrajectorySegmentDataInternal *temp_realloc = realloc(segments_collected_list, segments_collected_capacity * sizeof(TrajectorySegmentDataInternal));
            CHECK_ALLOC_GEN_CORE(temp_realloc, "core (realloc final segments_list)", "TrajectorySegmentDataInternal array", integration_stop_code, 3, cleanup_core_post_loop);
            segments_collected_list = temp_realloc;
        }
        segments_collected_list[num_segments_collected].K_pts = gsl_vector_alloc(current_segment_point_count); CHECK_GSL_ALLOC_VEC_CORE(segments_collected_list[num_segments_collected].K_pts, "core (final seg K)", "final segment K_pts", integration_stop_code, 3, cleanup_core_post_loop);
        segments_collected_list[num_segments_collected].r_pts = gsl_vector_alloc(current_segment_point_count); CHECK_GSL_ALLOC_VEC_CORE(segments_collected_list[num_segments_collected].r_pts, "core (final seg r)", "final segment r_pts", integration_stop_code, 3, cleanup_core_post_loop);
        segments_collected_list[num_segments_collected].phi_pts = gsl_vector_alloc(current_segment_point_count); CHECK_GSL_ALLOC_VEC_CORE(segments_collected_list[num_segments_collected].phi_pts, "core (final seg phi)", "final segment phi_pts", integration_stop_code, 3, cleanup_core_post_loop);
        
        for(size_t k=0; k<current_segment_point_count; ++k) {
            gsl_vector_set(segments_collected_list[num_segments_collected].K_pts, k, gsl_vector_get(current_segment_K_temp, k));
            gsl_vector_set(segments_collected_list[num_segments_collected].r_pts, k, gsl_vector_get(current_segment_r_temp, k));
            gsl_vector_set(segments_collected_list[num_segments_collected].phi_pts, k, gsl_vector_get(current_segment_phi_temp, k));
        }
        num_segments_collected++;
    }
cleanup_core_post_loop: {} 

    if (full_traj_output && integration_stop_code != 3) {
        if (num_segments_collected == 0 && fabs(K_final_reached_integration - 0.0) < DBL_EPSILON &&
            (full_traj_output->K && full_traj_output->K->size >= 1 && fabs(gsl_vector_get(full_traj_output->K,0) - 0.0) < DBL_EPSILON)) {
            // Trivial case already handled
        } else if (num_segments_collected > 0 && K_final_reached_integration >= -EPSILON_GENERAL) {
            size_t total_raw_points = 0;
            for (size_t s = 0; s < num_segments_collected; ++s) {
                total_raw_points += segments_collected_list[s].K_pts->size;
            }

            if (total_raw_points >= 1) {
                int n_pts_interp_actual = (num_interp_points_for_full_traj > 0) ? num_interp_points_for_full_traj : 1;
                if (full_traj_output->K == NULL || full_traj_output->K->size != (size_t)n_pts_interp_actual) {
                    gsl_vector_free(full_traj_output->K); gsl_vector_free(full_traj_output->r);
                    gsl_vector_free(full_traj_output->phi); gsl_vector_free(full_traj_output->x); gsl_vector_free(full_traj_output->y);
                    full_traj_output->K = gsl_vector_alloc(n_pts_interp_actual); CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->K, "core (interp K)", "interp K", integration_stop_code, 3, cleanup_core_interp);
                    full_traj_output->r = gsl_vector_alloc(n_pts_interp_actual); CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->r, "core (interp r)", "interp r", integration_stop_code, 3, cleanup_core_interp);
                    full_traj_output->phi = gsl_vector_alloc(n_pts_interp_actual); CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->phi, "core (interp phi)", "interp phi", integration_stop_code, 3, cleanup_core_interp);
                    full_traj_output->x = gsl_vector_alloc(n_pts_interp_actual); CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->x, "core (interp x)", "interp x", integration_stop_code, 3, cleanup_core_interp);
                    full_traj_output->y = gsl_vector_alloc(n_pts_interp_actual); CHECK_GSL_ALLOC_VEC_CORE(full_traj_output->y, "core (interp y)", "interp y", integration_stop_code, 3, cleanup_core_interp);
                }

                K_all_raw = gsl_vector_alloc(total_raw_points); CHECK_GSL_ALLOC_VEC_CORE(K_all_raw, "core (K_all_raw)", "K_all_raw", integration_stop_code, 3, cleanup_core_interp_arrays);
                r_all_raw = gsl_vector_alloc(total_raw_points); CHECK_GSL_ALLOC_VEC_CORE(r_all_raw, "core (r_all_raw)", "r_all_raw", integration_stop_code, 3, cleanup_core_interp_arrays);
                phi_all_raw_temp = gsl_vector_alloc(total_raw_points); CHECK_GSL_ALLOC_VEC_CORE(phi_all_raw_temp, "core (phi_all_raw_temp)", "phi_all_raw_temp", integration_stop_code, 3, cleanup_core_interp_arrays);

                size_t current_concat_idx = 0;
                for (size_t s = 0; s < num_segments_collected; ++s) {
                    for (size_t i = 0; i < segments_collected_list[s].K_pts->size; ++i) {
                        gsl_vector_set(K_all_raw, current_concat_idx, gsl_vector_get(segments_collected_list[s].K_pts, i));
                        gsl_vector_set(r_all_raw, current_concat_idx, gsl_vector_get(segments_collected_list[s].r_pts, i));
                        gsl_vector_set(phi_all_raw_temp, current_concat_idx, gsl_vector_get(segments_collected_list[s].phi_pts, i));
                        current_concat_idx++;
                    }
                }

                phi_unwrapped_all = gsl_vector_alloc(total_raw_points); CHECK_GSL_ALLOC_VEC_CORE(phi_unwrapped_all, "core (phi_unwrapped_all)", "phi_unwrapped_all", integration_stop_code, 3, cleanup_core_interp_arrays);
                unwrap_phi_values(K_all_raw, phi_all_raw_temp, phi_unwrapped_all);

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
                } else {
                    const gsl_interp_type *spline_type_cspline_ptr = gsl_interp_cspline;
                    const gsl_interp_type *spline_type_linear_ptr = gsl_interp_linear;
                    const gsl_interp_type *spline_type = (total_raw_points >= gsl_interp_type_min_size(spline_type_cspline_ptr)) ? spline_type_cspline_ptr : spline_type_linear_ptr;

                    spline_r_interp = gsl_spline_alloc(spline_type, total_raw_points); CHECK_ALLOC_GEN_CORE(spline_r_interp, "core (spline_r)", "gsl_spline", integration_stop_code, 3, cleanup_core_interp_splines);
                    spline_phi_interp = gsl_spline_alloc(spline_type, total_raw_points); CHECK_ALLOC_GEN_CORE(spline_phi_interp, "core (spline_phi)", "gsl_spline", integration_stop_code, 3, cleanup_core_interp_splines);
                    acc_r_interp = gsl_interp_accel_alloc(); CHECK_ALLOC_GEN_CORE(acc_r_interp, "core (acc_r)", "gsl_interp_accel", integration_stop_code, 3, cleanup_core_interp_splines);
                    acc_phi_interp = gsl_interp_accel_alloc(); CHECK_ALLOC_GEN_CORE(acc_phi_interp, "core (acc_phi)", "gsl_interp_accel", integration_stop_code, 3, cleanup_core_interp_splines);

                    gsl_spline_init(spline_r_interp, K_all_raw->data, r_all_raw->data, total_raw_points);
                    gsl_spline_init(spline_phi_interp, K_all_raw->data, phi_unwrapped_all->data, total_raw_points);

                    double K_interp_start = gsl_vector_get(K_all_raw, 0);
                    double K_interp_end = K_final_reached_integration;

                    for (size_t i = 0; i < full_traj_output->K->size; ++i) {
                        double K_val_interp = (full_traj_output->K->size > 1) ?
                            (K_interp_start + (double)i * (K_interp_end - K_interp_start) / (double)(full_traj_output->K->size - 1)) : K_interp_end;
                        if (K_interp_end < K_interp_start && full_traj_output->K->size > 1) { K_val_interp = K_interp_start; }
                        K_val_interp = fmax(K_interp_start, fmin(K_interp_end, K_val_interp));

                        double r_val = gsl_spline_eval(spline_r_interp, K_val_interp, acc_r_interp);
                        double phi_unwrapped_val = gsl_spline_eval(spline_phi_interp, K_val_interp, acc_phi_interp);
                        double phi_val_norm = normalize_phi(phi_unwrapped_val);

                        gsl_vector_set(full_traj_output->K, i, K_val_interp);
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
            }
        cleanup_core_interp:; 
        }
        full_traj_output->crossings_y_at_x_targets = crossings_collector_ptr;
        full_traj_output->num_x_targets = num_x_targets_val;
        crossings_collector_ptr = NULL; 
        if(crossings_collector_capacities) { free(crossings_collector_capacities); crossings_collector_capacities = NULL; }
        if (full_traj_output->error_code == 0 && integration_stop_code == 3) { full_traj_output->error_code = -2; }
    }

    if (crossings_output) {
        crossings_output->crossings_y_at_x_targets = crossings_collector_ptr;
        crossings_output->num_x_targets = num_x_targets_val;
        crossings_output->error_code = (integration_stop_code == 3) ? -2 : 0;
        crossings_collector_ptr = NULL; 
        if(crossings_collector_capacities) { free(crossings_collector_capacities); crossings_collector_capacities = NULL; }
    }

cleanup_core:
    if(driver) { gsl_odeiv2_driver_free(driver); driver = NULL; }
    if(current_segment_K_temp) { gsl_vector_free(current_segment_K_temp); current_segment_K_temp = NULL; }
    if(current_segment_r_temp) { gsl_vector_free(current_segment_r_temp); current_segment_r_temp = NULL; }
    if(current_segment_phi_temp) { gsl_vector_free(current_segment_phi_temp); current_segment_phi_temp = NULL; }
    if (segments_collected_list) {
        for (size_t i = 0; i < num_segments_collected; ++i) {
            if(segments_collected_list[i].K_pts) { gsl_vector_free(segments_collected_list[i].K_pts); }
            if(segments_collected_list[i].r_pts) { gsl_vector_free(segments_collected_list[i].r_pts); }
            if(segments_collected_list[i].phi_pts) { gsl_vector_free(segments_collected_list[i].phi_pts); }
        }
        free(segments_collected_list); segments_collected_list = NULL;
    }
    if (crossings_collector_ptr) { 
        for(size_t i=0; i < num_x_targets_val; ++i) {
            if (crossings_collector_ptr[i]) { gsl_vector_free(crossings_collector_ptr[i]); }
        }
        free(crossings_collector_ptr); crossings_collector_ptr = NULL;
    }
    if (crossings_collector_capacities) { free(crossings_collector_capacities); crossings_collector_capacities = NULL; }

    // Free any remaining GSL objects from interpolation if an error occurred mid-interpolation
    if(K_all_raw) gsl_vector_free(K_all_raw); if(r_all_raw) gsl_vector_free(r_all_raw);
    if(phi_all_raw_temp) gsl_vector_free(phi_all_raw_temp); if(phi_unwrapped_all) gsl_vector_free(phi_unwrapped_all);
    if(spline_r_interp) gsl_spline_free(spline_r_interp); if(spline_phi_interp) gsl_spline_free(spline_phi_interp);
    if(acc_r_interp) gsl_interp_accel_free(acc_r_interp); if(acc_phi_interp) gsl_interp_accel_free(acc_phi_interp);


    if (core_error_code != 0) {
        if (full_traj_output && full_traj_output->error_code == 0) { full_traj_output->error_code = core_error_code; }
        if (crossings_output && crossings_output->error_code == 0) { crossings_output->error_code = core_error_code; }
    }
    return (integration_stop_code == 3 || core_error_code != 0) ? -1 : 0;
}

// --- Public API wrappers for compute_trajectory and _crossings_only ---
PhotonTrajectory* compute_trajectory(
    double r_0, double phi_0, double M, double psi, double r_max,
    double x_stop_val, bool x_stop_active,
    const gsl_vector *x_targets,
    double t_end, int num_interp_points
) {
    PhotonTrajectory *result = calloc(1, sizeof(PhotonTrajectory));
    if (!result) { perror("compute_trajectory: calloc PhotonTrajectory failed"); return NULL; }
    result->error_code = 0;
    result->K = NULL; result->r = NULL; result->phi = NULL; result->x = NULL; result->y = NULL;
    result->crossings_y_at_x_targets = NULL; result->num_x_targets = 0;

    int n_interp = (num_interp_points > 0) ? num_interp_points : 1;

    int core_status = integrate_photon_trajectory_core(
        r_0, phi_0, M, psi, r_max, x_stop_val, x_stop_active, x_targets, t_end,
        1e-10, 1e-10, result, NULL, n_interp);

    if (core_status != 0 && result->error_code == 0) {
        result->error_code = core_status;
    }
    if (result->K == NULL && result->error_code == 0 && n_interp > 0) {
        if(result->error_code == 0) { result->error_code = -10; }
    }
    return result;
}

void free_photon_trajectory(PhotonTrajectory *traj) {
    if (!traj) { return; }
    if(traj->K) { gsl_vector_free(traj->K); }
    if(traj->r) { gsl_vector_free(traj->r); }
    if(traj->phi) { gsl_vector_free(traj->phi); }
    if(traj->x) { gsl_vector_free(traj->x); }
    if(traj->y) { gsl_vector_free(traj->y); }
    if (traj->crossings_y_at_x_targets) {
        for (size_t i = 0; i < traj->num_x_targets; ++i) {
            if(traj->crossings_y_at_x_targets[i]) { gsl_vector_free(traj->crossings_y_at_x_targets[i]); }
        }
        free(traj->crossings_y_at_x_targets);
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
    if (!result) { perror("compute_trajectory_crossings_only: calloc TrajectoryCrossings failed"); return NULL; }
    result->error_code = 0;
    result->crossings_y_at_x_targets = NULL; result->num_x_targets = 0;

    int core_status = integrate_photon_trajectory_core(
        r_0, phi_0, M, psi, r_max, x_stop_val, x_stop_active, x_targets, t_end,
        1e-8, 1e-10, NULL, result, 0 );

    if (core_status != 0 && result->error_code == 0) {
        result->error_code = core_status;
    }
    return result;
}

void free_trajectory_crossings(TrajectoryCrossings *crossings) {
    if (!crossings) { return; }
    if (crossings->crossings_y_at_x_targets) {
        for (size_t i = 0; i < crossings->num_x_targets; ++i) {
            if(crossings->crossings_y_at_x_targets[i]) { gsl_vector_free(crossings->crossings_y_at_x_targets[i]); }
        }
        free(crossings->crossings_y_at_x_targets);
    }
    free(crossings);
}


// --- Image Mapping Functions ---
ImageMapResult image_map(
    double Y, double Z, double x_0_plane, double x_1_window, double x_2_observer, double a_scale) {
    ImageMapResult im_result = { .miss = 1, .y_window = NAN, .y_image = NAN, .error_code = 0 };

    double a_prime = (fabs(x_2_observer - x_1_window) < EPSILON_GENERAL) ? a_scale :
                     ((x_2_observer - x_0_plane) / (x_2_observer - x_1_window) * a_scale);
    double r_max_integration = sqrt(2.0 * a_prime * a_prime + x_0_plane * x_0_plane);
    if (r_max_integration < (2.0 * 1.0 + 100.0)) { r_max_integration = (2.0 * 1.0 + 100.0); }

    double rho_sq = Y * Y + Z * Z;
    double observer_to_window_dist_sq = (x_2_observer - x_1_window) * (x_2_observer - x_1_window) + rho_sq;
    if (observer_to_window_dist_sq <= EPSILON_GENERAL) { im_result.error_code = 1; return im_result; }

    double arg_acos = sqrt(rho_sq) / sqrt(observer_to_window_dist_sq);
    arg_acos = fmax(-1.0, fmin(1.0, arg_acos));
    double psi_calc = -acos(arg_acos);

    gsl_vector *x_targets_map = gsl_vector_alloc(2);
    if (!x_targets_map) { im_result.error_code = -1; return im_result; }
    gsl_vector_set(x_targets_map, 0, x_0_plane);
    gsl_vector_set(x_targets_map, 1, x_1_window);

    double x_stop_map = x_0_plane - 1.0;

    TrajectoryCrossings *crossings = compute_trajectory_crossings_only(
        x_2_observer, 0.0, 1.0, psi_calc, r_max_integration,
        x_stop_map, true, x_targets_map, DEFAULT_T_END);

    gsl_vector_free(x_targets_map);

    if (!crossings) { im_result.error_code = -1; return im_result; }
    if (crossings->error_code != 0) {
        im_result.error_code = crossings->error_code;
        free_trajectory_crossings(crossings);
        return im_result;
    }

    if (crossings->num_x_targets == 2 &&
        crossings->crossings_y_at_x_targets &&
        crossings->crossings_y_at_x_targets[0] && crossings->crossings_y_at_x_targets[0]->size > 0 &&
        crossings->crossings_y_at_x_targets[1] && crossings->crossings_y_at_x_targets[1]->size > 0) {
        im_result.y_image = gsl_vector_get(crossings->crossings_y_at_x_targets[0], 0);
        im_result.y_window = gsl_vector_get(crossings->crossings_y_at_x_targets[1], 0);
        im_result.miss = 0;
    } else {
        im_result.miss = 1;
    }
    free_trajectory_crossings(crossings);
    return im_result;
}

ImageMapResult image_map_radial(
    double r_rho, double x_0_plane, double x_1_window, double x_2_observer, double a_scale_factor) {
    return image_map(r_rho, 0.0, x_0_plane, x_1_window, x_2_observer, a_scale_factor);
}


// --- Results Generation Functions ---
char** results_cartesian(
    double x_0, double x_1, double x_2, double a_param, int n_param,
    size_t chunk_size_photons, int* out_num_files_created) {

    if (out_num_files_created) { *out_num_files_created = 0; }
    PhotonMapDataPoint *results_chunk_buffer = malloc(chunk_size_photons * sizeof(PhotonMapDataPoint));
    if (!results_chunk_buffer) {
        fprintf(stderr, "results_cartesian: Failed to allocate chunk buffer.\n");
        return NULL;
    }
    size_t current_chunk_fill = 0;
    char **saved_files_list = NULL;
    int saved_files_count = 0; 
    int saved_files_capacity = 0;
    int chunk_index = 0;

    printf("Starting Results_Cartesian generation (x0=%.2f, x1=%.2f, x2=%.2f, a=%.2f, n=%d)...\n", x_0, x_1, x_2, a_param, n_param);
    time_t start_time_sec = time(NULL);

    char x0_str[64], x1_str[64], x2_str[64], n_str_fn[64];
    format_coord_str_for_filename(x0_str, sizeof(x0_str), x_0);
    format_coord_str_for_filename(x1_str, sizeof(x1_str), x_1);
    format_coord_str_for_filename(x2_str, sizeof(x2_str), x_2);
    snprintf(n_str_fn, sizeof(n_str_fn), "%d", n_param);
    char base_save_path[MAX_FILENAME_LEN];
    snprintf(base_save_path, MAX_FILENAME_LEN, "Results_Cartesian_%s_%s_%s_%s", x0_str, x1_str, x2_str, n_str_fn);

    for (int i = 0; i <= n_param; ++i) {
        double y_sample = (n_param > 0) ? (a_param * i / n_param) : ( (i > 0) ? a_param : 0.0);
        for (int j = 0; j <= i; ++j) {
            double z_sample = (n_param > 0) ? (a_param * j / n_param) : ( (j > 0) ? a_param : 0.0);

            if (j % 10000 == 0 && (i % (n_param/10 + 1) == 0 || i == n_param) ) {
                printf("  y=%.3f, z=%.3f (Chunk %d, Fill %zu/%zu)\n", y_sample, z_sample, chunk_index, current_chunk_fill, chunk_size_photons);
            }

            double phi_prime_base = (fabs(y_sample) < EPSILON_GENERAL && fabs(z_sample) < EPSILON_GENERAL) ? 0.0 : atan2(z_sample, y_sample);
            ImageMapResult map_output = image_map(y_sample, z_sample, x_0, x_1, x_2, a_param);

            if (map_output.miss || map_output.error_code != 0) { continue; }

            double y_w = map_output.y_window;
            double y_i = map_output.y_image;
            double phi_values_to_use[8];
            int num_phi_sym_points;

            if (fabs(y_sample - z_sample) < EPSILON_GENERAL || fabs(z_sample) < EPSILON_GENERAL) {
                num_phi_sym_points = 4;
                phi_values_to_use[0] = phi_prime_base;
                phi_values_to_use[1] = phi_prime_base + M_PI / 2.0;
                phi_values_to_use[2] = phi_prime_base + M_PI;
                phi_values_to_use[3] = phi_prime_base + 3.0 * M_PI / 2.0;
            } else {
                num_phi_sym_points = 8;
                phi_values_to_use[0] = phi_prime_base;
                phi_values_to_use[1] = phi_prime_base + M_PI / 2.0;
                phi_values_to_use[2] = phi_prime_base + M_PI;
                phi_values_to_use[3] = phi_prime_base + 3.0 * M_PI / 2.0;
                phi_values_to_use[4] = M_PI / 2.0 - phi_prime_base;
                phi_values_to_use[5] = M_PI - phi_prime_base;
                phi_values_to_use[6] = 3.0 * M_PI / 2.0 - phi_prime_base;
                phi_values_to_use[7] = 2.0 * M_PI - phi_prime_base;
            }

            for (int k_phi = 0; k_phi < num_phi_sym_points; ++k_phi) {
                double phi_p_val_norm = normalize_phi(phi_values_to_use[k_phi]);
                results_chunk_buffer[current_chunk_fill].y_window_cart_x = y_w * cos(phi_p_val_norm);
                results_chunk_buffer[current_chunk_fill].y_window_cart_y = y_w * sin(phi_p_val_norm);
                results_chunk_buffer[current_chunk_fill].y_image_cart_x  = y_i * cos(phi_p_val_norm);
                results_chunk_buffer[current_chunk_fill].y_image_cart_y  = y_i * sin(phi_p_val_norm);
                current_chunk_fill++;

                if (current_chunk_fill >= chunk_size_photons) {
                    char chunk_filename[MAX_FILENAME_LEN];
                    snprintf(chunk_filename, MAX_FILENAME_LEN, "%s_chunk_%03d.bin", base_save_path, chunk_index);
                    if (save_photon_data_chunk(chunk_filename, results_chunk_buffer, current_chunk_fill) == 0) {
                        if (add_string_to_list(&saved_files_list, &saved_files_count, &saved_files_capacity, chunk_filename) != 0) {
                            fprintf(stderr, "results_cartesian: Failed to add filename to list.\n");
                            free(results_chunk_buffer); free_string_array(saved_files_list, saved_files_count); return NULL;
                        }
                    } else { fprintf(stderr, "results_cartesian: Failed to save chunk %d.\n", chunk_index); }
                    current_chunk_fill = 0;
                    chunk_index++;
                }
            }
        }
    }

    if (current_chunk_fill > 0) {
        char chunk_filename[MAX_FILENAME_LEN];
        snprintf(chunk_filename, MAX_FILENAME_LEN, "%s_chunk_%03d.bin", base_save_path, chunk_index);
        if (save_photon_data_chunk(chunk_filename, results_chunk_buffer, current_chunk_fill) == 0) {
            if (add_string_to_list(&saved_files_list, &saved_files_count, &saved_files_capacity, chunk_filename) != 0) {
                 fprintf(stderr, "results_cartesian: Failed to add final filename to list.\n");
                 free(results_chunk_buffer); free_string_array(saved_files_list, saved_files_count); return NULL;
            }
        } else { fprintf(stderr, "results_cartesian: Failed to save final chunk %d.\n", chunk_index); }
    }

    free(results_chunk_buffer);
    if (out_num_files_created) { *out_num_files_created = saved_files_count; }
    add_string_to_list(&saved_files_list, &saved_files_count, &saved_files_capacity, NULL); 

    printf("\nFinished Results_Cartesian. Time: %ld s. Files created: %d\n", time(NULL) - start_time_sec, *out_num_files_created);
    return saved_files_list;
}

char** results_radial(
    double x_0, double x_1, double x_2, double R_max_sample, int n_param,
    size_t chunk_size_photons, int* out_num_files_created) {

    if (out_num_files_created) { *out_num_files_created = 0; }
    PhotonMapDataPoint *results_chunk_buffer = malloc(chunk_size_photons * sizeof(PhotonMapDataPoint));
    if (!results_chunk_buffer) {
        fprintf(stderr, "results_radial: Failed to allocate chunk buffer.\n");
        return NULL;
    }
    size_t current_chunk_fill = 0;
    char **saved_files_list = NULL;
    int saved_files_count = 0;
    int saved_files_capacity = 0;
    int chunk_index = 0;

    printf("Starting Results_Radial generation (x0=%.2f, x1=%.2f, x2=%.2f, Rmax=%.2f, n=%d)...\n", x_0, x_1, x_2, R_max_sample, n_param);
    time_t start_time_sec = time(NULL);

    double a_scale_for_map = 1.0 * (x_2 - x_1);
    if (fabs(a_scale_for_map) < EPSILON_GENERAL) { a_scale_for_map = 1.0; }

    char x0_str[64], x1_str[64], x2_str[64], n_str_fn[64];
    format_coord_str_for_filename(x0_str, sizeof(x0_str), x_0);
    format_coord_str_for_filename(x1_str, sizeof(x1_str), x_1);
    format_coord_str_for_filename(x2_str, sizeof(x2_str), x_2);
    snprintf(n_str_fn, sizeof(n_str_fn), "%d", n_param);
    char base_save_path[MAX_FILENAME_LEN];
    snprintf(base_save_path, MAX_FILENAME_LEN, "Results_rainbow_thesis_%s_%s_%s_%s", x0_str, x1_str, x2_str, n_str_fn);

    for (int i = 0; i <= n_param; ++i) {
        double r_s_val = (n_param > 0) ? (R_max_sample * i / n_param) : ( (i > 0) ? R_max_sample : 0.0);
        int k_angular_steps = (int)fmax(1.0, (5.0 * 1e4 * r_s_val));

        if (i % (n_param/10 +1) == 0 || i == n_param) {
             printf("  r_sample=%.4f, k_angular_steps=%d (Chunk %d)\n", r_s_val, k_angular_steps, chunk_index);
        }

        ImageMapResult map_output = image_map_radial(r_s_val, x_0, x_1, x_2, a_scale_for_map);
        if (map_output.miss || map_output.error_code != 0) { continue; }

        double y_w = map_output.y_window;
        double y_i = map_output.y_image;

        for (int j = 0; j <= k_angular_steps; ++j) {
            double phi_prime_val = (k_angular_steps > 0) ? ((double)j / k_angular_steps * 2.0 * M_PI) : 0.0;
            double cos_phi = cos(phi_prime_val);
            double sin_phi = sin(phi_prime_val);

            results_chunk_buffer[current_chunk_fill].y_window_cart_x = y_w * cos_phi;
            results_chunk_buffer[current_chunk_fill].y_window_cart_y = y_w * sin_phi;
            results_chunk_buffer[current_chunk_fill].y_image_cart_x  = y_i * cos_phi;
            results_chunk_buffer[current_chunk_fill].y_image_cart_y  = y_i * sin_phi;
            current_chunk_fill++;

            if (current_chunk_fill >= chunk_size_photons) {
                char chunk_filename[MAX_FILENAME_LEN];
                snprintf(chunk_filename, MAX_FILENAME_LEN, "%s_chunk_%03d.bin", base_save_path, chunk_index);
                if (save_photon_data_chunk(chunk_filename, results_chunk_buffer, current_chunk_fill) == 0) {
                    if (add_string_to_list(&saved_files_list, &saved_files_count, &saved_files_capacity, chunk_filename) != 0) {
                         fprintf(stderr, "results_radial: Failed to add filename to list.\n");
                         free(results_chunk_buffer); free_string_array(saved_files_list, saved_files_count); return NULL;
                    }
                } else { fprintf(stderr, "results_radial: Failed to save chunk %d.\n", chunk_index); }
                current_chunk_fill = 0;
                chunk_index++;
            }
        }
    }

    if (current_chunk_fill > 0) {
        char chunk_filename[MAX_FILENAME_LEN];
        snprintf(chunk_filename, MAX_FILENAME_LEN, "%s_chunk_%03d.bin", base_save_path, chunk_index);
        if (save_photon_data_chunk(chunk_filename, results_chunk_buffer, current_chunk_fill) == 0) {
            if (add_string_to_list(&saved_files_list, &saved_files_count, &saved_files_capacity, chunk_filename) != 0) {
                 fprintf(stderr, "results_radial: Failed to add final filename to list.\n");
                 free(results_chunk_buffer); free_string_array(saved_files_list, saved_files_count); return NULL;
            }
        } else { fprintf(stderr, "results_radial: Failed to save final chunk %d.\n", chunk_index); }
    }

    free(results_chunk_buffer);
    if (out_num_files_created) { *out_num_files_created = saved_files_count; }
    add_string_to_list(&saved_files_list, &saved_files_count, &saved_files_capacity, NULL);
    printf("\nFinished Results_Radial. Time: %ld s. Files created: %d\n", time(NULL) - start_time_sec, *out_num_files_created);
    return saved_files_list;
}


ResultsRadialLightRingOutput* results_radial_light_ring(
    double x_0, double x_1, double x_2, int n_param,
    size_t chunk_size_photons) {

    ResultsRadialLightRingOutput *output = calloc(1, sizeof(ResultsRadialLightRingOutput));
    if (!output) { perror("results_radial_light_ring: calloc output failed"); return NULL; }
    output->r_sample_misses = gsl_vector_alloc(INITIAL_CROSSING_POINTS_CAPACITY);
    output->r_sample_hits = gsl_vector_alloc(INITIAL_CROSSING_POINTS_CAPACITY);
    output->saved_files_list = NULL;
    output->num_files_created = 0;
    output->capacity_saved_files = 0;

    if (!output->r_sample_misses || !output->r_sample_hits) {
        free_results_radial_light_ring_output(output); return NULL;
    }
    gsl_vector_set_zero(output->r_sample_misses); output->r_sample_misses->size = 0;
    gsl_vector_set_zero(output->r_sample_hits); output->r_sample_hits->size = 0;
    size_t misses_capacity = INITIAL_CROSSING_POINTS_CAPACITY;
    size_t hits_capacity = INITIAL_CROSSING_POINTS_CAPACITY;

    PhotonMapDataPoint *results_chunk_buffer = malloc(chunk_size_photons * sizeof(PhotonMapDataPoint));
    if (!results_chunk_buffer) {
        fprintf(stderr, "results_radial_light_ring: Failed to allocate chunk buffer.\n");
        free_results_radial_light_ring_output(output); return NULL;
    }
    size_t current_chunk_fill = 0;
    int chunk_index = 0;

    printf("Starting Results_Radial_Light_Ring (x0=%.2f, x1=%.2f, x2=%.2f, n=%d)...\n", x_0, x_1, x_2, n_param);
    time_t start_time_sec = time(NULL);

    double a_scale_lr = 1.0 * (x_2 - x_1);
    if (fabs(a_scale_lr) < EPSILON_GENERAL) { a_scale_lr = 1.0; }

    char x0_str[64], x1_str[64], x2_str[64], n_str_fn[64];
    format_coord_str_for_filename(x0_str, sizeof(x0_str), x_0);
    format_coord_str_for_filename(x1_str, sizeof(x1_str), x_1);
    format_coord_str_for_filename(x2_str, sizeof(x2_str), x_2);
    snprintf(n_str_fn, sizeof(n_str_fn), "%d", n_param);
    char base_save_path[MAX_FILENAME_LEN];
    snprintf(base_save_path, MAX_FILENAME_LEN, "Light_Ring_Segement_0th_large_%s_%s_%s_%s", x0_str, x1_str, x2_str, n_str_fn);

    const double r_start_lr_const = 0.2543706590950274;
    const double r_delta_lr_const = 0.0008610946154524735;

    for (int i = 0; i <= n_param; ++i) {
        double r_s_val = r_start_lr_const + ((n_param > 0) ? ((double)i / n_param * r_delta_lr_const) : ( (i > 0) ? r_delta_lr_const : 0.0));
        int k_angular_steps_lr = (int)fmax(1.0, (6.0 * 1e4 * r_s_val));

        if (i % (n_param/10 +1) == 0 || i == n_param) {
             printf("  r_sample=%.7f, k_angular_steps=%d (Chunk %d)\n", r_s_val, k_angular_steps_lr, chunk_index);
        }

        ImageMapResult map_output = image_map_radial(r_s_val, x_0, x_1, x_2, a_scale_lr);
        if (map_output.miss || map_output.error_code != 0) {
            gsl_vector_dynamic_append(&output->r_sample_misses, r_s_val, &misses_capacity);
            continue;
        }
        gsl_vector_dynamic_append(&output->r_sample_hits, r_s_val, &hits_capacity);

        double y_w = map_output.y_window;
        double y_i = map_output.y_image;
        const double phi_prime_start_angle_lr = 3.0/2.0 * M_PI;
        const double phi_prime_angular_width_lr = M_PI / 180.0 / 10.0;

        for (int j = 0; j <= k_angular_steps_lr; ++j) {
            double phi_prime_val = phi_prime_start_angle_lr + ((k_angular_steps_lr > 0) ? ((double)j / k_angular_steps_lr * phi_prime_angular_width_lr) : 0.0);
            double cos_phi = cos(phi_prime_val);
            double sin_phi = sin(phi_prime_val);

            results_chunk_buffer[current_chunk_fill].y_window_cart_x = y_w * cos_phi;
            results_chunk_buffer[current_chunk_fill].y_window_cart_y = y_w * sin_phi;
            results_chunk_buffer[current_chunk_fill].y_image_cart_x  = y_i * cos_phi;
            results_chunk_buffer[current_chunk_fill].y_image_cart_y  = y_i * sin_phi;
            current_chunk_fill++;

            if (current_chunk_fill >= chunk_size_photons) {
                char chunk_filename[MAX_FILENAME_LEN];
                snprintf(chunk_filename, MAX_FILENAME_LEN, "%s_chunk_%03d.bin", base_save_path, chunk_index);
                if (save_photon_data_chunk(chunk_filename, results_chunk_buffer, current_chunk_fill) == 0) {
                    if (add_string_to_list(&output->saved_files_list, &output->num_files_created, &output->capacity_saved_files, chunk_filename) != 0) {
                         fprintf(stderr, "results_radial_light_ring: Failed to add filename to list.\n");
                         free(results_chunk_buffer); free_results_radial_light_ring_output(output); return NULL;
                    }
                } else { fprintf(stderr, "results_radial_light_ring: Failed to save chunk %d.\n", chunk_index); }
                current_chunk_fill = 0;
                chunk_index++;
            }
        }
    }

    if (current_chunk_fill > 0) {
        char chunk_filename[MAX_FILENAME_LEN];
        snprintf(chunk_filename, MAX_FILENAME_LEN, "%s_chunk_%03d.bin", base_save_path, chunk_index);
        if (save_photon_data_chunk(chunk_filename, results_chunk_buffer, current_chunk_fill) == 0) {
            if (add_string_to_list(&output->saved_files_list, &output->num_files_created, &output->capacity_saved_files, chunk_filename) != 0) {
                 fprintf(stderr, "results_radial_light_ring: Failed to add final filename to list.\n");
                 free(results_chunk_buffer); free_results_radial_light_ring_output(output); return NULL;
            }
        } else { fprintf(stderr, "results_radial_light_ring: Failed to save final chunk %d.\n", chunk_index); }
    }

    free(results_chunk_buffer);
    add_string_to_list(&output->saved_files_list, &output->num_files_created, &output->capacity_saved_files, NULL);
    printf("\nFinished Results_Radial_Light_Ring. Time: %ld s. Files created: %d\n", time(NULL) - start_time_sec, output->num_files_created);
    return output;
}

void free_results_radial_light_ring_output(ResultsRadialLightRingOutput* output) {
    if (!output) { return; }
    if (output->saved_files_list) {
        free_string_array(output->saved_files_list, output->num_files_created); 
        output->saved_files_list = NULL;
    }
    if(output->r_sample_misses) { gsl_vector_free(output->r_sample_misses); output->r_sample_misses = NULL;}
    if(output->r_sample_hits) { gsl_vector_free(output->r_sample_hits); output->r_sample_hits = NULL; }
    free(output);
}

void free_string_array(char** arr, int count) {
    if (!arr) { return; }
    for (int i = 0; i < count; ++i) { 
        if (arr[i]) { free(arr[i]); }
    }
    free(arr);
}


// --- Image Rendering from Photon Data ---
int map_photons(
    const char *image_source_path, PPMImage *image_source_data_param,
    char **photon_chunk_files, int num_chunk_files_param,
    const char *save_path,
    const double *dest_logical_bounds_param, double pixels_per_logical_unit_param,
    const int *output_shape_param, RGBColor default_color_param, bool flip_z_axis_render) {

    int mapped_h = 1, mapped_w = 1;
    bool using_manual_dest_bounds = false;
    double dest_map_bound_min_y = 0.0, dest_map_bound_max_y = 0.0;
    double dest_map_bound_min_z = 0.0, dest_map_bound_max_z = 0.0;
    PPMImage *source_image_loaded = NULL;
    const PPMImage *source_image_to_use = NULL;
    unsigned char *mapped_image_pixels = NULL;
    double *sum_r_array = NULL, *sum_g_array = NULL, *sum_b_array = NULL;
    int64_t *count_array_pixels = NULL;
    int final_status = -1;

    if (dest_logical_bounds_param != NULL) {
        using_manual_dest_bounds = true;
        if (pixels_per_logical_unit_param <= 0) {
            fprintf(stderr, "map_photons: pixels_per_logical_unit must be positive for manual bounds mode.\n");
            goto cleanup_map_photons;
        }
        dest_map_bound_min_y = dest_logical_bounds_param[0]; dest_map_bound_max_y = dest_logical_bounds_param[1];
        dest_map_bound_min_z = dest_logical_bounds_param[2]; dest_map_bound_max_z = dest_logical_bounds_param[3];
        if (dest_map_bound_min_y >= dest_map_bound_max_y || dest_map_bound_min_z >= dest_map_bound_max_z) {
            fprintf(stderr, "map_photons: Invalid dest_logical_bounds (min >= max).\n"); goto cleanup_map_photons;
        }
        mapped_w = (int)ceil((dest_map_bound_max_y - dest_map_bound_min_y) * pixels_per_logical_unit_param);
        mapped_h = (int)ceil((dest_map_bound_max_z - dest_map_bound_min_z) * pixels_per_logical_unit_param);
    } else if (output_shape_param != NULL) {
        using_manual_dest_bounds = false;
        mapped_h = output_shape_param[0]; mapped_w = output_shape_param[1];
        if (mapped_h <= 0 || mapped_w <= 0) {
            fprintf(stderr, "map_photons: output_shape dimensions must be positive.\n"); goto cleanup_map_photons;
        }
    } else {
        fprintf(stderr, "map_photons: Must provide either dest_logical_bounds or output_shape.\n"); goto cleanup_map_photons;
    }
    if (mapped_w <= 0) { mapped_w = 1; }
    if (mapped_h <= 0) { mapped_h = 1; }


    if (image_source_data_param) {
        source_image_to_use = image_source_data_param;
    } else if (image_source_path) {
        source_image_loaded = load_ppm_image(image_source_path);
        if (!source_image_loaded) {
            fprintf(stderr, "map_photons: Failed to load source image from '%s'.\n", image_source_path);
            goto cleanup_map_photons;
        }
        source_image_to_use = source_image_loaded;
    } else {
        fprintf(stderr, "map_photons: No source image provided (path or data).\n"); goto cleanup_map_photons;
    }
    if (!source_image_to_use || source_image_to_use->width <= 0 || source_image_to_use->height <= 0 || source_image_to_use->channels != 3) {
        fprintf(stderr, "map_photons: Invalid source image properties.\n"); goto cleanup_map_photons;
    }
    int orig_h = source_image_to_use->height; int orig_w = source_image_to_use->width;
    double source_pixel_aspect = (double)orig_w / orig_h;

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

    int num_chunks_to_process = 0;
    if (photon_chunk_files) {
        if (num_chunk_files_param >= 0) {
            num_chunks_to_process = num_chunk_files_param;
        } else {
            for (num_chunks_to_process = 0; photon_chunk_files[num_chunks_to_process] != NULL; ++num_chunks_to_process);
        }
    }

    if (num_chunks_to_process == 0) {
        printf("map_photons: No photon chunk files provided. Saving default background image.\n");
        PPMImage temp_out_img = {mapped_image_pixels, mapped_w, mapped_h, 3};
        if (save_ppm_image(save_path, &temp_out_img) == 0) { final_status = 0; }
        goto cleanup_map_photons;
    }

    double global_max_abs_y0_scan = EPSILON_GENERAL, global_max_abs_z0_scan = EPSILON_GENERAL;
    double global_max_abs_y1_scan = EPSILON_GENERAL, global_max_abs_z1_scan = EPSILON_GENERAL;
    bool found_any_valid_photons = false;

    for (int ci = 0; ci < num_chunks_to_process; ++ci) {
        size_t num_pts_in_chunk;
        PhotonMapDataPoint* chunk_data = load_photon_data_chunk(photon_chunk_files[ci], &num_pts_in_chunk);
        if (chunk_data && num_pts_in_chunk > 0) {
            found_any_valid_photons = true;
            for (size_t i = 0; i < num_pts_in_chunk; ++i) {
                if (gsl_isnan(chunk_data[i].y_window_cart_x) || gsl_isnan(chunk_data[i].y_window_cart_y) ||
                    gsl_isnan(chunk_data[i].y_image_cart_x)  || gsl_isnan(chunk_data[i].y_image_cart_y)) { continue; }
                if (!using_manual_dest_bounds) {
                    if (fabs(chunk_data[i].y_window_cart_x) > global_max_abs_y0_scan) { global_max_abs_y0_scan = fabs(chunk_data[i].y_window_cart_x); }
                    if (fabs(chunk_data[i].y_window_cart_y) > global_max_abs_z0_scan) { global_max_abs_z0_scan = fabs(chunk_data[i].y_window_cart_y); }
                }
                if (fabs(chunk_data[i].y_image_cart_x) > global_max_abs_y1_scan) { global_max_abs_y1_scan = fabs(chunk_data[i].y_image_cart_x); }
                if (fabs(chunk_data[i].y_image_cart_y) > global_max_abs_z1_scan) { global_max_abs_z1_scan = fabs(chunk_data[i].y_image_cart_y); }
            }
        }
        free(chunk_data);
    }

    if (!found_any_valid_photons) {
        fprintf(stderr, "map_photons: No valid (non-NaN) photons found in any chunk during scan.\n");
        PPMImage temp_out_img = {mapped_image_pixels, mapped_w, mapped_h, 3};
        if (save_ppm_image(save_path, &temp_out_img) == 0) { final_status = 0; }
        goto cleanup_map_photons;
    }

    double source_map_bound_y1_extent = global_max_abs_y1_scan;
    double source_map_bound_z1_extent = global_max_abs_z1_scan;
    if (source_pixel_aspect >= 1.0) {
         source_map_bound_z1_extent = fmax(source_map_bound_z1_extent, source_map_bound_y1_extent / source_pixel_aspect);
    } else {
         source_map_bound_y1_extent = fmax(source_map_bound_y1_extent, source_map_bound_z1_extent * source_pixel_aspect);
    }
    double source_denom_y1 = fmax(2.0 * source_map_bound_y1_extent, EPSILON_GENERAL);
    double source_denom_z1 = fmax(2.0 * source_map_bound_z1_extent, EPSILON_GENERAL);

    double dest_denom_y0, dest_denom_z0;
    if (using_manual_dest_bounds) {
        dest_denom_y0 = fmax(dest_map_bound_max_y - dest_map_bound_min_y, EPSILON_GENERAL);
        dest_denom_z0 = fmax(dest_map_bound_max_z - dest_map_bound_min_z, EPSILON_GENERAL);
    } else {
        double output_pixel_aspect = (double)mapped_w / mapped_h;
        double logical_dest_y0_extent = global_max_abs_y0_scan;
        double logical_dest_z0_extent = global_max_abs_z0_scan;
        if (output_pixel_aspect >= 1.0) {
            logical_dest_z0_extent = fmax(logical_dest_z0_extent, logical_dest_y0_extent / output_pixel_aspect);
        } else {
            logical_dest_y0_extent = fmax(logical_dest_y0_extent, logical_dest_z0_extent * output_pixel_aspect);
        }
        dest_map_bound_min_y = -logical_dest_y0_extent; dest_map_bound_max_y = logical_dest_y0_extent;
        dest_map_bound_min_z = -logical_dest_z0_extent; dest_map_bound_max_z = logical_dest_z0_extent;
        dest_denom_y0 = fmax(2.0 * logical_dest_y0_extent, EPSILON_GENERAL);
        dest_denom_z0 = fmax(2.0 * logical_dest_z0_extent, EPSILON_GENERAL);
    }

    for (int ci = 0; ci < num_chunks_to_process; ++ci) {
        printf("\rProcessing chunk %d/%d: %s...", ci + 1, num_chunks_to_process, photon_chunk_files[ci]); fflush(stdout);
        size_t num_pts_in_chunk;
        PhotonMapDataPoint* chunk_data = load_photon_data_chunk(photon_chunk_files[ci], &num_pts_in_chunk);
        if (!chunk_data || num_pts_in_chunk == 0) { free(chunk_data); continue; }

        for (size_t i = 0; i < num_pts_in_chunk; ++i) {
            double y0_f = chunk_data[i].y_window_cart_x; double z0_f = chunk_data[i].y_window_cart_y;
            double y1_f = chunk_data[i].y_image_cart_x;  double z1_f = chunk_data[i].y_image_cart_y;

            if (gsl_isnan(y0_f) || gsl_isnan(z0_f) || gsl_isnan(y1_f) || gsl_isnan(z1_f)) { continue; }
            if (fabs(y1_f) > source_map_bound_y1_extent + EPSILON_GENERAL || fabs(z1_f) > source_map_bound_z1_extent + EPSILON_GENERAL) { continue; }
            if (y0_f < dest_map_bound_min_y - EPSILON_GENERAL || y0_f > dest_map_bound_max_y + EPSILON_GENERAL ||
                z0_f < dest_map_bound_min_z - EPSILON_GENERAL || z0_f > dest_map_bound_max_z + EPSILON_GENERAL) { continue; }

            int src_col_idx = (int)round(((y1_f + source_map_bound_y1_extent) / source_denom_y1) * (orig_w - 1));
            int src_row_idx = (int)round(((z1_f + source_map_bound_z1_extent) / source_denom_z1) * (orig_h - 1));
            src_col_idx = (src_col_idx < 0) ? 0 : (src_col_idx >= orig_w ? orig_w - 1 : src_col_idx);
            src_row_idx = (src_row_idx < 0) ? 0 : (src_row_idx >= orig_h ? orig_h - 1 : src_row_idx);

            int dest_col_idx = (int)round(((y0_f - dest_map_bound_min_y) / dest_denom_y0) * (mapped_w - 1));
            int dest_row_idx = (flip_z_axis_render) ?
                               (int)round(((dest_map_bound_max_z - z0_f) / dest_denom_z0) * (mapped_h - 1)) :
                               (int)round(((z0_f - dest_map_bound_min_z) / dest_denom_z0) * (mapped_h - 1));
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
    printf("\n");

    for (size_t i = 0; i < num_output_pixels; ++i) {
        if (count_array_pixels[i] > 0) {
            mapped_image_pixels[i*3 + 0] = (unsigned char)fmax(0, fmin(255, round(sum_r_array[i] / count_array_pixels[i])));
            mapped_image_pixels[i*3 + 1] = (unsigned char)fmax(0, fmin(255, round(sum_g_array[i] / count_array_pixels[i])));
            mapped_image_pixels[i*3 + 2] = (unsigned char)fmax(0, fmin(255, round(sum_b_array[i] / count_array_pixels[i])));
        }
    }

    PPMImage final_img = {mapped_image_pixels, mapped_w, mapped_h, 3};
    if (save_ppm_image(save_path, &final_img) == 0) {
        final_status = 0;
    } else {
        fprintf(stderr, "map_photons: Failed to save final image to '%s'.\n", save_path);
    }

cleanup_map_photons:
    if(source_image_loaded) { free_ppm_image(source_image_loaded); }
    if(mapped_image_pixels) { free(mapped_image_pixels); }
    if(sum_r_array) { free(sum_r_array); }
    if(sum_g_array) { free(sum_g_array); }
    if(sum_b_array) { free(sum_b_array); }
    if(count_array_pixels) { free(count_array_pixels); }
    return final_status;
}
// schwarzschild_tracer.h
#ifndef SCHWARZSCHILD_TRACER_H
#define SCHWARZSCHILD_TRACER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h> // For bool type

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sys.h> // For gsl_isnan etc.
#include <float.h>    // For DBL_MAX, DBL_EPSILON
#include <ctype.h>    // For isspace in PPM loader

// --- Constants ---
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define DEFAULT_T_END 1e6
#define DEFAULT_NUM_INTERP_POINTS 1000
#define DEFAULT_CHUNK_SIZE_PHOTONS 1000000 // Number of PhotonMapDataPoint elements
#define EPSILON_GENERAL 1e-9         // General small epsilon for comparisons
#define EVENT_DETECTION_TOLERANCE 1e-9 // Tolerance for GSL root finder in event detection
#define MAX_FILENAME_LEN 2048
#define INITIAL_RAW_POINTS_CAPACITY 256 // For collecting raw trajectory points within a segment
#define INITIAL_SEGMENTS_CAPACITY 10    // For collecting trajectory segments
#define INITIAL_CROSSING_POINTS_CAPACITY 8 // For collecting crossing points per target
#define INITIAL_SAVED_FILES_CAPACITY 10 // For lists of saved filenames

// --- Struct Definitions ---

typedef struct {
    double M;
    double b;
    int *sign_dr_dk;
} ODEParams;

typedef struct {
    double M;
    double b;
    double r_max_event;
    double x_stop_event;
    double x_target_event;
    ODEParams *ode_p_event;
} EventFunctionParams;

typedef struct {
    gsl_vector *x;
    gsl_vector *y;
    gsl_vector *r;
    gsl_vector *phi;
    gsl_vector *K;
    gsl_vector **crossings_y_at_x_targets;
    size_t num_x_targets;
    int error_code;
} PhotonTrajectory;

typedef struct {
    gsl_vector **crossings_y_at_x_targets;
    size_t num_x_targets;
    int error_code;
} TrajectoryCrossings;

typedef struct {
    int miss;
    double y_window;
    double y_image;
    int error_code;
} ImageMapResult;

typedef struct {
    double y_window_cart_x;
    double y_window_cart_y;
    double y_image_cart_x;
    double y_image_cart_y;
} PhotonMapDataPoint;

typedef struct {
    unsigned char r, g, b;
} RGBColor;

typedef struct {
    unsigned char *data;
    int width;
    int height;
    int channels;
} PPMImage;

typedef struct {
    char** saved_files_list;
    gsl_vector* r_sample_misses;
    gsl_vector* r_sample_hits;
    int num_files_created; // Actual count of file strings
    int capacity_saved_files;
} ResultsRadialLightRingOutput;


// --- Function Declarations ---

PhotonTrajectory* compute_trajectory(
    double r_0, double phi_0, double M, double psi, double r_max,
    double x_stop_val, bool x_stop_active,
    const gsl_vector *x_targets,
    double t_end, int num_interp_points
);
void free_photon_trajectory(PhotonTrajectory *traj);

TrajectoryCrossings* compute_trajectory_crossings_only(
    double r_0, double phi_0, double M, double psi, double r_max,
    double x_stop_val, bool x_stop_active,
    const gsl_vector *x_targets,
    double t_end
);
void free_trajectory_crossings(TrajectoryCrossings *crossings);

ImageMapResult image_map(
    double Y, double Z, double x_0_plane, double x_1_window, double x_2_observer, double a_scale
);

ImageMapResult image_map_radial(
    double r_rho, double x_0_plane, double x_1_window, double x_2_observer, double a_scale_factor
);

char** results_cartesian(
    double x_0, double x_1, double x_2, double a_param, int n_param,
    size_t chunk_size_photons, int* out_num_files_created
);

char** results_radial(
    double x_0, double x_1, double x_2, double R_max_sample, int n_param,
    size_t chunk_size_photons, int* out_num_files_created
);

ResultsRadialLightRingOutput* results_radial_light_ring(
    double x_0, double x_1, double x_2, int n_param,
    size_t chunk_size_photons
);
void free_results_radial_light_ring_output(ResultsRadialLightRingOutput* output);
void free_string_array(char** arr, int count);

int map_photons(
    const char *image_source_path,
    PPMImage *image_source_data,
    char **photon_chunk_files,
    int num_chunk_files,
    const char *save_path,
    const double *dest_logical_bounds,
    double pixels_per_logical_unit,
    const int *output_shape,
    RGBColor default_color,
    bool flip_z_axis_render
);

int save_photon_data_chunk(const char *filename, const PhotonMapDataPoint* data, size_t num_points);
PhotonMapDataPoint* load_photon_data_chunk(const char *filename, size_t *num_points_read);

PPMImage* load_ppm_image(const char *filename);
int save_ppm_image(const char *filename, const PPMImage *image);
void free_ppm_image(PPMImage *image);

#endif // SCHWARZSCHILD_TRACER_H
// schwarzschild_tracer.h
#ifndef SCHWARZSCHILD_TRACER_H
#define SCHWARZSCHILD_TRACER_H

// Standard C Headers (that are generally useful for this library)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h> // For bool type
#include <ctype.h>   // For isspace etc. (used in PPM loader, good to have here)
#include <float.h>   // For DBL_MAX, DBL_EPSILON

// GSL Headers
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h> // Kept if any public interface uses it, otherwise can be removed
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sys.h>    // For gsl_isnan etc.

// --- Constants ---
#ifndef M_PI // Guard for M_PI, as math.h might define it with _USE_MATH_DEFINES on some systems
#define M_PI 3.14159265358979323846
#endif
#define MAX_EVENT_TYPES_CORE 100       
#define MAX_K_PER_GSL_SUPER_STEP 0.01
#define DEFAULT_T_END 1e6
#define DEFAULT_NUM_INTERP_POINTS 1000
#define DEFAULT_CHUNK_SIZE_PHOTONS 1000000
#define EPSILON_GENERAL 1e-9
#define EVENT_DETECTION_TOLERANCE 1e-9
#define MAX_FILENAME_LEN 4096
#define INITIAL_RAW_POINTS_CAPACITY 256
#define INITIAL_SEGMENTS_CAPACITY 10
#define INITIAL_CROSSING_POINTS_CAPACITY 8
#define INITIAL_SAVED_FILES_CAPACITY 10


// --- Public Struct Definitions ---

// For ODE solver in trajectory calculations
typedef struct {
    double M;
    double b;
    int *sign_dr_dk; // Pointer to a mutable sign
} ODEParams;

// For event functions in trajectory calculations
typedef struct {
    double M; // Can be different from ODEParams M if events are general
    double b; // Can be different from ODEParams b
    double r_max_event;
    double x_stop_event;
    double x_target_event;
    ODEParams *ode_p_event; // Pointer to the ODEParams for context if needed
} EventFunctionParams;

// For full photon trajectory output
typedef struct {
    gsl_vector *x;
    gsl_vector *y;
    gsl_vector *r;
    gsl_vector *phi;
    gsl_vector *K; // Affine parameter
    gsl_vector **crossings_y_at_x_targets; // Array of GSL vectors for y-crossings
    size_t num_x_targets;
    int error_code; // 0 for success, negative for errors
} PhotonTrajectory;

// For trajectory crossings only output
typedef struct {
    gsl_vector **crossings_y_at_x_targets;
    size_t num_x_targets;
    int error_code;
} TrajectoryCrossings;

// For image_map function results
typedef struct {
    int miss;        // 1 if miss, 0 if hit
    double y_window; // y-coordinate on window plane
    double y_image;  // y-coordinate on source/image plane
    int error_code;
} ImageMapResult;

// Data point structure for photon mapping (saved to/loaded from chunks)
typedef struct {
    double y_window_cart_x;
    double y_window_cart_y;
    double y_image_cart_x;
    double y_image_cart_y;
} PhotonMapDataPoint;

// For PPM images
typedef struct {
    unsigned char r, g, b;
} RGBColor;

typedef struct {
    unsigned char *data;
    int width;
    int height;
    int channels; // Should always be 3 for P6
} PPMImage;

// For results_radial_light_ring output
typedef struct {
    char** saved_files_list;         // List of chunk filenames created
    gsl_vector* r_sample_misses;    // r_rho values that missed
    gsl_vector* r_sample_hits;      // r_rho values that hit
    int num_files_created;          // Actual count of file strings in saved_files_list
    int capacity_saved_files;       // Current capacity of saved_files_list array
} ResultsRadialLightRingOutput;


// --- Public Function Declarations ---

// From trajectory_module.c (or trajectory_api.c in earlier plans)

double normalize_phi(double phi);
int gsl_vector_dynamic_append(gsl_vector **vec_ptr, double value_param, size_t *current_capacity_ptr);

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
    const gsl_vector *x_targets, double t_end
);
void free_trajectory_crossings(TrajectoryCrossings *crossings);

// From image_mapping_utils.c
ImageMapResult image_map(
    double Y, double Z, double x_0_plane, double x_1_window, double x_2_observer, double a_scale
);
ImageMapResult image_map_radial(
    double r_rho, double x_0_plane, double x_1_window, double x_2_observer, double a_scale_factor
);

// From results_generation.c
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
void free_string_array(char** arr, int count); // General utility

// From photon_rendering.c
int map_photons(
    const char *image_source_path,
    PPMImage *image_source_data, // Parameter name can be simpler in .h
    char **photon_chunk_files,
    int num_chunk_files,         // Parameter name can be simpler in .h
    const char *save_path,
    double source_logical_width_input,
    int target_output_pixel_width,
    const double *observer_window_bounds,
    RGBColor default_color            // Parameter name can be simpler in .h
);

// From photon_chunk_io.c
int save_photon_data_chunk(const char *filename, const PhotonMapDataPoint* data, size_t num_points);
PhotonMapDataPoint* load_photon_data_chunk(const char *filename, size_t *num_points_read);

// From ppm_image_io.c
PPMImage* load_ppm_image(const char *filename);
int save_ppm_image(const char *filename, const PPMImage *image);
void free_ppm_image(PPMImage *image);

#endif // SCHWARZSCHILD_TRACER_H
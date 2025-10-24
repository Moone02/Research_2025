"""
Generates the complete C project for simulating massive particle geodesics.

This script follows the nrpy "gold standard" format. It sets up all necessary
parameters, retrieves symbolic recipes from the dedicated 'geodesics' equations
module, registers all required C functions, and assembles the final project
directory with a Makefile.

The generated C code simulates the N-body dynamics of a "collisionless gas"
(e.g., stars, gas clouds) orbiting a black hole, relying on the external
GNU Scientific Library (GSL) for robust ODE integration.

This script also includes a post-processing step to convert the raw binary
snapshot files into query-optimized .kdtree.bin files for use by the
photon geodesic integrator.

Author: Dalton J. Moone
"""

# Step 0.a: Import standard Python modules
import argparse
import os
import shutil
import glob
import time
from functools import partial
from multiprocessing import Pool, cpu_count
from typing import List, Tuple

# Step 0.b: Import third-party modules
import numpy as np
import sympy as sp
import numba
from tqdm import tqdm

# Step 0.c: Import nrpy core modules
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.params as par
import nrpy.indexedexp as ixp 
import nrpy.helpers.generic as gh

# Step 0.d: Import nrpy BHaH infrastructure modules
from nrpy.infrastructures.BHaH import (
    BHaH_defines_h,
    cmdline_input_and_parfiles,
    Makefile_helpers as Makefile,
)

# Step 0.e: Import the new, validated geodesics equations module.
#           This module encapsulates all the symbolic physics logic.
from nrpy.equations.general_relativity import geodesics as geo


# ##############################################################################
# PART 1: C CODE GENERATION (THE BUILDERS)
# ##############################################################################
#
# This section contains the Python functions that generate and register all the
# necessary C code for the project. Each function is a "builder" that takes
# symbolic expressions (the "recipes") from the 'geodesics' module and uses
# nrpy's C code generation tools to produce an optimized C function, which is
# then registered with the CFunction_dict for final project assembly.
#

# --------------------------------------------------------------------------
# Section 1.A: Tier 4 - Low-Level Workers & Dispatchers
# --------------------------------------------------------------------------


def g4dd_kerr_schild(g4DD_expressions: List[List[sp.Expr]]) -> None:
    """
    Generate and register the C worker function for the Kerr-Schild metric.

    This function takes the symbolic expressions for the 10 unique components
    of the Kerr-Schild metric tensor and generates an optimized C function
    that computes them at a given position.

    Args:
        g4DD_expressions: A list of lists containing the symbolic expressions
                          for the Kerr-Schild metric components.
    """
    list_of_g4DD_syms = [g4DD_expressions[i][j] for i in range(4) for j in range(i, 4)]
    list_of_g4DD_C_vars = [f"metric->g{i}{j}" for i in range(4) for j in range(i, 4)]
    cfc.register_CFunction(
        includes=["BHaH_defines.h"],
        desc="Computes the 10 unique components of the Kerr metric in Cartesian Kerr-Schild coords.",
        name="g4DD_kerr_schild",
        params="""const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const double y[4],
                  metric_struct *restrict metric""",
        body=ccg.c_codegen(
            list_of_g4DD_syms, list_of_g4DD_C_vars, enable_cse=True, cse_varprefix="KS"
        ),
        include_CodeParameters_h=True,
    )


def g4dd_schwarzschild_cartesian(g4DD_expressions: List[List[sp.Expr]]) -> None:
    """
    Generate and register the C worker for standard Schwarzschild metric.

    Args:
        g4DD_expressions: A list of lists containing the symbolic expressions
                          for the Schwarzschild metric components.
    """
    list_of_g4DD_syms = [g4DD_expressions[i][j] for i in range(4) for j in range(i, 4)]
    list_of_g4DD_C_vars = [f"metric->g{i}{j}" for i in range(4) for j in range(i, 4)]
    cfc.register_CFunction(
        includes=["BHaH_defines.h"],
        desc="Computes the 10 unique components of the Schwarzschild metric in Cartesian coords.",
        name="g4DD_schwarzschild_cartesian",
        params="""const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const double y[4],
                  metric_struct *restrict metric""",
        body=ccg.c_codegen(
            list_of_g4DD_syms,
            list_of_g4DD_C_vars,
            enable_cse=True,
            cse_varprefix="Schw",
        ),
        include_CodeParameters_h=True,
    )


def con_kerr_schild(Gamma4UDD_expressions: List[List[List[sp.Expr]]]) -> None:
    """
    Generate and register the C worker for the Kerr-Schild Christoffel symbols.

    Args:
        Gamma4UDD_expressions: A rank-3 list of lists of lists containing the
                               symbolic expressions for the Christoffel symbols.
    """
    list_of_Gamma_syms = [
        Gamma4UDD_expressions[i][j][k]
        for i in range(4)
        for j in range(4)
        for k in range(j, 4)
    ]
    conn_Gamma4UDD = ixp.declarerank3("conn->Gamma4UDD", dimension=4)
    list_of_Gamma_C_vars = [
        str(conn_Gamma4UDD[i][j][k])
        for i in range(4)
        for j in range(4)
        for k in range(j, 4)
    ]
    cfc.register_CFunction(
        includes=["BHaH_defines.h"],
        desc="Computes the 40 unique Christoffel symbols for the Kerr metric in Kerr-Schild coords.",
        name="con_kerr_schild",
        params="""const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const double y[4],
                  connection_struct *restrict conn""",
        body=ccg.c_codegen(
            list_of_Gamma_syms,
            list_of_Gamma_C_vars,
            enable_cse=True,
            cse_varprefix="conn",
        ),
        include_CodeParameters_h=True,
    )


def con_schwarzschild_cartesian(Gamma4UDD_expressions: List[List[List[sp.Expr]]]) -> None:
    """
    Generate and register the C worker for Schwarzschild Christoffel symbols.

    Args:
        Gamma4UDD_expressions: A rank-3 list of lists of lists containing the
                               symbolic expressions for the Christoffel symbols.
    """
    list_of_Gamma_syms = [
        Gamma4UDD_expressions[i][j][k]
        for i in range(4)
        for j in range(4)
        for k in range(j, 4)
    ]
    conn_Gamma4UDD = ixp.declarerank3("conn->Gamma4UDD", dimension=4)
    list_of_Gamma_C_vars = [
        str(conn_Gamma4UDD[i][j][k])
        for i in range(4)
        for j in range(4)
        for k in range(j, 4)
    ]
    cfc.register_CFunction(
        includes=["BHaH_defines.h"],
        desc="Computes the unique Christoffel symbols for the Schwarzschild metric in Cartesian coords.",
        name="con_schwarzschild_cartesian",
        params="""const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const double y[4],
                  connection_struct *restrict conn""",
        body=ccg.c_codegen(
            list_of_Gamma_syms,
            list_of_Gamma_C_vars,
            enable_cse=True,
            cse_varprefix="conn_schw",
        ),
        include_CodeParameters_h=True,
    )


def g4dd_metric() -> None:
    """
    Generate and register the C dispatcher for metric calculations.

    This function generates a C function `g4DD_metric()` that contains a switch
    statement. Based on the runtime `metric->type`, it calls the appropriate
    specialized worker function (e.g., `g4DD_kerr_schild`).
    """
    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        desc="Dispatcher to compute the 4-metric g_munu for the chosen analytic metric.",
        name="g4DD_metric",
        params="""const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric,
                  const double y[8],
                  metric_struct *restrict metric_out""",
        body=r"""
    const double y_pos[4] = {y[0], y[1], y[2], y[3]};
    switch(metric->type) {
        case Schwarzschild:
        case Kerr:
            g4DD_kerr_schild(commondata, params, y_pos, metric_out);
            break;
        case Schwarzschild_Standard:
            g4DD_schwarzschild_cartesian(commondata, params, y_pos, metric_out);
            break;
        default:
            fprintf(stderr, "Error: MetricType %d not supported in g4DD_metric().\n", metric->type);
            exit(1);
    }
    """,
    )


def connections() -> None:
    """
    Generate and register the C dispatcher for Christoffel symbol calculations.

    Similar to g4DD_metric(), this function generates a C dispatcher that calls
    the appropriate worker function based on the runtime `metric->type`.
    """
    cfc.register_CFunction(
        includes=[
            "BHaH_defines.h",
            "BHaH_function_prototypes.h",
            "stdio.h",
            "stdlib.h",
        ],
        desc="Dispatcher to compute Christoffel symbols for the chosen analytic metric.",
        name="connections",
        params="""const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric,
                  const double y[8],
                  connection_struct *restrict conn""",
        body=r"""
    const double y_pos[4] = {y[0], y[1], y[2], y[3]};
    switch(metric->type) {
        case Schwarzschild:
        case Kerr:
            con_kerr_schild(commondata, params, y_pos, conn);
            break;
        case Schwarzschild_Standard:
            con_schwarzschild_cartesian(commondata, params, y_pos, conn);
            break;
        default:
            fprintf(stderr, "Error: MetricType %d not supported in connections().\n", metric->type);
            exit(1);
    }
    """,
    )


# --------------------------------------------------------------------------
# Section 1.B: Tier 3 - Core Subsystems & Engines
# --------------------------------------------------------------------------


def calculate_ode_rhs_massive(geodesic_rhs_expressions: List[sp.Expr]) -> None:
    """
    Generate the C engine for the massive particle geodesic ODE right-hand sides.

    Args:
        geodesic_rhs_expressions: List of 8 symbolic expressions for the RHS.
    """
    cfc.register_CFunction(
        includes=["BHaH_defines.h"],
        desc="Calculates the RHS of the 8 massive particle geodesic ODEs.",
        name="calculate_ode_rhs_massive",
        params="""const double y[8],
                  const connection_struct *restrict conn,
                  double rhs_out[8]""",
        body=ccg.c_codegen(
            geodesic_rhs_expressions,
            [f"rhs_out[{i}]" for i in range(8)],
            enable_cse=True,
        ),
    )


def check_conservation_massive(
    kerr_E_expr: sp.Expr,
    kerr_L_exprs: List[sp.Expr],
    kerr_Q_expr: sp.Expr,
    schw_E_expr: sp.Expr,
    schw_L_exprs: List[sp.Expr],
    schw_Q_expr: sp.Expr,
) -> None:
    """
    Generate the C engine for computing conserved quantities (E, L_i, Q).

    This engine contains two distinct C code blocks, one for Kerr (using the
    full Carter constant) and one for Schwarzschild (where Q simplifies to L^2),
    and dispatches to the correct one at runtime.

    Args:
        kerr_E_expr: Symbolic expression for energy in Kerr.
        kerr_L_exprs: List of symbolic expressions for [Lx, Ly, Lz] in Kerr.
        kerr_Q_expr: Symbolic expression for the Carter Constant in Kerr.
        schw_E_expr: Symbolic expression for energy in Schwarzschild.
        schw_L_exprs: List of symbolic expressions for [Lx, Ly, Lz] in Schwarzschild.
        schw_Q_expr: Symbolic expression for L^2 in Schwarzschild.
    """
    output_vars = ["*E", "*Lx", "*Ly", "*Lz", "*Q"]
    kerr_expressions = [kerr_E_expr] + kerr_L_exprs + [kerr_Q_expr]
    schw_expressions = [schw_E_expr] + schw_L_exprs + [schw_Q_expr]

    body_kerr = ccg.c_codegen(
        kerr_expressions, output_vars, enable_cse=True, include_braces=False
    )
    body_schw = ccg.c_codegen(
        schw_expressions, output_vars, enable_cse=True, include_braces=False
    )

    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        desc="Computes conserved quantities (E, L_i, Q) for a given massive particle state vector.",
        name="check_conservation_massive",
        params="""const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric_params_in,
                  const double y[8],
                  double *E, double *Lx, double *Ly, double *Lz, double *Q""",
        body = f"""
    // Unpack parameters from commondata struct that are needed symbolically.
    const REAL a_spin = commondata->a_spin;

    // Allocate a metric struct on the HEAP to get a pointer.
    // This is necessary because the symbolic expressions for conserved
    // quantities were generated assuming pointer access (metric->g...).
    metric_struct* metric = (metric_struct*)malloc(sizeof(metric_struct));
    if (metric == NULL) {{
        fprintf(stderr, "Error: Failed to allocate memory for metric_struct in check_conservation.\\n");
        exit(1);
    }}
    
    // Call the dispatcher to fill the struct with metric components.
    // Note that we now pass the pointer 'metric' directly, not its address '&metric'.
    g4DD_metric(commondata, params, metric_params_in, y, metric);

    // Dispatch to the correct block of code based on metric type.
    if (metric_params_in->type == Kerr) {{
        {body_kerr}
    }} else {{ // Both Schwarzschild and Schwarzschild_Standard
        {body_schw}
    }}
    
    // Free the heap-allocated memory to prevent memory leaks.
    free(metric);
    """,
    )


def ode_gsl_wrapper_massive() -> None:
    """
    Generate the C wrapper function to interface with the GSL ODE solver.

    This function acts as a bridge, unpacking a generic `void*` pointer
    to access simulation data, calling our project-specific dispatchers and
    engines, and returning the computed derivatives to the GSL solver.
    """
    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h", "gsl/gsl_errno.h"],
        desc="GSL wrapper for the massive particle geodesic ODEs.",
        cfunc_type="int",
        name="ode_gsl_wrapper_massive",
        params="double t, const double y[8], double f[8], void *params",
        body=r"""
    (void)t; // Mark proper time 't' as unused to avoid compiler warnings.
    gsl_params *gsl_parameters = (gsl_params *)params;
    metric_struct g4DD;
    connection_struct conn;
    g4DD_metric(gsl_parameters->commondata, gsl_parameters->params, gsl_parameters->metric, y, &g4DD);
    connections(gsl_parameters->commondata, gsl_parameters->params, gsl_parameters->metric, y, &conn);
    calculate_ode_rhs_massive(y, &conn, f);
    return GSL_SUCCESS;
    """,
    )


# --------------------------------------------------------------------------
# Section 1.C: Tier 2 - Initial Condition & Integration Orchestrators
# --------------------------------------------------------------------------


def calculate_ut_uphi_from_r(ut_expr: sp.Expr, uphi_expr: sp.Expr) -> None:
    """
    Generate the C helper for computing u^t and u^phi from a radius.

    Args:
        ut_expr: The symbolic expression for u^t.
        uphi_expr: The symbolic expression for u^phi.
    """
    cfc.register_CFunction(
        includes=["BHaH_defines.h"],
        desc="Computes u^t and u^phi for a circular orbit at a given radius using a numerically stable method.",
        name="calculate_ut_uphi_from_r",
        params="""const double r_initial,
                  const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  double *ut, double *uphi""",
        body=ccg.c_codegen(
            [ut_expr, uphi_expr], ["*ut", "*uphi"], enable_cse=True, cse_varprefix="ic"
        ),
        include_CodeParameters_h=True,
    )


def set_initial_conditions_massive() -> None:
    """
    Generate the C engine to set the full initial 8-component state vector.

    This orchestrator takes a particle's initial position, calls the
    `calculate_ut_uphi_from_r` helper, and then transforms the resulting u^phi
    into Cartesian 4-velocity components u^x and u^y.
    """
    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h", "<math.h>"],
        desc="Sets the initial 8-component state vector for a massive particle.",
        name="set_initial_conditions_massive",
        params="""const particle_initial_state_t *restrict initial_state,
                  const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  double y_out[8]""",
        body=r"""
    const double x = initial_state->pos[1];
    const double y = initial_state->pos[2];
    const double z = initial_state->pos[3];
    const double r = sqrt(x*x + y*y + z*z);
    double ut, uphi;
    calculate_ut_uphi_from_r(r, commondata, params, &ut, &uphi);
    if (isnan(ut)) {
        for(int i=0; i<8; i++) y_out[i] = NAN;
        return;
    }
    y_out[0] = initial_state->pos[0];
    y_out[1] = x;
    y_out[2] = y;
    y_out[3] = z;
    y_out[4] = ut;
    y_out[5] = -y * uphi;
    y_out[6] =  x * uphi;
    y_out[7] = 0.0;
    """,
    )


def integrate_single_particle() -> None:
    """
    Generate the C orchestrator for integrating a single particle in production.

    This function evolves a particle from a start time to an end time using the
    high-level GSL driver, which ensures the particle state is returned at
    the exact end time. It also includes termination checks.
    """
    cfc.register_CFunction(
        includes=[
            "BHaH_defines.h",
            "BHaH_function_prototypes.h",
            "gsl/gsl_errno.h",
            "gsl/gsl_odeiv2.h",
            "<math.h>",
        ],
        desc="Integrates a single massive particle path between two coordinate times using the GSL driver.",
        cfunc_type="int",
        name="integrate_single_particle",
        params="""const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric,
                  const double t_start, const double t_end,
                  double y_in_out[8]""",
        body=r"""
    gsl_params gsl_parameters = {commondata, params, metric};
    gsl_odeiv2_system sys = {ode_gsl_wrapper_massive, NULL, 8, &gsl_parameters};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-11, 1e-11);
    double t = t_start;
    int status = gsl_odeiv2_driver_apply(d, &t, t_end, y_in_out);
    gsl_odeiv2_driver_free(d);
    if (status != GSL_SUCCESS) return 1;
    const double r_sq = y_in_out[1]*y_in_out[1] + y_in_out[2]*y_in_out[2] + y_in_out[3]*y_in_out[3];
    const double r_horizon = commondata->M_scale * (1.0 + sqrt(1.0 - commondata->a_spin*commondata->a_spin));
    if (r_sq < r_horizon*r_horizon || r_sq > commondata->r_escape*commondata->r_escape || fabs(y_in_out[4]) > commondata->ut_max) {
        return 1;
    }
    return 0;
    """,
        include_CodeParameters_h=True,
    )


def integrate_single_particle_DEBUG() -> None:
    """
    Generate the C orchestrator for a single-particle debug run.

    This function traces the full path of a single particle, writing its state
    at each internal step of the integrator to a text file for validation.
    """
    cfc.register_CFunction(
        includes=[
            "BHaH_defines.h",
            "BHaH_function_prototypes.h",
            "gsl/gsl_errno.h",
            "gsl/gsl_odeiv2.h",
        ],
        desc="DEBUG integrator for a single massive particle. Writes the full trajectory to a text file.",
        name="integrate_single_particle_DEBUG",
        params="""const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric,
                  const double start_y[8],
                  double final_y_state[8]""",
        body=r"""
    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step * step = gsl_odeiv2_step_alloc(T, 8);
    gsl_odeiv2_control * control = gsl_odeiv2_control_yp_new(1e-14, 1e-14);
    gsl_odeiv2_evolve * evol = gsl_odeiv2_evolve_alloc(8);
    gsl_params gsl_parameters = {commondata, params, metric};
    gsl_odeiv2_system sys = {ode_gsl_wrapper_massive, NULL, 8, &gsl_parameters};
    double y_c[8];
    double t = 0.0, dt = 0.01;
    for (int j = 0; j < 8; j++) { y_c[j] = start_y[j]; }
    FILE *fp = fopen("massive_particle_path.txt", "w");
    if (fp == NULL) {
        fprintf(stderr, "Error: Could not open massive_particle_path.txt for writing.\n");
        exit(1);
    }
    fprintf(fp, "# ProperTime_tau\tCoordTime_t\tx\ty\tz\tu^t\tu^x\tu^y\tu^z\n");
    printf("Starting debug trace for single massive particle...\n");
    printf("Step | Proper Time (τ) | Coord Time (t) |      x     |      y     |      z     |      u^t   \n");
    printf("-------------------------------------------------------------------------------------------\n");
    for (int i = 0; i < 2000000; i++) {
        int status = gsl_odeiv2_evolve_apply(evol, control, step, &sys, &t, 1e10, &dt, y_c);
        fprintf(fp, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                t, y_c[0], y_c[1], y_c[2], y_c[3], y_c[4], y_c[5], y_c[6], y_c[7]);
        if (i % 500 == 0) {
            printf("%4d | %15.4e | %14.4f | %10.4f | %10.4f | %10.4f | %10.4f\n",
                   i, t, y_c[0], y_c[1], y_c[2], y_c[3], y_c[4]);
        }
        const double r_sq = y_c[1]*y_c[1] + y_c[2]*y_c[2] + y_c[3]*y_c[3];
        const double r_horizon = commondata->M_scale * (1.0 + sqrt(1.0 - commondata->a_spin*commondata->a_spin));
        if (status != GSL_SUCCESS) { printf("Termination: GSL ERROR (status = %d)\n", status); break; }
        if (r_sq < r_horizon*r_horizon) { printf("Termination: Fell below event horizon (r=%.2f)\n", sqrt(r_sq)); break; }
        if (r_sq > commondata->r_escape*commondata->r_escape) { printf("Termination: Escaped to r > %.1f\n", commondata->r_escape); break; }
        if (fabs(y_c[4]) > commondata->ut_max) { printf("Termination: Runaway u^t > %.1e\n", commondata->ut_max); break; }
        if (fabs(y_c[0]) > commondata->t_max_integration) { printf("Termination: Exceeded max integration time t > %.1f\n", commondata->t_max_integration); break; }
    }
    for(int j=0; j<8; j++) { final_y_state[j] = y_c[j]; }
    fclose(fp);
    gsl_odeiv2_evolve_free(evol);
    gsl_odeiv2_control_free(control);
    gsl_odeiv2_step_free(step);
    """,
        include_CodeParameters_h=True,
    )


# --------------------------------------------------------------------------
# Section 1.D: Tier 1 & 0 - High-Level Orchestrators & main()
# --------------------------------------------------------------------------


def generate_disk_initial_conditions() -> None:
    """
    Generate a C orchestrator to create a simple, axisymmetric Keplerian disk.

    This C function populates a large array of initial state vectors by looping
    through a grid of radii `r` and azimuthal angles `phi`.
    """
    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        desc="Generates the complete y[8] initial state for all particles in a Keplerian disk.",
        cfunc_type="int",
        name="generate_disk_initial_conditions",
        params="""const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  double *restrict y_initial_states""",
        body=r"""
    int particle_count = 0;
    const double dr = (commondata->disk_num_r > 1) ? (commondata->disk_r_max - commondata->disk_r_min) / (commondata->disk_num_r - 1) : 0;
    for (int i = 0; i < commondata->disk_num_r; i++) {
        const double r = commondata->disk_r_min + i * dr;
        const int num_phi_at_r = (commondata->disk_num_phi > 1) ? (int)(commondata->disk_num_phi * (r / commondata->disk_r_max)) : 1;
        if (num_phi_at_r == 0) continue;
        const double dphi = 2.0 * M_PI / num_phi_at_r;
        double ut_at_r, uphi_at_r;
        calculate_ut_uphi_from_r(r, commondata, params, &ut_at_r, &uphi_at_r);
        for (int j = 0; j < num_phi_at_r; j++) {
            const double phi = j * dphi;
            double *y = &y_initial_states[particle_count * 8];
            y[0] = 0.0; y[1] = r * cos(phi); y[2] = r * sin(phi); y[3] = 0.0;
            y[4] = ut_at_r; y[5] = -y[2] * uphi_at_r; y[6] = y[1] * uphi_at_r; y[7] = 0.0;
            particle_count++;
        }
    }
    return particle_count;
    """,
    )


def generate_spiral_galaxy_initial_conditions() -> None:
    """
    Generate a C orchestrator to create a "grand design" spiral galaxy.

    This function uses a random placement strategy based on a logarithmic
    spiral formula to create a spiral galaxy pattern.
    """
    cfc.register_CFunction(
        includes=[
            "BHaH_defines.h",
            "BHaH_function_prototypes.h",
            "<math.h>",
            "<stdlib.h>",
        ],
        desc="Generates the complete y[8] initial state for all particles in a spiral galaxy disk.",
        cfunc_type="int",
        name="generate_spiral_galaxy_initial_conditions",
        params="""const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  double *restrict y_initial_states""",
        body=r"""
    int particle_count = 0;
    const int num_particles_total = commondata->disk_num_r * commondata->disk_num_phi;
    const int num_arms = commondata->spiral_galaxy_num_arms;
    const double arm_tightness = commondata->spiral_galaxy_arm_tightness;
    srand(42);
    for (int i = 0; i < num_particles_total; i++) {
        const double r = commondata->disk_r_min + (commondata->disk_r_max - commondata->disk_r_min) * sqrt((double)rand() / RAND_MAX);
        const double theta_base = (1.0 / arm_tightness) * log(r / commondata->disk_r_min);
        const int arm_index = rand() % num_arms;
        const double arm_offset = (2.0 * M_PI / num_arms) * arm_index;
        const double phi_spread = (M_PI / num_arms) * 0.2 * ((double)rand() / RAND_MAX - 0.5);
        const double phi = theta_base + arm_offset + phi_spread;
        double ut_at_r, uphi_at_r;
        calculate_ut_uphi_from_r(r, commondata, params, &ut_at_r, &uphi_at_r);
        double *y = &y_initial_states[particle_count * 8];
        y[0] = 0.0; y[1] = r * cos(phi); y[2] = r * sin(phi); y[3] = 0.0;
        y[4] = ut_at_r; y[5] = -y[2] * uphi_at_r; y[6] = y[1] * uphi_at_r; y[7] = 0.0;
        particle_count++;
    }
    return particle_count;
    """,
    )


def generate_barred_flocculent_spiral_ic() -> None:
    """
    Generate the C orchestrator for a barred, clumpy ("flocculent") spiral.

    This function uses rejection sampling to create a complex geometry with a
    dense central bulge, a rectangular bar, and clumpy spiral arms.
    """
    cfc.register_CFunction(
        includes=[
            "BHaH_defines.h",
            "BHaH_function_prototypes.h",
            "<math.h>",
            "<stdlib.h>",
            "<time.h>",
        ],
        desc="Generates the initial state for all particles in a barred flocculent spiral galaxy.",
        cfunc_type="int",
        name="generate_barred_flocculent_spiral_ic",
        params="""const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  double *restrict y_initial_states""",
        body=r"""
    unsigned int seed = 42;
    int particle_count = 0;
    const int num_particles_total = commondata->disk_num_r * commondata->disk_num_phi;
    const double r_min = commondata->disk_r_min, r_max = commondata->disk_r_max;
    const double bar_len = commondata->bar_length, bar_width = commondata->bar_length * commondata->bar_aspect_ratio;
    const double bulge_rad = commondata->bulge_radius;
    const int num_arms = commondata->spiral_galaxy_num_arms;
    const double arm_tightness = commondata->spiral_galaxy_arm_tightness;
    for (int i = 0; i < num_particles_total; i++) {
        double r, phi, x, y;
        seed += i;
        while (1) {
            double random_val = (double)rand_r(&seed) / RAND_MAX;
            r = sqrt(random_val * (r_max*r_max - r_min*r_min) + r_min*r_min);
            phi = 2.0 * M_PI * ((double)rand_r(&seed) / RAND_MAX);
            x = r * cos(phi); y = r * sin(phi);
            double acceptance_prob = 0.0;
            if (r < bulge_rad) {
                acceptance_prob = commondata->bulge_density_factor * commondata->arm_particle_density;
            } else if (fabs(x) < bar_len / 2.0 && fabs(y) < bar_width / 2.0) {
                acceptance_prob = commondata->bar_density_factor * commondata->arm_particle_density;
            } else {
                if (((double)rand_r(&seed) / RAND_MAX) < commondata->arm_particle_density) {
                    double theta_base = (1.0 / arm_tightness) * log(r / r_min);
                    double delta_phi_raw = phi - theta_base;
                    double angle_between_arms = 2.0 * M_PI / num_arms;
                    double delta_phi_folded = fmod(delta_phi_raw, angle_between_arms);
                    if (delta_phi_folded > 0.5 * angle_between_arms) delta_phi_folded -= angle_between_arms;
                    else if (delta_phi_folded < -0.5 * angle_between_arms) delta_phi_folded += angle_between_arms;
                    double arm_profile = exp(-0.5 * SQR(delta_phi_folded / (arm_tightness * 0.5)));
                    double clump_modulation = cos(commondata->arm_clumpiness_factor * (theta_base - phi));
                    clump_modulation *= cos(num_arms * phi * commondata->arm_clump_size);
                    acceptance_prob = arm_profile + 0.5 * clump_modulation;
                }
            }
            if (((double)rand_r(&seed) / RAND_MAX) < acceptance_prob) break;
        }
        double ut_at_r, uphi_at_r;
        calculate_ut_uphi_from_r(r, commondata, params, &ut_at_r, &uphi_at_r);
        double *y_state = &y_initial_states[particle_count * 8];
        y_state[0] = 0.0; y_state[1] = x; y_state[2] = y; y_state[3] = 0.0;
        y_state[4] = ut_at_r; y_state[5] = -y * uphi_at_r; y_state[6] = x * uphi_at_r; y_state[7] = 0.0;
        particle_count++;
    }
    return particle_count;
    """,
    )


def run_mass_integrator_production() -> None:
    """
    Generate the top-level orchestrator for a full production run.

    This C function sets up the particle ensemble, then enters the main time
    loop. At each step, it saves a snapshot and then calls
    `integrate_single_particle` for each active particle in a parallel loop.
    """
    cfc.register_CFunction(
        includes=[
            "BHaH_defines.h",
            "BHaH_function_prototypes.h",
            "<math.h>",
            "<sys/stat.h>",
            "<string.h>",
            "<stdio.h>",
            "<stdlib.h>",
        ],
        desc="Orchestrates the full production run for the mass integrator.",
        name="run_mass_integrator_production",
        params="""const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric""",
        body=r"""
    const int max_particles = commondata->disk_num_r * commondata->disk_num_phi;
    double *y_initial_states = (double *)malloc(sizeof(double) * 8 * max_particles);
    if (!y_initial_states) { fprintf(stderr, "Error: Failed to allocate memory for initial states.\n"); exit(1); }
    int num_particles = 0;
    printf("Generating initial conditions of type: %s\n", commondata->initial_conditions_type);
    if (strcmp(commondata->initial_conditions_type, "SpiralGalaxy") == 0) {
        num_particles = generate_spiral_galaxy_initial_conditions(commondata, params, y_initial_states);
    } else if (strcmp(commondata->initial_conditions_type, "BarredFlocculentSpiral") == 0) {
        num_particles = generate_barred_flocculent_spiral_ic(commondata, params, y_initial_states);
    } else {
        if (strcmp(commondata->initial_conditions_type, "KeplerianDisk") != 0) {
            printf("Warning: Unrecognized type '%s'. Defaulting to KeplerianDisk.\n", commondata->initial_conditions_type);
        }
        num_particles = generate_disk_initial_conditions(commondata, params, y_initial_states);
    }
    printf("Generated %d particles.\n", num_particles);
    mass_particle_state_t *particle_states = (mass_particle_state_t *)malloc(sizeof(mass_particle_state_t) * num_particles);
    if (!particle_states) { fprintf(stderr, "Error: Failed to allocate memory for particle states.\n"); exit(1); }
    for (int i=0; i<num_particles; i++) {
        double *y_start = &y_initial_states[i*8];
        particle_states[i].id = i;
        particle_states[i].pos[0] = y_start[1]; particle_states[i].pos[1] = y_start[2]; particle_states[i].pos[2] = y_start[3];
        particle_states[i].u[0] = y_start[4]; particle_states[i].u[1] = y_start[5]; particle_states[i].u[2] = y_start[6]; particle_states[i].u[3] = y_start[7];
        const double r = sqrt(y_start[1]*y_start[1] + y_start[2]*y_start[2]);
        particle_states[i].lambda_rest = commondata->disk_lambda_rest_at_r_min * pow(r / commondata->disk_r_min, 0.75);
        particle_states[i].j_intrinsic = (float)pow(r / commondata->disk_r_min, -3.0);
    }
    free(y_initial_states);
    mkdir(commondata->output_folder, 0755);
    int snapshot_count = 0;
    for (double current_t = 0; current_t <= commondata->t_final; current_t += commondata->snapshot_every_t) {
        char filename[200];
        snprintf(filename, 200, "%s/mass_blueprint_t_%04d.bin", commondata->output_folder, snapshot_count);
        printf("Saving snapshot: %s (t=%.2f)\n", filename, current_t);
        FILE *fp_out = fopen(filename, "wb");
        if (!fp_out) { exit(1); }
        int active_particles = 0;
        for(int i=0; i<num_particles; i++) if (!isnan(particle_states[i].pos[0])) active_particles++;
        fwrite(&active_particles, sizeof(int), 1, fp_out);
        for (int i=0; i<num_particles; i++) if (!isnan(particle_states[i].pos[0])) fwrite(&particle_states[i], sizeof(mass_particle_state_t), 1, fp_out);
        fclose(fp_out);
        snapshot_count++;
        if (current_t >= commondata->t_final) break;
        const double t_next_snapshot = current_t + commondata->snapshot_every_t;
        #pragma omp parallel for
        for (int i = 0; i < num_particles; i++) {
            if (isnan(particle_states[i].pos[0])) continue;
            double y_particle[8];
            y_particle[0] = current_t;
            y_particle[1] = particle_states[i].pos[0]; y_particle[2] = particle_states[i].pos[1]; y_particle[3] = particle_states[i].pos[2];
            y_particle[4] = particle_states[i].u[0]; y_particle[5] = particle_states[i].u[1]; y_particle[6] = particle_states[i].u[2]; y_particle[7] = particle_states[i].u[3];
            int status = integrate_single_particle(commondata, params, metric, y_particle[0], t_next_snapshot, y_particle);
            if (status != 0) {
                particle_states[i].pos[0] = NAN;
            } else {
                particle_states[i].pos[0] = y_particle[1]; particle_states[i].pos[1] = y_particle[2]; particle_states[i].pos[2] = y_particle[3];
                particle_states[i].u[0] = y_particle[4]; particle_states[i].u[1] = y_particle[5]; particle_states[i].u[2] = y_particle[6]; particle_states[i].u[3] = y_particle[7];
            }
        }
    }
    free(particle_states);
    """,
    )


def main_c_function() -> None:
    """
    Generate the main() C function, the entry point for the executable.

    This function acts as the master orchestrator. It initializes parameters,
    determines the run mode (debug or production), and dispatches to the
    appropriate top-level orchestrator. It also handles the post-processing
    step of calling this script to generate k-d trees.
    """
    debug_body = r"""
        particle_initial_state_t initial_state;
        initial_state.id = 0;
        const char *filename = "particle_debug_initial_conditions.txt";
        FILE *fp_in = fopen(filename, "r");
        if (fp_in == NULL) {
            printf("File '%s' not found. Creating it with default values.\n", filename);
            fp_in = fopen(filename, "w");
            if (fp_in == NULL) { fprintf(stderr, "Error: Could not create '%s'.\n", filename); return 1; }
            fprintf(fp_in, "# Format: t_initial pos_x pos_y pos_z u_x u_y u_z\n");
            fprintf(fp_in, "0.0 10.0  0.0  0.0   0.0  0.363232  0.0\n");
            fclose(fp_in);
            fp_in = fopen(filename, "r");
            if (fp_in == NULL) { fprintf(stderr, "Error: Could not re-open '%s'.\n", filename); return 1; }
        }
        char line[256];
        while (fgets(line, sizeof(line), fp_in)) {
            if (line[0] != '#') {
                sscanf(line, "%lf %lf %lf %lf %lf %lf %lf",
                       &initial_state.pos[0], &initial_state.pos[1], &initial_state.pos[2], &initial_state.pos[3],
                       &initial_state.u_spatial[0], &initial_state.u_spatial[1], &initial_state.u_spatial[2]);
                break;
            }
        }
        fclose(fp_in);
        printf("--- Single Particle Debug Run ---\n");
        printf("  pos = (t=%.4f, x=%.4f, y=%.4f, z=%.4f)\n", initial_state.pos[0], initial_state.pos[1], initial_state.pos[2], initial_state.pos[3]);
        double y_start[8], y_final[8];
        set_initial_conditions_massive(&initial_state, &commondata, &params, y_start);
        printf("\nInitial State Vector (y_start):\n");
        printf("  t=%.2f, x=%.2f, y=%.2f, z=%.2f\n", y_start[0], y_start[1], y_start[2], y_start[3]);
        printf("  u^t=%.4f, u^x=%.4f, u^y=%.4f, u^z=%.4f\n\n", y_start[4], y_start[5], y_start[6], y_start[7]);
        if(commondata.perform_conservation_check) {
            double E_i, Lx_i, Ly_i, Lz_i, Q_i, E_f, Lx_f, Ly_f, Lz_f, Q_f;
            check_conservation_massive(&commondata, &params, &metric, y_start, &E_i, &Lx_i, &Ly_i, &Lz_i, &Q_i);
            integrate_single_particle_DEBUG(&commondata, &params, &metric, y_start, y_final);
            check_conservation_massive(&commondata, &params, &metric, y_final, &E_f, &Lx_f, &Ly_f, &Lz_f, &Q_f);
        } else {
            integrate_single_particle_DEBUG(&commondata, &params, &metric, y_start, y_final);
        }
        printf("\nDebug run finished. Trajectory saved to 'massive_particle_path.txt'.\n");
    """
    cfc.register_CFunction(
        includes=[
            "BHaH_defines.h",
            "BHaH_function_prototypes.h",
            "<string.h>",
            "<stdio.h>",
            "<stdlib.h>",
        ],
        desc="Main entry point for the massive particle geodesic integrator.",
        cfunc_type="int",
        name="main",
        params="int argc, const char *argv[]",
        body=f"""
    commondata_struct commondata;
    params_struct params;
    metric_params metric;
    commondata_struct_set_to_default(&commondata);
    cmdline_input_and_parfile_parser(&commondata, argc, argv);
    metric.type = (commondata.a_spin == 0.0) ? Schwarzschild : Kerr;
    if (commondata.run_in_debug_mode) {{
        {debug_body}
    }} else {{
        printf("----------------------------------------\\n");
        printf("Massive Particle Integrator\\n");
        printf("----------------------------------------\\n");
        printf("Metric Settings: %s (a=%.2f, M=%.2f)\\n", (metric.type == Kerr) ? "Kerr" : "Schwarzschild", commondata.a_spin, commondata.M_scale);
        printf("Integration Settings: Max Time = %.1f M, Snapshot Every = %.1f M\\n", commondata.t_final, commondata.snapshot_every_t);
        printf("Initial Conditions: Type = %s, Particles = %d x %d, Radius = %.2f to %.2f M\\n", commondata.initial_conditions_type, commondata.disk_num_r, commondata.disk_num_phi, commondata.disk_r_min, commondata.disk_r_max);
        printf("----------------------------------------\\n\\n");
        run_mass_integrator_production(&commondata, &params, &metric);
        printf("\\nProduction run finished successfully.\\n");
        if (commondata.generate_kdtree_files) {{
            printf("\\n--- Starting Post-Processing: Generating K-d Tree files ---\\n");
            char command[1024];
            snprintf(command, sizeof(command), "python3 ../../mass_geodesic_integrator.py --post-process --outdir ..");
            int status = system(command);
            if (status != 0) {{
                fprintf(stderr, "Error: Post-processing script failed with exit code %d.\\n", status);
            }} else {{
                printf("--- Post-Processing Complete ---\\n");
            }}
        }}
    }}
    return 0;
    """,
    )


# ##############################################################################
# PART 2: POST-PROCESSING PYTHON LOGIC (K-D TREE GENERATION)
# ##############################################################################
#
# The following functions are executed by this script when it is called
# in post-processing mode (i.e., with the --post-process flag). They are
# responsible for converting the .bin snapshot files into .kdtree.bin files.
#


@numba.jit(nopython=True)
def build_implicit_kdtree_kernel(
    particle_pos_array: np.ndarray, num_particles: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Numba-compiled kernel to build the implicit k-d tree.

    Args:
        particle_pos_array: A NumPy array of shape (N, 3) with particle positions.
        num_particles: The number of particles.

    Returns:
        A tuple (reordering_indices, node_metadata).
    """
    reordering_indices = np.zeros(num_particles, dtype=np.int32)
    node_metadata = np.full(num_particles, -1, dtype=np.int32)
    q_target_idx = np.zeros(num_particles, dtype=np.int32)
    q_slice_start = np.zeros(num_particles, dtype=np.int32)
    q_slice_end = np.zeros(num_particles, dtype=np.int32)
    q_head = q_tail = 0
    indices_array = np.arange(num_particles, dtype=np.int32)
    if num_particles > 0:
        q_target_idx[q_tail], q_slice_start[q_tail], q_slice_end[q_tail] = (
            0,
            0,
            num_particles,
        )
        q_tail += 1
    while q_head < q_tail:
        target_idx, start, end = (
            q_target_idx[q_head],
            q_slice_start[q_head],
            q_slice_end[q_head],
        )
        q_head += 1
        current_indices = indices_array[start:end]
        current_positions = particle_pos_array[current_indices]
        min_coord, max_coord = current_positions[0], current_positions[0]
        for i in range(1, len(current_positions)):
            min_coord = np.minimum(min_coord, current_positions[i])
            max_coord = np.maximum(max_coord, current_positions[i])
        split_axis = np.argmax(max_coord - min_coord)
        sorted_indices = current_indices[
            np.argsort(particle_pos_array[current_indices, split_axis])
        ]
        indices_array[start:end] = sorted_indices
        median_offset = len(sorted_indices) // 2
        reordering_indices[target_idx] = sorted_indices[median_offset]
        node_metadata[target_idx] = split_axis
        left_child_idx, right_child_idx = 2 * target_idx + 1, 2 * target_idx + 2
        if median_offset > 0 and left_child_idx < num_particles:
            q_target_idx[q_tail], q_slice_start[q_tail], q_slice_end[q_tail] = (
                left_child_idx,
                start,
                start + median_offset,
            )
            q_tail += 1
        if median_offset + 1 < len(sorted_indices) and right_child_idx < num_particles:
            q_target_idx[q_tail], q_slice_start[q_tail], q_slice_end[q_tail] = (
                right_child_idx,
                start + median_offset + 1,
                end,
            )
            q_tail += 1
    return reordering_indices, node_metadata


def build_and_save_kdtree_snapshot_fast(
    raw_snapshot_file: str, output_dir: str
) -> None:
    """
    Driver function that handles file I/O and calls the fast Numba kernel.

    Args:
        raw_snapshot_file: Path to the input mass_blueprint .bin file.
        output_dir: Directory to save the output .kdtree.bin file.
    """
    raw_dtype = np.dtype(
        [
            ("id", np.int32),
            ("pos", "f8", (3,)),
            ("u", "f8", (4,)),
            ("lambda_rest", "f8"),
            ("j_intrinsic", "f4"),
        ]
    )
    with open(raw_snapshot_file, "rb") as f:
        num_particles = np.fromfile(f, dtype=np.int32, count=1)[0]
        if num_particles == 0:
            print(
                f"  Snapshot {os.path.basename(raw_snapshot_file)} is empty, skipping."
            )
            return
        particle_data = np.fromfile(f, dtype=raw_dtype, count=num_particles)
    particle_positions = np.ascontiguousarray(particle_data["pos"])
    reordering_map, node_metadata = build_implicit_kdtree_kernel(
        particle_positions, num_particles
    )
    reordered_particles = particle_data[reordering_map]
    output_filename = os.path.join(
        output_dir, os.path.basename(raw_snapshot_file).replace(".bin", ".kdtree.bin")
    )
    with open(output_filename, "wb") as f:
        f.write(np.uint64(num_particles))
        f.write(np.uint64(3))
        f.write(node_metadata.tobytes())
        f.write(reordered_particles.tobytes())


def run_kdtree_preprocessor_parallel(
    snapshot_input_dir: str, snapshot_output_dir: str
) -> None:
    """
    Find all raw .bin snapshot files and process them in parallel.

    Args:
        snapshot_input_dir: The directory containing the raw .bin files.
        snapshot_output_dir: The directory where .kdtree.bin files will be saved.
    """
    print("--- Starting K-d Tree Pre-processing (Parallel) ---")
    print(f"Input directory:  '{snapshot_input_dir}'")
    print(f"Shared Output directory: '{snapshot_output_dir}'")
    os.makedirs(snapshot_output_dir, exist_ok=True)
    snapshot_files = sorted(
        glob.glob(os.path.join(snapshot_input_dir, "mass_blueprint_t_*.bin"))
    )
    if not snapshot_files:
        print(
            "\nWARNING: No raw snapshot files found. Did you run the mass_integrator C code first?"
        )
        return
    num_processes = cpu_count()
    print(
        f"\nFound {len(snapshot_files)} files. Processing in parallel using {num_processes} CPU cores..."
    )
    task_function = partial(
        build_and_save_kdtree_snapshot_fast, output_dir=snapshot_output_dir
    )
    start_time = time.time()
    with Pool(processes=num_processes) as pool:
        list(
            tqdm(
                pool.imap_unordered(task_function, snapshot_files),
                total=len(snapshot_files),
            )
        )
    end_time = time.time()
    print(
        f"\n--- All snapshots processed successfully in {end_time - start_time:.2f} seconds. ---"
    )
    print(
        f"The query-ready .kdtree.bin files are now in the shared directory: '{snapshot_output_dir}'."
    )


# ##############################################################################
# PART 3: FINAL PROJECT ASSEMBLY
# ##############################################################################


def register_custom_structures_and_params() -> None:
    """
    Generate C code for all custom structs and enums, then register them.

    This function defines the C `typedef`s for all custom data structures
    and registers this block of C code to be written into `BHaH_defines.h`.
    """
    metric_components = [f"g{nu}{mu}" for nu in range(4) for mu in range(nu, 4)]
    metric_struct_str = (
        "typedef struct { double " + "; double ".join(metric_components) + "; } metric_struct;"
    )
    connection_components = [
        f"Gamma4UDD{i}{j}{k}" for i in range(4) for j in range(4) for k in range(j, 4)
    ]
    connections_struct_str = (
        "typedef struct { double "
        + "; double ".join(connection_components)
        + "; } connection_struct;"
    )
    other_structs = r"""
typedef enum { Schwarzschild, Kerr, Schwarzschild_Standard } Metric_t;
typedef struct { Metric_t type; } metric_params;
typedef struct {
    const commondata_struct *commondata;
    const params_struct *params;
    const metric_params *metric;
} gsl_params;
typedef struct {
    int id; double pos[4]; double u_spatial[3];
} particle_initial_state_t;
typedef struct {
    int id; double pos[3]; double u[4]; double lambda_rest; float j_intrinsic;
} __attribute__((packed)) mass_particle_state_t;
"""
    BHaH_defines_h.register_BHaH_defines(
        "data_structures",
        f"{metric_struct_str}\n{connections_struct_str}\n{other_structs}",
    )

# ##############################################################################
# PART 4: MAIN EXECUTION BLOCK
# ##############################################################################

if __name__ == "__main__":
    # Step 4.A: Set up command-line argument parser.
    # This allows the script to be run in one of two modes:
    # 1. Default mode: Generates the entire C code project.
    # 2. Post-processing mode (--post-process): Skips C code generation and
    #    only runs the Python-based k-d tree pre-processor on existing data.
    parser = argparse.ArgumentParser(
        description="Generate C code for the mass geodesic integrator or run post-processing."
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="project",
        help="The parent directory where the C project and processed data will be generated.",
    )
    parser.add_argument(
        "--post-process",
        action="store_true",
        help="Run only the post-processing (k-d tree generation) step on existing snapshot files.",
    )
    args = parser.parse_args()

    if args.post_process:
        # --- POST-PROCESSING-ONLY MODE ---
        # In this mode, the script only generates the k-d tree files from existing
        # .bin snapshots. It does not generate any C code.
        project_dir = args.outdir
        mass_project_name = "mass_integrator"
        snapshot_input_dir = os.path.join(project_dir, mass_project_name, "output")
        snapshot_output_dir = os.path.join(project_dir, "processed_snapshots")
        run_kdtree_preprocessor_parallel(snapshot_input_dir, snapshot_output_dir)

    else:
        # --- FULL C-CODE GENERATION MODE ---
        # This is the standard mode, which generates the complete C project.

        # Step 4.B: Set up project directories and core nrpy parameters.
        project_name = "mass_integrator"
        project_dir = os.path.join(args.outdir, project_name)

        # Clean the project directory to ensure a fresh build.
        shutil.rmtree(project_dir, ignore_errors=True)
        os.makedirs(project_dir, exist_ok=True)

        # Set the BHaH infrastructure for nrpy.
        par.set_parval_from_str("Infrastructure", "BHaH")

        # Step 4.C: Register all CodeParameters for the mass integrator.
        # This dictionary-based approach is clean and maintainable.
        param_groups = {
            "REAL": {
                "M_scale": 1.0,
                "a_spin": 0.9,
                "t_max_integration": 2000.0,
                "r_escape": 1500.0,
                "ut_max": 1e3,
                "disk_lambda_rest_at_r_min": 656.3,
                "disk_r_min": 6.0,
                "disk_r_max": 25.0,
                "snapshot_every_t": 10.0,
                "t_final": 2000.0,
                "bar_length": 5.0,
                "bar_aspect_ratio": 0.25,
                "bulge_radius": 1.5,
                "arm_particle_density": 0.3,
                "arm_clumpiness_factor": 8.0,
                "arm_clump_size": 0.5,
                "bar_density_factor": 2.0,
                "bulge_density_factor": 3.0,
                "spiral_galaxy_arm_tightness": 0.2,
            },
            "int": {
                "disk_num_r": 100,
                "disk_num_phi": 200,
                "spiral_galaxy_num_arms": 2,
            },
            "bool": {
                "generate_kdtree_files": True,
                "perform_conservation_check": True,
                "run_in_debug_mode": True,
            },
            "char[100]": {
                "output_folder": "output",
                "initial_conditions_type": "KeplerianDisk",
            },
        }
        for c_type, params in param_groups.items():
            for name, default in params.items():
                par.register_CodeParameter(
                    c_type, __name__, name, default, commondata=True, add_to_parfile=True
                )

        # Step 4.D: Execute the symbolic workflow by retrieving expressions from the module.
        print(" -> Generating symbolic expressions for geodesic equations...")
        kerr_eqs = geo.Geodesic_Equations["KerrSchild_massive"]
        schw_eqs = geo.Geodesic_Equations["Schwarzschild_Cartesian_massive"]
        print("    ... Symbolic expressions generated successfully.")

        # Step 4.E: Generate C code for parameter handling.
        # This must be done after all parameters are registered but before any
        # C functions that use them are generated, as they create the
        # set_CodeParameters.h header file.
        print(" -> Generating C code for parameter handling...")
        # We need to import CPs here because these are BHaH infrastructure-specific functions.
        from nrpy.infrastructures.BHaH import CodeParameters as CPs
        CPs.write_CodeParameters_h_files(project_dir=project_dir)
        CPs.register_CFunctions_params_commondata_struct_set_to_default()
        cmdline_input_and_parfiles.generate_default_parfile(
            project_dir=project_dir, project_name=project_name
        )
        cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
            project_name=project_name,
            cmdline_inputs=list(param_groups["REAL"].keys())
            + list(param_groups["int"].keys())
            + list(param_groups["bool"].keys())
            + list(param_groups["char[100]"].keys()),
        )
        print("    ... Parameter handling code generated.")

        # Step 4.F: Call all C-generating Python functions ("builders") to register C code.
        print(" -> Registering C functions for project...")
        g4dd_kerr_schild(kerr_eqs.g4DD)
        g4dd_schwarzschild_cartesian(schw_eqs.g4DD)
        con_kerr_schild(kerr_eqs.Gamma4UDD)
        con_schwarzschild_cartesian(schw_eqs.Gamma4UDD)
        g4dd_metric()
        connections()
        calculate_ode_rhs_massive(kerr_eqs.geodesic_rhs)
        check_conservation_massive(
            kerr_eqs.E_expr,
            kerr_eqs.L_exprs,
            kerr_eqs.Q_expr,
            schw_eqs.E_expr,
            schw_eqs.L_exprs,
            schw_eqs.Q_expr,
        )
        ode_gsl_wrapper_massive()
        calculate_ut_uphi_from_r(kerr_eqs.ut_stable_expr, kerr_eqs.uphi_stable_expr)
        set_initial_conditions_massive()
        generate_disk_initial_conditions()
        generate_spiral_galaxy_initial_conditions()
        generate_barred_flocculent_spiral_ic()
        integrate_single_particle()
        integrate_single_particle_DEBUG()
        run_mass_integrator_production()
        main_c_function()
        print("    ... All C functions registered.")

        # Step 4.G: Assemble the final C project on disk.
        print(" -> Assembling final C project directory...")
        register_custom_structures_and_params()
        BHaH_defines_h.output_BHaH_defines_h(project_dir=project_dir)

        # The BHaH infrastructure requires some helper files for SIMD intrinsics.
        # This function copies them from the nrpy library into the project.
        import nrpy.helpers.generic as gh
        gh.copy_files(
            package="nrpy.helpers",
            filenames_list=["simd_intrinsics.h"],
            project_dir=project_dir,
            subdirectory="intrinsics",
        )

        # Finally, generate the Makefile.
        Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
            project_dir=project_dir,
            project_name=project_name,
            exec_or_library_name="mass_integrator",
            addl_CFLAGS=["-Wall -Wextra -g $(shell gsl-config --cflags) -fopenmp"],
            addl_libraries=["$(shell gsl-config --libs) -fopenmp"],
        )
        print("    ... Project assembly complete.")

        # Step 4.H: Print a final confirmation message.
        print(f"\nFinished! A C project has been generated in '{project_dir}/'")
        print(f"To build, navigate to this directory in your terminal and type 'make'.")
        print(f"To run, type './mass_integrator'.")

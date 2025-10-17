"""
Generates the complete C project for simulating massive particle geodesics.

This script follows the nrpy "gold standard" format. It sets up all necessary
parameters, executes the symbolic workflow by calling recipes from the shared
geodesic_recipes module, registers all required C functions, and assembles the
final project directory with a Makefile.

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
from typing import Dict, List, Tuple

# Step 0.b: Import third-party modules
import numpy as np
import sympy as sp
import numba
from tqdm import tqdm

# Step 0.c: Import NRPy+ core modules
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.generic as gh
import nrpy.indexedexp as ixp
import nrpy.params as par

# Step 0.d: Import NRPy+ BHaH infrastructure modules
from nrpy.infrastructures.BHaH import (
    BHaH_defines_h,
    CodeParameters as CPs,
    cmdline_input_and_parfiles,
    Makefile_helpers as Makefile,
)

# Step 0.e: Import the shared symbolic recipes module
from nrpy.helpers import geodesic_recipes as georec


# ##############################################################################
# PART 1: SYMBOLIC RECIPES & C-CODE GENERATION (THE BUILDERS)
# ##############################################################################

# --------------------------------------------------------------------------
# Section 1.A: Mass-Integrator-Specific Symbolic Recipes
# --------------------------------------------------------------------------

def ut_massive() -> sp.Expr:
    """
    Symbolically derive u^t for a MASSIVE particle from the 4-velocity.

    The derivation comes from solving the timelike normalization condition:
        g_{μν} u^μ u^ν = -1
    This expands to a quadratic equation for u^t (u^0):
        g_00(u^0)^2 + 2g_0i u^0 u^i + g_ij u^i u^j = -1
    This function solves for u^0 using the quadratic formula.

    Returns:
        A SymPy expression for u^t in terms of metric components and
        spatial 4-velocity components.
    """
    # The symbolic recipe will use the standard variable names u0, u1, etc.
    # The C-generating function will map y[4], y[5], etc. to these.
    u0, u1, u2, u3 = sp.symbols("u0 u1 u2 u3", real=True)
    uU = [u0, u1, u2, u3]

    # The recipe uses the standard name "metric" for the struct.
    g4DD = ixp.declarerank2("metric->g", sym="sym01", dimension=4)

    # This is the quadratic equation for u^0, derived from g_μν u^μ u^ν = -1
    # We solve for u^0.
    # A = g_00
    # B = 2 * g_0i * u^i
    # C = g_ij * u^i * u^j + 1
    sum_g0i_ui = sp.sympify(0)
    for i in range(1, 4):
        sum_g0i_ui += g4DD[0][i] * uU[i]

    sum_gij_ui_uj = sp.sympify(0)
    for i in range(1, 4):
        for j in range(i, 4):
            term = g4DD[i][j] * uU[i] * uU[j]
            if i != j:
                term *= 2
            sum_gij_ui_uj += term

    # The discriminant of the quadratic formula for u^0: B^2 - 4AC
    discriminant = (2 * sum_g0i_ui) ** 2 - 4 * g4DD[0][0] * (sum_gij_ui_uj + 1)

    # We choose the negative root for a forward-in-time particle outside the horizon.
    # u^0 = (-B - sqrt(B^2 - 4AC)) / 2A
    answer = (-2 * sum_g0i_ui - sp.sqrt(discriminant)) / (2 * g4DD[0][0])
    return answer


def symbolic_ut_uphi_from_r_stable(
    M_scale: sp.Symbol, a_spin: sp.Symbol
) -> Tuple[sp.Expr, sp.Expr]:
    """
    Symbolically derive u^t and u^phi for a circular, equatorial orbit.

    This function uses a NUMERICALLY STABLE three-step method to avoid the
    catastrophic cancellation that can occur in the general formula when
    orbits are very close to the black hole.

    The method is as follows:
    1. Calculate Angular Velocity Ω = dφ/dt for a stable prograde circular orbit:
       Ω = sqrt(M) / (r^(3/2) + a*sqrt(M))
    2. Solve for u^t using the 4-velocity normalization condition, which gives:
       (u^t)^2 * (g_tt + 2*g_tφ*Ω + g_φφ*Ω^2) = -1
    3. Calculate u^φ = Ω * u^t.

    Args:
        M_scale: The symbolic mass of the black hole.
        a_spin: The symbolic dimensionless spin parameter of the black hole.

    Returns:
        A tuple containing the symbolic expressions for [u^t, u^phi].
    """
    # Define symbolic variables with final C names.
    r = sp.symbols("r_initial", real=True)

    # Step 1: Calculate the stable angular velocity Omega = d(phi)/dt.
    Omega = sp.sqrt(M_scale) / (r ** sp.Rational(3, 2) + a_spin * sp.sqrt(M_scale))

    # Step 2: Define the required Kerr metric components in Boyer-Lindquist
    # coordinates, specialized to the equatorial plane (theta = pi/2).
    g_tt = -(1 - 2 * M_scale / r)
    g_tphi = -2 * a_spin * M_scale / r
    g_phiphi = r**2 + a_spin**2 + (2 * M_scale * a_spin**2) / r

    # Step 3: Solve for u^t using the 4-velocity normalization condition.
    ut_squared_inv_denom = g_tt + 2 * g_tphi * Omega + g_phiphi * Omega**2
    ut_squared = -1 / ut_squared_inv_denom
    ut = sp.sqrt(ut_squared)

    # Step 4: Calculate u^phi.
    uphi = Omega * ut

    return ut, uphi


def symbolic_conserved_quantities() -> Tuple[sp.Expr, List[sp.Expr], sp.Expr]:
    """
    Generate symbolic recipes for conserved quantities E, L_i, and Q.

    - Energy (E): Conserved due to time-independence (stationarity).
      E = -p_t = -g_{tμ} p^μ
    - Angular Momentum (L_i): The z-component is conserved due to axisymmetry.
      L_i = ε_{ijk} x^j p_k
    - Carter Constant (Q): Conserved due to a hidden symmetry of the Kerr metric.
      For massive particles (m=1):
      Q = p_θ^2 + cos^2(θ) * (a^2(1 - E^2) + L_z^2/sin^2(θ))

    Returns:
        A tuple (E_expr, L_exprs, Q_expr), where E_expr is the energy,
        L_exprs is a list [Lx, Ly, Lz], and Q_expr is the Carter Constant.
    """
    # Define coordinate and 4-momentum components using the y[i] convention.
    t, x, y, z = sp.symbols("y[0] y[1] y[2] y[3]", real=True)
    pt, px, py, pz = sp.symbols("y[4] y[5] y[6] y[7]", real=True)
    pU = [pt, px, py, pz]
    a = sp.Symbol("a_spin", real=True)

    # Define an abstract metric tensor to be filled by a C struct at runtime.
    g4DD = ixp.declarerank2("metric->g", sym="sym01", dimension=4)

    # --- 1. Calculate Energy E = -p_t ---
    p_t = sp.sympify(0)
    for mu in range(4):
        p_t += g4DD[0][mu] * pU[mu]
    E_expr = -p_t

    # --- 2. Calculate Angular Momentum L_i ---
    # First, find covariant momentum components p_k = g_{kμ} p^μ.
    p_down = ixp.zerorank1(dimension=4)
    for k in range(1, 4):  # We only need p_x, p_y, p_z for L_i
        for mu in range(4):
            p_down[k] += g4DD[k][mu] * pU[mu]
    p_x, p_y, p_z = p_down[1], p_down[2], p_down[3]

    # Then, calculate L_i = ε_{ijk} x^j p_k.
    L_x = y * p_z - z * p_y
    L_y = z * p_x - x * p_z
    L_z = x * p_y - y * p_x
    L_exprs = [L_x, L_y, L_z]

    # --- 3. Calculate Carter Constant Q ---
    # First, compute geometric terms in Cartesian coordinates.
    r_sq = x**2 + y**2 + z**2
    rho_sq = x**2 + y**2

    # Then, compute p_theta^2 directly from Cartesian components.
    xpx_plus_ypy = x * p_x + y * p_y
    p_theta_sq = (
        (z**2 * xpx_plus_ypy**2 / rho_sq)
        - (2 * z * p_z * xpx_plus_ypy)
        + (rho_sq * p_z**2)
    )

    # Assemble the final formula for Q.
    # Q = p_θ^2 + cos^2(θ) * (a^2(1 - E^2) + L_z^2/sin^2(θ))
    # where cos^2(θ) = z^2/r^2 and sin^2(θ) = rho^2/r^2
    second_term = (z**2 / r_sq) * (a**2 * (1 - E_expr**2) + L_z**2 * (r_sq / rho_sq))
    Q_formula = p_theta_sq + second_term

    # Handle the axial singularity (rho_sq -> 0), where Q should be 0.
    Q_expr = sp.Piecewise((0, rho_sq < 1e-12), (Q_formula, True))

    return E_expr, L_exprs, Q_expr


def geodesic_eom_rhs() -> Tuple[List[sp.Expr], List[sp.Expr]]:
    """
    Generate symbolic RHS expressions for the massive geodesic equations of motion.

    This function defines the system of 8 first-order ODEs:
    1. Position ODEs: dx^α/dτ = u^α
    2. Velocity ODEs: du^α/dτ = -Γ^α_{μν} u^μ u^ν

    It constructs the velocity RHS using a symmetry optimization to reduce the
    number of terms from 16 to 10 for each component.

    Returns:
        A tuple (pos_rhs, vel_rhs), where pos_rhs is a list of the 4 position
        RHS expressions and vel_rhs is a list of the 4 velocity RHS expressions.
    """
    # --- Position ODE: dx^α/dτ = u^α ---
    ut, ux, uy, uz = sp.symbols("y[4] y[5] y[6] y[7]", real=True)
    pos_rhs = [ut, ux, uy, uz]

    # --- Velocity ODE: du^α/dτ = -Γ^α_{μν} u^μ u^ν ---
    # Use a placeholder for the Christoffel symbols, which will be provided by a C struct.
    Gamma4UDD = ixp.declarerank3("conn->Gamma4UDD", dimension=4, sym="sym12")
    uU = [ut, ux, uy, uz]
    vel_rhs = ixp.zerorank1(dimension=4)
    for alpha in range(4):
        sum_term = sp.sympify(0)
        # Exploit symmetry Γ^α_{μν} = Γ^α_{νμ} to optimize the sum.
        for mu in range(4):
            for nu in range(mu, 4):
                term = Gamma4UDD[alpha][mu][nu] * uU[mu] * uU[nu]
                if mu != nu:
                    term *= 2
                sum_term += term
        vel_rhs[alpha] = -sum_term

    return pos_rhs, vel_rhs

# mass_geodesic_integrator.py (Continued)

# ##############################################################################
# PART 2: C CODE GENERATION (THE BUILDERS)
# ##############################################################################
#
# This section contains the Python functions that generate and register all the
# necessary C code for the project. Each function is a "builder" that takes
# symbolic expressions (the "recipes") and uses nrpy's C code generation
# tools to produce an optimized C function, which is then registered with the
# CFunction_dict for final project assembly.
#

# --------------------------------------------------------------------------
# Section 2.A: Tier 4 - Low-Level Workers & Dispatchers
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
    # Create a flattened list of the 10 unique symbolic metric components.
    list_of_g4DD_syms = [g4DD_expressions[i][j] for i in range(4) for j in range(i, 4)]

    # Create a corresponding list of C variable names for the output.
    list_of_g4DD_C_vars = [f"metric->g{i}{j}" for i in range(4) for j in range(i, 4)]

    # Define local variables for each argument of cfc.register_CFunction.
    includes = ["BHaH_defines.h"]
    desc = "Computes the 10 unique components of the Kerr metric in Cartesian Kerr-Schild coords."
    name = "g4DD_kerr_schild"
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const double y[4],
                  metric_struct *restrict metric"""
    body = ccg.c_codegen(
        list_of_g4DD_syms, list_of_g4DD_C_vars, enable_cse=True, cse_varprefix="KS"
    )

    # Register the C function with the nrpy framework.
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=True,
    )


def con_kerr_schild(Gamma4UDD_expressions: List[List[List[sp.Expr]]]) -> None:
    """
    Generate and register the C worker for the Kerr-Schild Christoffel symbols.

    This function takes the symbolic expressions for the 40 unique components
    of the Christoffel symbols and generates an optimized C function to compute them.

    Args:
        Gamma4UDD_expressions: A rank-3 list of lists of lists containing the
                               symbolic expressions for the Christoffel symbols.
    """
    # Create a flattened list of the 40 unique symbolic Christoffel components.
    list_of_Gamma_syms = [
        Gamma4UDD_expressions[i][j][k]
        for i in range(4)
        for j in range(4)
        for k in range(j, 4)
    ]

    # Create a corresponding list of C variable names for the output.
    conn_Gamma4UDD = ixp.declarerank3("conn->Gamma4UDD", dimension=4)
    list_of_Gamma_C_vars = [
        str(conn_Gamma4UDD[i][j][k])
        for i in range(4)
        for j in range(4)
        for k in range(j, 4)
    ]

    # Define local variables for each argument of cfc.register_CFunction.
    includes = ["BHaH_defines.h"]
    desc = "Computes the 40 unique Christoffel symbols for the Kerr metric in Kerr-Schild coords."
    name = "con_kerr_schild"
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const double y[4],
                  connection_struct *restrict conn"""
    body = ccg.c_codegen(
        list_of_Gamma_syms, list_of_Gamma_C_vars, enable_cse=True, cse_varprefix="conn"
    )

    # Register the C function.
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=True,
    )


def g4dd_metric() -> None:
    """
    Generate and register the C dispatcher for metric calculations.

    This function generates a C function `g4DD_metric()` that contains a switch
    statement. Based on the runtime `metric->type`, it calls the appropriate
    specialized worker function (e.g., `g4DD_kerr_schild`). This decouples the
    high-level logic from the specific metric implementation.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Dispatcher to compute the 4-metric g_munu for the chosen analytic metric."
    name = "g4DD_metric"
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric,
                  const double y[8],
                  metric_struct *restrict metric_out"""
    body = r"""
    // The state vector y_pos contains only the position coordinates.
    const double y_pos[4] = {y[0], y[1], y[2], y[3]};

    // This switch statement chooses which "worker" function to call
    // based on the metric type provided.
    switch(metric->type) {
        case Schwarzschild:
        case Kerr:
            // For Kerr or Schwarzschild in KS coords, call the unified Kerr-Schild C function.
            g4DD_kerr_schild(commondata, params, y_pos, metric_out);
            break;
        case Schwarzschild_Standard:
            // In this project, we use the Kerr-Schild implementation for a=0
            // instead of a separate standard Schwarzschild implementation.
            g4DD_kerr_schild(commondata, params, y_pos, metric_out);
            break;
        default:
            fprintf(stderr, "Error: MetricType %d not supported in g4DD_metric().\\n", metric->type);
            exit(1);
    }
    """
    cfc.register_CFunction(
        includes=includes, desc=desc, name=name, params=params, body=body
    )


def connections() -> None:
    """
    Generate and register the C dispatcher for Christoffel symbol calculations.

    Similar to g4DD_metric(), this function generates a C dispatcher that calls
    the appropriate worker function (e.g., `con_kerr_schild`) based on the
    runtime `metric->type`.
    """
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "stdio.h",
        "stdlib.h",
    ]
    desc = "Dispatcher to compute Christoffel symbols for the chosen analytic metric."
    name = "connections"
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric,
                  const double y[8],
                  connection_struct *restrict conn"""
    body = r"""
    // The state vector y_pos contains only the position coordinates.
    const double y_pos[4] = {y[0], y[1], y[2], y[3]};

    // This switch statement chooses which "worker" function to call
    // based on the metric type provided.
    switch(metric->type) {
        case Schwarzschild:
        case Kerr:
            con_kerr_schild(commondata, params, y_pos, conn);
            break;
        case Schwarzschild_Standard:
            con_kerr_schild(commondata, params, y_pos, conn);
            break;
        default:
            fprintf(stderr, "Error: MetricType %d not supported in connections().\\n", metric->type);
            exit(1);
    }
    """
    cfc.register_CFunction(
        includes=includes, desc=desc, name=name, params=params, body=body
    )


# --------------------------------------------------------------------------
# Section 2.B: Tier 3 - Core Subsystems & Engines
# --------------------------------------------------------------------------


def calculate_ode_rhs_massive(
    pos_rhs_expressions: List[sp.Expr], vel_rhs_expressions: List[sp.Expr]
) -> None:
    """
    Generate the C engine for the massive particle geodesic ODE right-hand sides.

    This function is a generic engine that implements the geodesic equations
    using pre-computed Christoffel symbols from the `connection_struct`.

    Args:
        pos_rhs_expressions: List of 4 symbolic expressions for dx^α/dτ.
        vel_rhs_expressions: List of 4 symbolic expressions for du^α/dτ.
    """
    all_rhs_expressions = pos_rhs_expressions + vel_rhs_expressions
    rhs_output_vars = [f"rhs_out[{i}]" for i in range(8)]

    includes = ["BHaH_defines.h"]
    desc = "Calculates the RHS of the 8 massive particle geodesic ODEs."
    name = "calculate_ode_rhs_massive"
    params = """const double y[8],
                  const connection_struct *restrict conn,
                  double rhs_out[8]"""
    body = ccg.c_codegen(all_rhs_expressions, rhs_output_vars, enable_cse=True)

    cfc.register_CFunction(
        includes=includes, desc=desc, name=name, params=params, body=body
    )


def check_conservation_massive(
    E_expr: sp.Expr, L_exprs: List[sp.Expr], Q_expr: sp.Expr
) -> None:
    """
    Generate the C engine for computing conserved quantities (E, L_i, Q).

    This engine is used for validation. It computes the conserved quantities
    from a particle's state vector, allowing for the measurement of numerical
    accuracy by checking their constancy over time.

    Args:
        E_expr: Symbolic expression for energy.
        L_exprs: List of symbolic expressions for [Lx, Ly, Lz].
        Q_expr: Symbolic expression for the Carter Constant.
    """
    list_of_expressions = [E_expr] + L_exprs + [Q_expr]
    output_vars = ["*E", "*Lx", "*Ly", "*Lz", "*Q"]
    body_C_code = ccg.c_codegen(
        list_of_expressions, output_vars, enable_cse=True, include_braces=False
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Computes conserved quantities (E, L_i, Q) for a given massive particle state vector."
    name = "check_conservation_massive"
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric_params_in,
                  const double y[8],
                  double *E, double *Lx, double *Ly, double *Lz, double *Q"""
    body = f"""
    // Unpack parameters from commondata struct that are needed symbolically.
    const REAL a_spin = commondata->a_spin;

    // Allocate a metric struct on the heap to prevent stack overflow.
    metric_struct* metric = (metric_struct*)malloc(sizeof(metric_struct));
    if (metric == NULL) {{
        fprintf(stderr, "Error: Failed to allocate memory for metric_struct in check_conservation.\\n");
        exit(1);
    }}
    
    // Call the dispatcher to fill the struct with metric components.
    g4DD_metric(commondata, params, metric_params_in, y, metric);

    // Execute the NRPy-generated code for conserved quantities.
    {body_C_code}
    
    // Free the heap-allocated memory.
    free(metric);
    """
    cfc.register_CFunction(
        includes=includes, desc=desc, name=name, params=params, body=body
    )


def ode_gsl_wrapper_massive() -> None:
    """
    Generate the C wrapper function to interface with the GSL ODE solver.

    This function generates a C function with the exact signature required by
    GSL's ODE driver. It acts as a bridge, unpacking a generic `void*` pointer
    to access simulation data, calling our project-specific dispatchers and
    engines, and returning the computed derivatives to the GSL solver.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "gsl/gsl_errno.h"]
    desc = "GSL wrapper for the massive particle geodesic ODEs."
    cfunc_type = "int"
    name = "ode_gsl_wrapper_massive"
    params = "double t, const double y[8], double f[8], void *params"
    body = r"""
    // Proper time 't' (the independent variable for GSL) is not explicitly
    // used in the RHS expressions, so we mark it as unused to avoid warnings.
    (void)t;

    // Unpack the carrier struct to access simulation parameters and metric choice.
    gsl_params *gsl_parameters = (gsl_params *)params;

    // Declare structs on the stack to hold metric and connection values.
    metric_struct g4DD;
    connection_struct conn;

    // Call dispatchers to compute the metric and Christoffel symbols at the current position y.
    // Note: The y array is 8D, but these functions only need the first 4 position components.
    g4DD_metric(gsl_parameters->commondata, gsl_parameters->params, gsl_parameters->metric, y, &g4DD);
    connections(gsl_parameters->commondata, gsl_parameters->params, gsl_parameters->metric, y, &conn);

    // Call the engine to compute the RHS of the ODEs and store the result in f[].
    calculate_ode_rhs_massive(y, &conn, f);

    // Signal success to the GSL solver.
    return GSL_SUCCESS;
    """
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )

# mass_geodesic_integrator.py (Continued)

# ##############################################################################
# PART 3: C CODE GENERATION (THE BUILDERS, CONTINUED)
# ##############################################################################

# --------------------------------------------------------------------------
# Section 3.A: Tier 2 - Initial Condition & Integration Orchestrators
# --------------------------------------------------------------------------

def calculate_ut_uphi_from_r(ut_expr: sp.Expr, uphi_expr: sp.Expr) -> None:
    """
    Generate the C helper for computing u^t and u^phi from a radius.

    This function implements the numerically stable recipe for finding the
    initial 4-velocity components for a stable, circular, equatorial orbit.

    Args:
        ut_expr: The symbolic expression for u^t.
        uphi_expr: The symbolic expression for u^phi.
    """
    includes = ["BHaH_defines.h"]
    desc = "Computes u^t and u^phi for a circular orbit at a given radius using a numerically stable method."
    name = "calculate_ut_uphi_from_r"
    params = """const double r_initial,
                  const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  double *ut, double *uphi"""
    body = ccg.c_codegen(
        [ut_expr, uphi_expr], ["*ut", "*uphi"], enable_cse=True, cse_varprefix="ic"
    )

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=True,
    )


def set_initial_conditions_massive() -> None:
    """
    Generate the C engine to set the full initial 8-component state vector.

    This orchestrator takes a particle's initial position and computes the full
    state vector for a stable, circular, equatorial orbit. It calls the
    `calculate_ut_uphi_from_r` helper and then transforms the resulting u^phi
    into Cartesian 4-velocity components u^x and u^y.
    """
    includes = ["BHaH_defines.h", "<math.h>"]
    desc = "Sets the initial 8-component state vector for a massive particle using a numerically stable method."
    name = "set_initial_conditions_massive"
    params = """const particle_initial_state_t *restrict initial_state,
                  const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  double y_out[8]"""
    body = r"""
    // Unpack parameters for clarity
    const double M = commondata->M_scale;
    const double a = commondata->a_spin;

    // Unpack initial position
    const double x = initial_state->pos[1];
    const double y = initial_state->pos[2];
    const double z = initial_state->pos[3];
    const double r = sqrt(x*x + y*y + z*z);

    // --- Numerically Stable Method to find u^t and u^phi ---
    double ut, uphi;
    calculate_ut_uphi_from_r(r, commondata, params, &ut, &uphi);

    // Check for unstable or invalid orbits where ut becomes NaN.
    if (isnan(ut)) {
        for(int i=0; i<8; i++) y_out[i] = NAN;
        return;
    }

    // --- Assemble the final 8-component state vector ---
    y_out[0] = initial_state->pos[0]; // t
    y_out[1] = x;
    y_out[2] = y;
    y_out[3] = z;
    y_out[4] = ut;
    // Transform u^phi to Cartesian u^x, u^y:
    // u^x = dx/dτ = d(r*cosφ)/dτ = -r*sin(φ)*(dφ/dτ) = -y*u^φ
    // u^y = dy/dτ = d(r*sinφ)/dτ =  r*cos(φ)*(dφ/dτ) =  x*u^φ
    y_out[5] = -y * uphi;
    y_out[6] =  x * uphi;
    y_out[7] = 0.0; // u^z is zero for equatorial orbits
    """
    cfc.register_CFunction(
        includes=includes, desc=desc, name=name, params=params, body=body
    )


def integrate_single_particle() -> None:
    """
    Generate the C orchestrator for integrating a single particle in production.

    This function is the workhorse of the main simulation loop. It evolves a
    particle from a start time `t_start` to an end time `t_end` using the
    high-level GSL driver, which ensures the particle state is returned at
    the exact `t_end`. It also includes termination checks.
    """
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "gsl/gsl_errno.h",
        "gsl/gsl_odeiv2.h",
        "<math.h>",
    ]
    desc = "Integrates a single massive particle path between two coordinate times using the GSL driver."
    cfunc_type = "int"
    name = "integrate_single_particle"
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric,
                  const double t_start, const double t_end,
                  double y_in_out[8]"""
    body = r"""
    // Define the GSL ODE system, pointing to our wrapper function.
    gsl_params gsl_parameters = {commondata, params, metric};
    gsl_odeiv2_system sys = {ode_gsl_wrapper_massive, NULL, 8, &gsl_parameters};

    // Set up the GSL driver with an adaptive RKF45 stepper.
    // Error tolerances are 1e-11 for both absolute and relative error.
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-11, 1e-11);

    double t = t_start;

    // The driver will take internal steps to reach t_end precisely.
    int status = gsl_odeiv2_driver_apply(d, &t, t_end, y_in_out);

    // Free memory before checking status.
    gsl_odeiv2_driver_free(d);

    if (status != GSL_SUCCESS) {
        return 1; // Return failure code
    }

    // Robustness check after the step for termination conditions.
    const double r_sq = y_in_out[1]*y_in_out[1] + y_in_out[2]*y_in_out[2] + y_in_out[3]*y_in_out[3];
    const double r_horizon = commondata->M_scale * (1.0 + sqrt(1.0 - commondata->a_spin*commondata->a_spin));

    if (r_sq < r_horizon*r_horizon || r_sq > commondata->r_escape*commondata->r_escape || fabs(y_in_out[4]) > commondata->ut_max) {
        return 1; // Return failure code
    }

    return 0; // Return success code
    """
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=True,
    )


def integrate_single_particle_DEBUG() -> None:
    """
    Generate the C orchestrator for a single-particle debug run.

    This function generates a C function that traces the full, detailed path
    of a single particle, writing its state at each internal step of the
    integrator to a text file for validation and visualization. It uses the
    low-level GSL `evolve` interface for fine-grained control.
    """
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "gsl/gsl_errno.h",
        "gsl/gsl_odeiv2.h",
    ]
    desc = "DEBUG integrator for a single massive particle. Writes the full trajectory to a text file."
    name = "integrate_single_particle_DEBUG"
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric,
                  const double start_y[8],
                  double final_y_state[8]"""
    body = r"""
    // GSL Setup using the low-level interface for step-by-step evolution.
    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step * step = gsl_odeiv2_step_alloc(T, 8);
    gsl_odeiv2_control * control = gsl_odeiv2_control_yp_new(1e-14, 1e-14);
    gsl_odeiv2_evolve * evol = gsl_odeiv2_evolve_alloc(8);
    gsl_params gsl_parameters = {commondata, params, metric};
    gsl_odeiv2_system sys = {ode_gsl_wrapper_massive, NULL, 8, &gsl_parameters};

    double y_c[8];
    double t = 0.0, dt = 0.01; // t is proper time τ
    for (int j = 0; j < 8; j++) { y_c[j] = start_y[j]; }

    // Setup output file.
    FILE *fp = fopen("massive_particle_path.txt", "w");
    if (fp == NULL) {
        fprintf(stderr, "Error: Could not open massive_particle_path.txt for writing.\\n");
        exit(1);
    }
    fprintf(fp, "# ProperTime_tau\tCoordTime_t\tx\ty\tz\tu^t\tu^x\tu^y\tu^z\n");

    printf("Starting debug trace for single massive particle...\n");
    printf("Step | Proper Time (τ) | Coord Time (t) |      x     |      y     |      z     |      u^t   \n");
    printf("-------------------------------------------------------------------------------------------\n");

    // Main Integration Loop.
    for (int i = 0; i < 2000000; i++) {
        int status = gsl_odeiv2_evolve_apply(evol, control, step, &sys, &t, 1e10, &dt, y_c);

        // Write full state to file.
        fprintf(fp, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                t, y_c[0], y_c[1], y_c[2], y_c[3], y_c[4], y_c[5], y_c[6], y_c[7]);

        if (i % 500 == 0) {
            printf("%4d | %15.4e | %14.4f | %10.4f | %10.4f | %10.4f | %10.4f\n",
                   i, t, y_c[0], y_c[1], y_c[2], y_c[3], y_c[4]);
        }

        const double r_sq = y_c[1]*y_c[1] + y_c[2]*y_c[2] + y_c[3]*y_c[3];
        const double r_horizon = commondata->M_scale * (1.0 + sqrt(1.0 - commondata->a_spin*commondata->a_spin));

        // Termination Conditions.
        if (status != GSL_SUCCESS) { printf("Termination: GSL ERROR (status = %d)\\n", status); break; }
        if (r_sq < r_horizon*r_horizon) { printf("Termination: Fell below event horizon (r=%.2f)\\n", sqrt(r_sq)); break; }
        if (r_sq > commondata->r_escape*commondata->r_escape) { printf("Termination: Escaped to r > %.1f\\n", commondata->r_escape); break; }
        if (fabs(y_c[4]) > commondata->ut_max) { printf("Termination: Runaway u^t > %.1e\\n", commondata->ut_max); break; }
        if (fabs(y_c[0]) > commondata->t_max_integration) { printf("Termination: Exceeded max integration time t > %.1f\\n", commondata->t_max_integration); break; }
    }

    // Copy final state to output and clean up.
    for(int j=0; j<8; j++) { final_y_state[j] = y_c[j]; }
    fclose(fp);
    gsl_odeiv2_evolve_free(evol);
    gsl_odeiv2_control_free(control);
    gsl_odeiv2_step_free(step);
    """
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=True,
    )

# mass_geodesic_integrator.py (Continued)

# ##############################################################################
# PART 4: C CODE GENERATION (BUILDERS CONCLUDED) & POST-PROCESSING LOGIC
# ##############################################################################

# --------------------------------------------------------------------------
# Section 4.A: Tier 1 & 2 - High-Level Orchestrators & main()
# --------------------------------------------------------------------------

def generate_disk_initial_conditions() -> None:
    """
    Generate a C orchestrator to create a simple, axisymmetric Keplerian disk.

    This C function populates a large array of initial state vectors by looping
    through a grid of radii `r` and azimuthal angles `phi`. At each point, it
    calls the `calculate_ut_uphi_from_r` helper to get the correct 4-velocity
    for a stable circular orbit and assembles the full 8-component state vector.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Generates the complete y[8] initial state for all particles in a Keplerian disk."
    cfunc_type = "int"
    name = "generate_disk_initial_conditions"
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  double *restrict y_initial_states"""
    body = r"""
    int particle_count = 0;
    const double dr = (commondata->disk_num_r > 1) ? (commondata->disk_r_max - commondata->disk_r_min) / (commondata->disk_num_r - 1) : 0;

    for (int i = 0; i < commondata->disk_num_r; i++) {
        const double r = commondata->disk_r_min + i * dr;

        // Adjust number of azimuthal particles based on radius to maintain density.
        const int num_phi_at_r = (commondata->disk_num_phi > 1) ? (int)(commondata->disk_num_phi * (r / commondata->disk_r_max)) : 1;
        if (num_phi_at_r == 0) continue;
        const double dphi = 2.0 * M_PI / num_phi_at_r;

        double ut_at_r, uphi_at_r;
        calculate_ut_uphi_from_r(r, commondata, params, &ut_at_r, &uphi_at_r);

        for (int j = 0; j < num_phi_at_r; j++) {
            const double phi = j * dphi;
            const double cos_phi = cos(phi);
            const double sin_phi = sin(phi);

            // Set a pointer to the start of the state vector for this particle.
            double *y = &y_initial_states[particle_count * 8];

            y[0] = 0.0; // t
            y[1] = r * cos_phi; // x
            y[2] = r * sin_phi; // y
            y[3] = 0.0; // z

            y[4] = ut_at_r;
            y[5] = -(r * sin_phi) * uphi_at_r; // u^x
            y[6] =  (r * cos_phi) * uphi_at_r; // u^y
            y[7] = 0.0; // u^z

            particle_count++;
        }
    }
    return particle_count;
    """
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )


def generate_spiral_galaxy_initial_conditions() -> None:
    """
    Generate a C orchestrator to create a "grand design" spiral galaxy.

    This function uses a random placement strategy. For each particle, it
    chooses a radius and then calculates an angle based on the logarithmic
    spiral formula, adding some scatter for thickness. It then calls
    `calculate_ut_uphi_from_r` to place the particle on a stable circular orbit.
    """
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<stdlib.h>",
    ]
    desc = "Generates the complete y[8] initial state for all particles in a spiral galaxy disk."
    cfunc_type = "int"
    name = "generate_spiral_galaxy_initial_conditions"
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  double *restrict y_initial_states"""
    body = r"""
    int particle_count = 0;
    const int num_particles_total = commondata->disk_num_r * commondata->disk_num_phi;

    const int num_arms = commondata->spiral_galaxy_num_arms;
    const double arm_tightness = commondata->spiral_galaxy_arm_tightness;

    // Seed the random number generator for reproducibility.
    srand(42);

    for (int i = 0; i < num_particles_total; i++) {
        // Randomly choose a radius, weighted to be uniform in area.
        const double r = commondata->disk_r_min + (commondata->disk_r_max - commondata->disk_r_min) * sqrt((double)rand() / RAND_MAX);

        // Calculate the base angle for this radius from the logarithmic spiral formula.
        const double theta_base = (1.0 / arm_tightness) * log(r / commondata->disk_r_min);

        // Randomly assign the particle to one of the spiral arms.
        const int arm_index = rand() % num_arms;
        const double arm_offset = (2.0 * M_PI / num_arms) * arm_index;

        // Add scatter to give the arms some thickness.
        const double phi_spread = (M_PI / num_arms) * 0.2 * ((double)rand() / RAND_MAX - 0.5);
        const double phi = theta_base + arm_offset + phi_spread;

        const double cos_phi = cos(phi);
        const double sin_phi = sin(phi);

        double ut_at_r, uphi_at_r;
        calculate_ut_uphi_from_r(r, commondata, params, &ut_at_r, &uphi_at_r);

        double *y = &y_initial_states[particle_count * 8];

        y[0] = 0.0;
        y[1] = r * cos_phi;
        y[2] = r * sin_phi;
        y[3] = 0.0;
        y[4] = ut_at_r;
        y[5] = -(r * sin_phi) * uphi_at_r;
        y[6] =  (r * cos_phi) * uphi_at_r;
        y[7] = 0.0;

        particle_count++;
    }
    return particle_count;
    """
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )


def generate_barred_flocculent_spiral_ic() -> None:
    """
    Generate the C orchestrator for a barred, clumpy ("flocculent") spiral.

    This function uses rejection sampling to create a complex geometry with a
    dense central bulge, a rectangular bar, and clumpy spiral arms. Despite
    the complex placement, each particle is still placed on a stable circular
    orbit by calling `calculate_ut_uphi_from_r`.
    """
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<stdlib.h>",
        "<time.h>",
    ]
    desc = "Generates the initial state for all particles in a barred flocculent spiral galaxy."
    cfunc_type = "int"
    name = "generate_barred_flocculent_spiral_ic"
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  double *restrict y_initial_states"""
    body = r"""
    // Seed the random number generator.
    unsigned int seed = 42;

    int particle_count = 0;
    const int num_particles_total = commondata->disk_num_r * commondata->disk_num_phi;

    // Unpack geometry parameters from commondata for clarity.
    const double r_min = commondata->disk_r_min;
    const double r_max = commondata->disk_r_max;
    const double bar_len = commondata->bar_length;
    const double bar_width = commondata->bar_length * commondata->bar_aspect_ratio;
    const double bulge_rad = commondata->bulge_radius;
    const int num_arms = commondata->spiral_galaxy_num_arms;
    const double arm_tightness = commondata->spiral_galaxy_arm_tightness;

    for (int i = 0; i < num_particles_total; i++) {
        double r, phi, x, y;
        seed += i; // Use a unique seed for each particle.

        // --- REJECTION SAMPLING LOOP ---
        while (1) {
            // Generate a trial particle position uniformly in the disk annulus.
            double random_val = (double)rand_r(&seed) / RAND_MAX;
            r = sqrt(random_val * (r_max*r_max - r_min*r_min) + r_min*r_min);
            phi = 2.0 * M_PI * ((double)rand_r(&seed) / RAND_MAX);
            x = r * cos(phi);
            y = r * sin(phi);

            double acceptance_prob = 0.0;

            // Region 1: Central Bulge
            if (r < bulge_rad) {
                acceptance_prob = commondata->bulge_density_factor * commondata->arm_particle_density;
            }
            // Region 2: Central Bar
            else if (fabs(x) < bar_len / 2.0 && fabs(y) < bar_width / 2.0) {
                acceptance_prob = commondata->bar_density_factor * commondata->arm_particle_density;
            }
            // Region 3: Flocculent Spiral Arms
            else {
                if (((double)rand_r(&seed) / RAND_MAX) < commondata->arm_particle_density) {
                    double theta_base = (1.0 / arm_tightness) * log(r / r_min);
                    double delta_phi_raw = phi - theta_base;
                    double angle_between_arms = 2.0 * M_PI / num_arms;
                    double delta_phi_folded = fmod(delta_phi_raw, angle_between_arms);

                    if (delta_phi_folded > 0.5 * angle_between_arms) {
                        delta_phi_folded -= angle_between_arms;
                    } else if (delta_phi_folded < -0.5 * angle_between_arms) {
                        delta_phi_folded += angle_between_arms;
                    }

                    double arm_profile = exp(-0.5 * SQR(delta_phi_folded / (arm_tightness * 0.5)));
                    double clump_modulation = cos(commondata->arm_clumpiness_factor * (theta_base - phi));
                    clump_modulation *= cos(num_arms * phi * commondata->arm_clump_size);
                    acceptance_prob = arm_profile + 0.5 * clump_modulation;
                }
            }

            if (((double)rand_r(&seed) / RAND_MAX) < acceptance_prob) {
                break; // Accept this particle's position.
            }
        } // End of rejection sampling while loop.

        // --- Velocity Assignment ---
        double ut_at_r, uphi_at_r;
        calculate_ut_uphi_from_r(r, commondata, params, &ut_at_r, &uphi_at_r);

        // --- Store the Full 8-Component State Vector ---
        double *y_state = &y_initial_states[particle_count * 8];
        y_state[0] = 0.0;
        y_state[1] = x;
        y_state[2] = y;
        y_state[3] = 0.0;
        y_state[4] = ut_at_r;
        y_state[5] = -y * uphi_at_r;
        y_state[6] =  x * uphi_at_r;
        y_state[7] = 0.0;

        particle_count++;
    }
    return particle_count;
    """
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )


def run_mass_integrator_production() -> None:
    """
    Generate the top-level orchestrator for a full production run.

    This C function sets up the particle ensemble by calling the appropriate
    initial condition generator, then enters the main time loop. At each step,
    it saves a snapshot and then calls `integrate_single_particle` for each
    active particle in a parallel OpenMP loop.
    """
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<sys/stat.h>",
        "<string.h>",
        "<stdio.h>",
        "<stdlib.h>",
    ]
    desc = "Orchestrates the full production run for the mass integrator."
    name = "run_mass_integrator_production"
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric"""
    body = r"""
    // Step 1: Allocate memory for the flat initial state array.
    const int max_particles = commondata->disk_num_r * commondata->disk_num_phi;
    double *y_initial_states = (double *)malloc(sizeof(double) * 8 * max_particles);
    if (y_initial_states == NULL) { fprintf(stderr, "Error: Failed to allocate memory for initial states.\\n"); exit(1); }
    int num_particles = 0;

    // Step 2: Dispatch to the correct initial conditions generator.
    printf("Generating initial conditions of type: %s\\n", commondata->initial_conditions_type);
    if (strcmp(commondata->initial_conditions_type, "SpiralGalaxy") == 0) {
        num_particles = generate_spiral_galaxy_initial_conditions(commondata, params, y_initial_states);
    } else if (strcmp(commondata->initial_conditions_type, "BarredFlocculentSpiral") == 0) {
        num_particles = generate_barred_flocculent_spiral_ic(commondata, params, y_initial_states);
    } else { // Default to KeplerianDisk
        if (strcmp(commondata->initial_conditions_type, "KeplerianDisk") != 0) {
            printf("Warning: Unrecognized initial_conditions_type '%s'. Defaulting to KeplerianDisk.\\n", commondata->initial_conditions_type);
        }
        num_particles = generate_disk_initial_conditions(commondata, params, y_initial_states);
    }
    printf("Generated %d particles.\\n", num_particles);

    // Step 3: Create and populate our primary data structure, an array of mass_particle_state_t.
    mass_particle_state_t *particle_states = (mass_particle_state_t *)malloc(sizeof(mass_particle_state_t) * num_particles);
    if (particle_states == NULL) { fprintf(stderr, "Error: Failed to allocate memory for particle states.\\n"); exit(1); }
    for (int i=0; i<num_particles; i++) {
        double *y_start = &y_initial_states[i*8];
        particle_states[i].id = i;
        particle_states[i].pos[0] = y_start[1];
        particle_states[i].pos[1] = y_start[2];
        particle_states[i].pos[2] = y_start[3];
        particle_states[i].u[0] = y_start[4];
        particle_states[i].u[1] = y_start[5];
        particle_states[i].u[2] = y_start[6];
        particle_states[i].u[3] = y_start[7];
        const double r = sqrt(y_start[1]*y_start[1] + y_start[2]*y_start[2]);
        particle_states[i].lambda_rest = commondata->disk_lambda_rest_at_r_min * pow(r / commondata->disk_r_min, 0.75);
        particle_states[i].j_intrinsic = (float)pow(r / commondata->disk_r_min, -3.0);
    }
    free(y_initial_states);

    // Step 4: Create the output directory.
    mkdir(commondata->output_folder, 0755);

    // Step 5: Main Time Evolution and Snapshotting Loop.
    int snapshot_count = 0;
    for (double current_t = 0; current_t <= commondata->t_final; current_t += commondata->snapshot_every_t) {
        char filename[200];
        snprintf(filename, 200, "%s/mass_blueprint_t_%04d.bin", commondata->output_folder, snapshot_count);
        printf("Saving snapshot: %s (t=%.2f)\\n", filename, current_t);

        FILE *fp_out = fopen(filename, "wb");
        if (fp_out == NULL) { exit(1); }

        int active_particles = 0;
        for(int i=0; i<num_particles; i++) {
            if (!isnan(particle_states[i].pos[0])) active_particles++;
        }
        fwrite(&active_particles, sizeof(int), 1, fp_out);
        for (int i=0; i<num_particles; i++) {
            if (!isnan(particle_states[i].pos[0])) {
                fwrite(&particle_states[i], sizeof(mass_particle_state_t), 1, fp_out);
            }
        }
        fclose(fp_out);
        snapshot_count++;

        if (current_t >= commondata->t_final) break;

        // Evolve all particles for one snapshot interval.
        const double t_next_snapshot = current_t + commondata->snapshot_every_t;
        #pragma omp parallel for
        for (int i = 0; i < num_particles; i++) {
            if (isnan(particle_states[i].pos[0])) continue;

            double y_particle[8];
            y_particle[0] = current_t;
            y_particle[1] = particle_states[i].pos[0];
            y_particle[2] = particle_states[i].pos[1];
            y_particle[3] = particle_states[i].pos[2];
            y_particle[4] = particle_states[i].u[0];
            y_particle[5] = particle_states[i].u[1];
            y_particle[6] = particle_states[i].u[2];
            y_particle[7] = particle_states[i].u[3];

            int status = integrate_single_particle(commondata, params, metric, y_particle[0], t_next_snapshot, y_particle);

            if (status != 0) {
                particle_states[i].pos[0] = NAN; // Mark particle as terminated.
            } else {
                particle_states[i].pos[0] = y_particle[1];
                particle_states[i].pos[1] = y_particle[2];
                particle_states[i].pos[2] = y_particle[3];
                particle_states[i].u[0] = y_particle[4];
                particle_states[i].u[1] = y_particle[5];
                particle_states[i].u[2] = y_particle[6];
                particle_states[i].u[3] = y_particle[7];
            }
        }
    }
    free(particle_states);
    """
    cfc.register_CFunction(
        includes=includes, desc=desc, name=name, params=params, body=body
    )


def main_c_function() -> None:
    """
    Generate the main() C function, the entry point for the executable.

    This function acts as the master orchestrator. It initializes parameters,
    determines the run mode (debug or production), and dispatches to the
    appropriate top-level orchestrator.

    UPDATED: After the production run, this function now optionally calls
    the Python script itself in post-processing mode to generate k-d trees.
    """
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<string.h>",
        "<stdio.h>",
        "<stdlib.h>",  # Needed for system()
    ]
    desc = "Main entry point for the massive particle geodesic integrator."
    cfunc_type = "int"
    name = "main"
    params = "int argc, const char *argv[]"

    # The body for debug mode is long, so we construct it first.
    debug_body = r"""
        particle_initial_state_t initial_state;
        initial_state.id = 0;

        const char *filename = "particle_debug_initial_conditions.txt";
        FILE *fp_in = fopen(filename, "r");

        if (fp_in == NULL) {
            printf("File '%s' not found. Creating it with default values.\\n", filename);
            fp_in = fopen(filename, "w");
            if (fp_in == NULL) { fprintf(stderr, "Error: Could not create '%s'.\\n", filename); return 1; }
            fprintf(fp_in, "# Format: t_initial pos_x pos_y pos_z u_x u_y u_z\\n");
            fprintf(fp_in, "0.0 10.0  0.0  0.0   0.0  0.363232  0.0\\n");
            fclose(fp_in);
            fp_in = fopen(filename, "r");
            if (fp_in == NULL) { fprintf(stderr, "Error: Could not re-open '%s'.\\n", filename); return 1; }
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

        printf("--- Single Particle Debug Run ---\\n");
        printf("  pos = (t=%.4f, x=%.4f, y=%.4f, z=%.4f)\\n", initial_state.pos[0], initial_state.pos[1], initial_state.pos[2], initial_state.pos[3]);

        double y_start[8], y_final[8];
        set_initial_conditions_massive(&initial_state, &commondata, &params, y_start);

        printf("\\nInitial State Vector (y_start):\\n");
        printf("  t=%.2f, x=%.2f, y=%.2f, z=%.2f\\n", y_start[0], y_start[1], y_start[2], y_start[3]);
        printf("  u^t=%.4f, u^x=%.4f, u^y=%.4f, u^z=%.4f\\n\\n", y_start[4], y_start[5], y_start[6], y_start[7]);

        if(commondata.perform_conservation_check) {
            double E_i, Lx_i, Ly_i, Lz_i, Q_i, E_f, Lx_f, Ly_f, Lz_f, Q_f;
            check_conservation_massive(&commondata, &params, &metric, y_start, &E_i, &Lx_i, &Ly_i, &Lz_i, &Q_i);
            integrate_single_particle_DEBUG(&commondata, &params, &metric, y_start, y_final);
            check_conservation_massive(&commondata, &params, &metric, y_final, &E_f, &Lx_f, &Ly_f, &Lz_f, &Q_f);
        } else {
            integrate_single_particle_DEBUG(&commondata, &params, &metric, y_start, y_final);
        }
        printf("\\nDebug run finished. Trajectory saved to 'massive_particle_path.txt'.\\n");
    """

    # This is the main template for the C main() function body.
    body = f"""
    // Step 1: Initialize structs and parameters.
    commondata_struct commondata;
    params_struct params; // This struct is currently unused but required by function signatures.
    metric_params metric;

    commondata_struct_set_to_default(&commondata);
    cmdline_input_and_parfile_parser(&commondata, argc, argv);

    // Step 2: Set Metric Type Enum. Kerr-Schild is used for both.
    metric.type = (commondata.a_spin == 0.0) ? Schwarzschild : Kerr;

    // Step 3: Check the run mode and execute the appropriate logic.
    if (commondata.run_in_debug_mode) {{
        /***********************************/
        /*** SINGLE-PARTICLE DEBUG MODE ***/
        /***********************************/
        {debug_body}
    }} else {{
        /***********************************/
        /*** FULL DISK PRODUCTION MODE ***/
        /***********************************/
        printf("----------------------------------------\\n");
        printf("Massive Particle Integrator\\n");
        printf("----------------------------------------\\n");
        printf("Metric Settings:\\n");
        printf("  Metric Type             = %s (a=%.2f, M=%.2f)\\n", (metric.type == Kerr) ? "Kerr" : "Schwarzschild", commondata.a_spin, commondata.M_scale);
        printf("\\nIntegration Settings:\\n");
        printf("  Max Integration Time    = %.1f M\\n", commondata.t_final);
        printf("  Snapshot Every          = %.1f M\\n", commondata.snapshot_every_t);
        printf("\\nInitial Conditions:\\n");
        printf("  Generator Type          = %s\\n", commondata.initial_conditions_type);
        printf("  Num Particles (r x phi) = %d x %d\\n", commondata.disk_num_r, commondata.disk_num_phi);
        printf("  Disk Radial Min/Max     = %.2f / %.2f M\\n", commondata.disk_r_min, commondata.disk_r_max);
        printf("----------------------------------------\\n\\n");

        run_mass_integrator_production(&commondata, &params, &metric);
        printf("\\nProduction run finished successfully.\\n");

        // *** NEW: POST-PROCESSING STEP ***
        // If enabled, call this script again in post-processing mode.
        if (commondata.generate_kdtree_files) {{
            printf("\\n--- Starting Post-Processing: Generating K-d Tree files ---\\n");
            // Construct the command to call the Python script.
            // We assume the python script is in the parent directory of the executable.
            char command[1024];
            snprintf(command, sizeof(command), "python3 ../mass_geodesic_integrator.py --post-process --outdir ..");

            // Execute the command.
            int status = system(command);
            if (status != 0) {{
                fprintf(stderr, "Error: Post-processing script failed with exit code %d.\\n", status);
            }} else {{
                printf("--- Post-Processing Complete ---\\n");
            }}
        }}
    }}

    return 0;
    """
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )


# ##############################################################################
# PART 4: POST-PROCESSING PYTHON LOGIC (K-D TREE GENERATION)
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

    This function contains the high-performance numerical logic. It takes
    particle positions and recursively partitions them, returning the reordering
    map and the metadata (split axes) for the tree.

    Args:
        particle_pos_array: A NumPy array of shape (N, 3) with particle positions.
        num_particles: The number of particles.

    Returns:
        A tuple (reordering_indices, node_metadata), where reordering_indices
        maps old indices to new, and node_metadata stores the split axis for each node.
    """
    # These will store the final results.
    reordering_indices = np.zeros(num_particles, dtype=np.int32)
    node_metadata = np.full(num_particles, -1, dtype=np.int32)

    # Simulate a queue for the breadth-first build using NumPy arrays.
    q_target_idx = np.zeros(num_particles, dtype=np.int32)
    q_slice_start = np.zeros(num_particles, dtype=np.int32)
    q_slice_end = np.zeros(num_particles, dtype=np.int32)
    q_head = 0
    q_tail = 0

    # The array of original indices that we will sort and partition.
    indices_array = np.arange(num_particles, dtype=np.int32)

    if num_particles > 0:
        # Enqueue the root item: target index 0, and the slice representing all particles.
        q_target_idx[q_tail] = 0
        q_slice_start[q_tail] = 0
        q_slice_end[q_tail] = num_particles
        q_tail += 1

    while q_head < q_tail:
        # Dequeue an item.
        target_idx = q_target_idx[q_head]
        start = q_slice_start[q_head]
        end = q_slice_end[q_head]
        q_head += 1

        # Determine the best axis to split on by finding the widest dimension.
        current_indices = indices_array[start:end]
        current_positions = particle_pos_array[current_indices]

        min_coord = current_positions[0]
        max_coord = current_positions[0]
        for i in range(1, len(current_positions)):
            pos = current_positions[i]
            min_coord = np.minimum(min_coord, pos)
            max_coord = np.maximum(max_coord, pos)

        spread = max_coord - min_coord
        split_axis = np.argmax(spread)

        # Partition the current slice of the index array based on the median.
        sorted_indices_on_axis = current_indices[
            np.argsort(particle_pos_array[current_indices, split_axis])
        ]
        indices_array[start:end] = sorted_indices_on_axis
        median_offset = len(sorted_indices_on_axis) // 2
        median_original_index = sorted_indices_on_axis[median_offset]

        # Place the pivot particle's metadata.
        reordering_indices[target_idx] = median_original_index
        node_metadata[target_idx] = split_axis

        # Enqueue the left and right children for the next level of the tree.
        left_child_idx = 2 * target_idx + 1
        right_child_idx = 2 * target_idx + 2

        if median_offset > 0 and left_child_idx < num_particles:
            q_target_idx[q_tail] = left_child_idx
            q_slice_start[q_tail] = start
            q_slice_end[q_tail] = start + median_offset
            q_tail += 1

        if (
            median_offset + 1 < len(sorted_indices_on_axis)
            and right_child_idx < num_particles
        ):
            q_target_idx[q_tail] = right_child_idx
            q_slice_start[q_tail] = start + median_offset + 1
            q_slice_end[q_tail] = end
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
    # Define the NumPy dtype that matches the C struct 'mass_particle_state_t'.
    raw_dtype = np.dtype(
        [
            ("id", np.int32),
            ("pos", "f8", (3,)),
            ("u", "f8", (4,)),
            ("lambda_rest", "f8"),
            ("j_intrinsic", "f4"),
        ]
    )

    # Step 1: Load the raw snapshot data.
    with open(raw_snapshot_file, "rb") as f:
        num_particles = np.fromfile(f, dtype=np.int32, count=1)[0]
        if num_particles == 0:
            print(f"  Snapshot {os.path.basename(raw_snapshot_file)} is empty, skipping.")
            return
        particle_data = np.fromfile(f, dtype=raw_dtype, count=num_particles)

    # Step 2: Call the fast JIT-compiled kernel with only the position data.
    particle_positions = np.ascontiguousarray(particle_data["pos"])
    reordering_map, node_metadata = build_implicit_kdtree_kernel(
        particle_positions, num_particles
    )

    # Use the returned map to reorder the full, original particle data array.
    reordered_particles = particle_data[reordering_map]

    # Step 3: Save the two parallel arrays to the output .kdtree.bin file.
    output_filename = os.path.join(
        output_dir, os.path.basename(raw_snapshot_file).replace(".bin", ".kdtree.bin")
    )
    with open(output_filename, "wb") as f:
        # Write header: number of particles and dimensions.
        f.write(np.uint64(num_particles))
        f.write(np.uint64(3))
        # Write the node metadata array (split axes).
        f.write(node_metadata.tobytes())
        # Write the reordered particle data payload.
        f.write(reordered_particles.tobytes())


def run_kdtree_preprocessor_parallel(
    snapshot_input_dir: str, snapshot_output_dir: str
) -> None:
    """
    Find all raw .bin snapshot files and process them in parallel.

    This function uses a multiprocessing pool to run the k-d tree build
    kernel on all available snapshot files concurrently.

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

    # Use functools.partial to create a task function with the output_dir argument fixed.
    task_function = partial(
        build_and_save_kdtree_snapshot_fast, output_dir=snapshot_output_dir
    )

    start_time = time.time()
    with Pool(processes=num_processes) as pool:
        # Use tqdm to create a progress bar for the parallel processing.
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
# PART 5: FINAL PROJECT ASSEMBLY
# ##############################################################################

def register_custom_structures_and_params() -> None:
    """
    Generate C code for all custom structs and enums, then register them.

    This function defines the C `typedef`s for all custom data structures
    (structs and enums) used by the project and registers this block of C code
    to be written into the master header file, `BHaH_defines.h`.
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
typedef enum { Schwarzschild, Kerr, Numerical, Schwarzschild_Standard } Metric_t;
typedef struct { Metric_t type; } metric_params;

// Carrier struct to pass simulation data through GSL's generic void* pointer.
typedef struct {
    const commondata_struct *commondata;
    const params_struct *params;
    const metric_params *metric;
} gsl_params;

// Struct for reading initial conditions for a single particle in debug mode.
typedef struct {
    int id;
    double pos[4];
    double u_spatial[3];
} particle_initial_state_t;

// Struct defining the binary format for the production-run snapshot files.
// `__attribute__((packed))` ensures no padding bytes are inserted by the compiler,
// guaranteeing a portable binary file format.
typedef struct {
    int id;
    double pos[3];          // x, y, z position
    double u[4];            // u^t, u^x, u^y, u^z
    double lambda_rest;     // Rest-frame emission wavelength (nm)
    float j_intrinsic;      // Rest-frame intrinsic emissivity (intensity)
} __attribute__((packed)) mass_particle_state_t;
"""
    Bdefines_h.register_BHaH_defines(
        "data_structures",
        f"{metric_struct_str}\n{connections_struct_str}\n{other_structs}",
    )

# ##############################################################################
# PART 6: MAIN EXECUTION BLOCK
# ##############################################################################

if __name__ == "__main__":
    # Step 6.A: Set up command-line argument parser.
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

        # Step 6.B: Set up project directories and core nrpy parameters.
        project_name = "mass_integrator"
        project_dir = os.path.join(args.outdir, project_name)

        # Clean the project directory to ensure a fresh build.
        shutil.rmtree(project_dir, ignore_errors=True)
        os.makedirs(project_dir, exist_ok=True)

        # Set the BHaH infrastructure for nrpy.
        par.set_parval_from_str("Infrastructure", "BHaH")

        # Step 6.C: Register all CodeParameters for the mass integrator.
        # --- Physical Parameters ---
        M_scale = CPs.register_CodeParameter(
            "REAL", __name__, "M_scale", 1.0, commondata=True, add_to_parfile=True
        )
        a_spin = CPs.register_CodeParameter(
            "REAL", __name__, "a_spin", 0.9, commondata=True, add_to_parfile=True
        )
        # --- Post-Processing Control ---
        CPs.register_CodeParameter(
            "bool",
            __name__,
            "generate_kdtree_files",
            True,
            commondata=True,
            add_to_parfile=True,
            description="If true, generate .kdtree.bin files after the simulation.",
        )
        # --- Integration & Termination ---
        CPs.register_CodeParameters(
            "REAL",
            __name__,
            ["t_max_integration", "r_escape", "ut_max"],
            [2000.0, 1500.0, 1e3],
            commondata=True,
            add_to_parfile=True,
        )
        # --- Debugging & Validation ---
        CPs.register_CodeParameters(
            "bool",
            __name__,
            ["perform_conservation_check", "run_in_debug_mode"],
            True,
            commondata=True,
            add_to_parfile=True,
        )
        # --- Disk & Output Parameters ---
        CPs.register_CodeParameters(
            "REAL",
            __name__,
            [
                "disk_lambda_rest_at_r_min",
                "disk_r_min",
                "disk_r_max",
                "snapshot_every_t",
                "t_final",
            ],
            [656.3, 6.0, 25.0, 10.0, 2000.0],
            commondata=True,
            add_to_parfile=True,
        )
        CPs.register_CodeParameters(
            "int",
            __name__,
            ["disk_num_r", "disk_num_phi"],
            [100, 200],
            commondata=True,
            add_to_parfile=True,
        )
        CPs.register_CodeParameter(
            "char[100]",
            __name__,
            "output_folder",
            "output",
            commondata=True,
            add_to_parfile=True,
        )
        # --- Galaxy Shape Parameters ---
        CPs.register_CodeParameters(
            "REAL",
            __name__,
            [
                "bar_length",
                "bar_aspect_ratio",
                "bulge_radius",
                "arm_particle_density",
                "arm_clumpiness_factor",
                "arm_clump_size",
                "bar_density_factor",
                "bulge_density_factor",
                "spiral_galaxy_arm_tightness",
            ],
            [5.0, 0.25, 1.5, 0.3, 8.0, 0.5, 2.0, 3.0, 0.2],
            commondata=True,
            add_to_parfile=True,
        )
        CPs.register_CodeParameter(
            "int",
            __name__,
            "spiral_galaxy_num_arms",
            2,
            commondata=True,
            add_to_parfile=True,
        )
        # --- Initial Conditions Type Selector ---
        CPs.register_CodeParameter(
            "char[100]",
            __name__,
            "initial_conditions_type",
            "KeplerianDisk",
            commondata=True,
            add_to_parfile=True,
        )

        # Step 6.D: Execute the symbolic workflow.
        g4DD_kerr, xx_kerr = georec.define_kerr_metric_Cartesian_Kerr_Schild(
            M_scale, a_spin
        )
        g4DD_dD_kerr = georec.derivative_g4DD(g4DD_kerr, xx_kerr)
        Gamma4UDD_kerr = georec.four_connections(g4DD_kerr, g4DD_dD_kerr)

        ut_expr_stable, uphi_expr_stable = symbolic_ut_uphi_from_r_stable(
            M_scale, a_spin
        )
        E_expr, L_exprs, Q_expr = symbolic_conserved_quantities()
        pos_rhs, vel_rhs = geodesic_eom_rhs()

        # Step 6.E: Call all C-generating Python functions to register C code.
        g4dd_kerr_schild(g4DD_kerr)
        con_kerr_schild(Gamma4UDD_kerr)
        g4dd_metric()
        connections()
        calculate_ode_rhs_massive(pos_rhs, vel_rhs)
        check_conservation_massive(E_expr, L_exprs, Q_expr)
        ode_gsl_wrapper_massive()
        calculate_ut_uphi_from_r(ut_expr_stable, uphi_expr_stable)
        set_initial_conditions_massive()
        generate_disk_initial_conditions()
        generate_spiral_galaxy_initial_conditions()
        generate_barred_flocculent_spiral_ic()
        integrate_single_particle()
        integrate_single_particle_DEBUG()
        run_mass_integrator_production()
        main_c_function()

        # Step 6.F: Assemble the final C project on disk.
        register_custom_structures_and_params()
        CPs.write_CodeParameters_h_files(project_dir=project_dir)
        CPs.register_CFunctions_params_commondata_struct_set_to_default()
        cmdline_input_and_parfiles.generate_default_parfile(
            project_dir=project_dir, project_name=project_name
        )
        cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
            project_name=project_name,
            cmdline_inputs=[
                "M_scale",
                "a_spin",
                "generate_kdtree_files",
                "t_max_integration",
                "r_escape",
                "ut_max",
                "perform_conservation_check",
                "run_in_debug_mode",
                "initial_conditions_type",
                "spiral_galaxy_num_arms",
                "spiral_galaxy_arm_tightness",
                "bar_length",
                "bar_aspect_ratio",
                "bulge_radius",
                "arm_particle_density",
                "arm_clumpiness_factor",
                "arm_clump_size",
                "bar_density_factor",
                "bulge_density_factor",
            ],
        )
        Bdefines_h.output_BHaH_defines_h(project_dir=project_dir)

        addl_CFLAGS = ["-Wall -Wextra -g $(shell gsl-config --cflags) -fopenmp"]
        addl_libraries = ["$(shell gsl-config --libs) -fopenmp"]

        Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
            project_dir=project_dir,
            project_name=project_name,
            exec_or_library_name="mass_integrator",
            addl_CFLAGS=addl_CFLAGS,
            addl_libraries=addl_libraries,
        )

        # Step 6.G: Print a final confirmation message.
        print(f"\nFinished! A C project has been generated in '{project_dir}/'")
        print(f"To build, navigate to this directory in your terminal and type 'make'.")
        print(f"To run, type './mass_integrator'.")
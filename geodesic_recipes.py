# nrpy/helpers/geodesic_recipes.py

"""
Provides a library of shared, pure-math symbolic recipes for geodesic motion.

This module contains Python functions that take SymPy expressions as input
(e.g., a metric tensor) and return new SymPy expressions (e.g., Christoffel symbols).
It has no C-code generation logic and is designed to be imported by dedicated
project generator scripts for both massive and massless particles.

Author: Dalton J. Moone
"""

from typing import List, Tuple

import sympy as sp
import nrpy.indexedexp as ixp


def define_kerr_metric_Cartesian_Kerr_Schild(
    M_scale: sp.Symbol, a_spin: sp.Symbol
) -> Tuple[List[List[sp.Expr]], List[sp.Symbol]]:
    """
    Define the Kerr metric tensor in Cartesian Kerr-Schild coordinates.

    This coordinate system is regular everywhere, including at the event horizon,
    making it highly advantageous for numerical work. The metric is constructed
    using the formula:
        g_{μν} = η_{μν} + 2H l_μ l_ν
    where η_μν is the Minkowski metric, l_μ is a null vector, and H is a scalar
    function. A key feature is that if the spin parameter 'a' is set to zero,
    this metric automatically and exactly reduces to the Schwarzschild metric.

    Args:
        M_scale: The symbolic mass of the black hole.
        a_spin: The symbolic dimensionless spin parameter of the black hole.

    Returns:
        A tuple (g4DD, xx), where g4DD is the symbolic 4x4 metric tensor
        and xx is the list of symbolic coordinate variables [t, x, y, z].
    """
    # Define the symbolic coordinates using the 'y[i]' convention for the integrator.
    # This ensures the generated C code will correctly map to the state vector array.
    t, x, y, z = sp.symbols("y[0] y[1] y[2] y[3]", real=True)
    xx: List[sp.Symbol] = [t, x, y, z]

    # Define intermediate quantities based on Cartesian coordinates.
    r2 = x**2 + y**2 + z**2
    r = sp.sqrt(r2)

    # Define the Kerr-Schild null vector l_μ.
    l_down = ixp.zerorank1(dimension=4)
    l_down[0] = sp.sympify(1)  # Explicitly convert to a SymPy object for robustness.
    l_down[1] = (r * x + a_spin * y) / (r2 + a_spin**2)
    l_down[2] = (r * y - a_spin * x) / (r2 + a_spin**2)
    l_down[3] = z / r

    # Define the scalar function H.
    H = (M_scale * r**3) / (r**4 + a_spin**2 * z**2)

    # Construct the Kerr-Schild metric g_μν = η_μν + 2H * l_μ * l_ν,
    # where η_μν is the Minkowski metric diag(-1, 1, 1, 1).
    g4DD: List[List[sp.Expr]] = ixp.zerorank2(dimension=4)
    for mu in range(4):
        for nu in range(4):
            # Start with the flat Minkowski metric component.
            eta_mu_nu = sp.sympify(0)
            if mu == nu:
                eta_mu_nu = sp.sympify(1)
            if mu == 0 and nu == 0:
                eta_mu_nu = sp.sympify(-1)

            # Add the Kerr-Schild correction term.
            g4DD[mu][nu] = eta_mu_nu + 2 * H * l_down[mu] * l_down[nu]

    return g4DD, xx


def define_schwarzschild_metric_cartesian(
    M_scale: sp.Symbol,
) -> Tuple[List[List[sp.Expr]], List[sp.Symbol]]:
    """
    Define the Schwarzschild metric tensor directly in Cartesian coordinates.

    This function uses the standard textbook formula and serves as an important
    cross-validation case for the more general Kerr-Schild implementation.
    The metric components are:
    - g_tt = -(1 - 2M/r)
    - g_ij = δ_ij + (2M/r) * (x_i * x_j / r^2)

    Args:
        M_scale: The symbolic mass of the black hole.

    Returns:
        A tuple (g4DD, xx), where g4DD is the symbolic 4x4 metric tensor
        and xx is the list of symbolic coordinate variables [t, x, y, z].
    """
    # Define Cartesian coordinates.
    t, x, y, z = sp.symbols("y[0] y[1] y[2] y[3]", real=True)
    xx: List[sp.Symbol] = [t, x, y, z]

    # Define r in terms of Cartesian coordinates.
    r = sp.sqrt(x**2 + y**2 + z**2)

    # Initialize the metric tensor with zeros.
    g4DD: List[List[sp.Expr]] = ixp.zerorank2(dimension=4)

    # Set the time-time component, g_tt.
    g4DD[0][0] = -(1 - 2 * M_scale / r)

    # Set the spatial components g_ij = δ_ij + (2M/r) * (x_i * x_j / r^2).
    x_i = [x, y, z]
    for i in range(3):
        for j in range(3):
            delta_ij = sp.sympify(1) if i == j else sp.sympify(0)
            # The indices for the 4D metric g4DD are offset by 1 from the 3D spatial indices.
            g4DD[i + 1][j + 1] = delta_ij + (2 * M_scale / r) * (
                x_i[i] * x_i[j] / (r**2)
            )

    # Explicitly set the time-space components to zero as SymPy objects.
    for i in range(1, 4):
        g4DD[0][i] = g4DD[i][0] = sp.sympify(0)

    return g4DD, xx


def derivative_g4DD(
    g4DD: List[List[sp.Expr]], xx: List[sp.Symbol]
) -> List[List[List[sp.Expr]]]:
    """
    Compute the symbolic first derivatives of the metric tensor, g_{μν,α}.

    This function iterates through all components of the metric tensor and
    symbolically calculates the partial derivative of each component with
    respect to each coordinate using SymPy's sp.diff(). The resulting
    quantity is defined as:
        g_{μν,α} ≡ ∂g_{μν} / ∂x^α

    Args:
        g4DD: The symbolic 4x4 metric tensor.
        xx: The list of symbolic coordinate variables [t, x, y, z].

    Returns:
        A rank-3 symbolic tensor containing the derivatives g_{μν,α}.
    """
    g4DD_dD: List[List[List[sp.Expr]]] = ixp.zerorank3(dimension=4)
    for nu in range(4):
        for mu in range(4):
            for alpha in range(4):
                g4DD_dD[nu][mu][alpha] = sp.diff(g4DD[nu][mu], xx[alpha])
    return g4DD_dD


def four_connections(
    g4DD: List[List[sp.Expr]], g4DD_dD: List[List[List[sp.Expr]]]
) -> List[List[List[sp.Expr]]]:
    """
    Compute and simplify Christoffel symbols from the metric and its derivatives.

    This function implements the standard formula for the Christoffel symbols
    of the second kind, Γ^{δ}_{μν}:
        Γ^{δ}_{μν} = (1/2) g^{δα} (g_{να,μ} + g_{μα,ν} - g_{μν,α})

    It first computes the inverse metric g^{μν} and then applies the formula.
    The result is simplified using sp.trigsimp(), which is highly effective
    for the Kerr-Schild metric.

    Args:
        g4DD: The symbolic 4x4 metric tensor.
        g4DD_dD: The symbolic rank-3 tensor of metric derivatives.

    Returns:
        A rank-3 symbolic tensor containing the Christoffel symbols Γ^{δ}_{μν}.
    """
    # Initialize the output tensor for the Christoffel symbols.
    Gamma4UDD: List[List[List[sp.Expr]]] = ixp.zerorank3(dimension=4)

    # Compute the inverse metric tensor g^{μν}. The determinant is not needed here.
    g4UU, _ = ixp.symm_matrix_inverter4x4(g4DD)

    # Iterate over all indices to compute each component of the Christoffel symbol.
    for mu in range(4):
        for nu in range(4):
            for delta in range(4):
                # Initialize a sum for the Einstein summation over the dummy index alpha.
                sum_term = sp.sympify(0)
                for alpha in range(4):
                    sum_term += (
                        g4DD_dD[nu][alpha][mu]
                        + g4DD_dD[mu][alpha][nu]
                        - g4DD_dD[mu][nu][alpha]
                    )
                
                # Multiply by (1/2)g^{δα} and assign to the final tensor.
                Gamma4UDD[delta][mu][nu] = sp.Rational(1, 2) * g4UU[delta][alpha] * sum_term
                
                # Use sp.trigsimp() to simplify the resulting expression. This is a key
                # optimization for metrics with trigonometric functions, like Kerr-Schild.
                Gamma4UDD[delta][mu][nu] = sp.trigsimp(Gamma4UDD[delta][mu][nu])

    return Gamma4UDD
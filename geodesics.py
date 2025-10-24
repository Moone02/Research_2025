"""
Construct symbolic expressions for geodesic equations in various spacetime backgrounds.

This module provides a class-based structure for generating the symbolic
expressions needed for integrating massive particle geodesics. It adheres to the
nrpy.equations "gold standard" by encapsulating equations within a class,
using a caching dictionary to prevent redundant symbolic computations, and
including a comprehensive validation suite.

The commenting style is designed to be pedagogical, explaining each mathematical
step for an audience familiar with advanced undergraduate physics, mirroring the
standard set by modules like BSSN_RHSs.py.

Author: Dalton J. Moone
"""
# Step 0.a: Import standard Python modules
import doctest
import os
import sys
from typing import Dict, List, Tuple

# Step 0.b: Import third-party modules
import sympy as sp

# Step 0.c: Import NRPy+ core modules
import nrpy.indexedexp as ixp
import nrpy.validate_expressions.validate_expressions as ve

# Step 0.d: Define global symbolic parameters for mass and spin.
# These are used by all functions within the class.
M_scale = sp.Symbol("M_scale", real=True)
a_spin = sp.Symbol("a_spin", real=True)


class GeodesicEquations:
    """
    Generate and store symbolic expressions for geodesic motion.

    This class is instantiated with a specific spacetime and particle type.
    It then drives a symbolic workflow to compute all necessary quantities,
    such as the metric, Christoffel symbols, equations of motion, and
    conserved quantities, storing them as instance attributes.
    """

    def __init__(
        self,
        spacetime: str,
        particle_type: str = "massive",
    ) -> None:
        """
        Initialize and compute all symbolic expressions for geodesic equations.

        Args:
            spacetime: The spacetime to use (e.g., "KerrSchild", "Schwarzschild_Cartesian").
            particle_type: The type of particle, currently supports "massive".
        """
        self.spacetime = spacetime
        self.particle_type = particle_type

        # Step 1: Define the metric based on the chosen spacetime.
        if self.spacetime == "KerrSchild":
            self.g4DD, self.xx = self._define_kerr_metric_Cartesian_Kerr_Schild()
        elif self.spacetime == "Schwarzschild_Cartesian":
            self.g4DD, self.xx = self._define_schwarzschild_metric_cartesian()
        else:
            raise ValueError(f"Spacetime '{self.spacetime}' is not supported.")

        # Step 2: Compute metric derivatives and Christoffel symbols.
        self.g4DD_dD = self._derivative_g4DD(self.g4DD, self.xx)
        self.Gamma4UDD = self._four_connections(self.g4DD, self.g4DD_dD)

        # Step 3: Generate particle-specific equations.
        if self.particle_type == "massive":
            self.geodesic_rhs = self._geodesic_eom_rhs()
            self.ut_from_vel_expr = self._ut_massive()
            (
                self.ut_stable_expr,
                self.uphi_stable_expr,
            ) = self._symbolic_ut_uphi_from_r_stable()
            (
                self.E_expr,
                self.L_exprs,
                self.Q_expr,
            ) = self._symbolic_conserved_quantities()
        else:
            raise ValueError(f"Particle type '{self.particle_type}' is not supported.")

    @staticmethod
    def _define_kerr_metric_Cartesian_Kerr_Schild() -> (
        Tuple[List[List[sp.Expr]], List[sp.Symbol]]
    ):
        """
        Define the Kerr metric in Cartesian Kerr-Schild coordinates.

        This coordinate system is regular everywhere, including at the event horizon,
        making it ideal for numerical simulations. The metric is constructed as a
        perturbation of the flat Minkowski metric: g_munu = eta_munu + 2H * l_mu * l_nu.
        If the spin parameter 'a' is set to zero, this metric exactly reduces to
        the Schwarzschild metric in the same coordinate system.

        Returns:
            A tuple (g4DD, xx), where g4DD is the symbolic 4x4 metric tensor
            and xx is the list of symbolic coordinate variables (t, x, y, z).
        """
        t, x, y, z = sp.symbols("y[0] y[1] y[2] y[3]", real=True)
        xx = [t, x, y, z]
        r2 = x**2 + y**2 + z**2
        r = sp.sqrt(r2)
        l_down = ixp.zerorank1(dimension=4)
        l_down[0] = sp.sympify(1)
        l_down[1] = (r * x + a_spin * y) / (r2 + a_spin**2)
        l_down[2] = (r * y - a_spin * x) / (r2 + a_spin**2)
        l_down[3] = z / r
        H = (M_scale * r**3) / (r**4 + a_spin**2 * z**2)
        eta4DD = ixp.zerorank2(dimension=4)
        eta4DD[0][0] = -1
        eta4DD[1][1] = eta4DD[2][2] = eta4DD[3][3] = 1
        g4DD = ixp.zerorank2(dimension=4)
        for mu in range(4):
            for nu in range(4):
                g4DD[mu][nu] = eta4DD[mu][nu] + 2 * H * l_down[mu] * l_down[nu]
        return g4DD, xx

    @staticmethod
    def _define_schwarzschild_metric_cartesian() -> (
        Tuple[List[List[sp.Expr]], List[sp.Symbol]]
    ):
        """
        Define the Schwarzschild metric in standard Cartesian coordinates.

        This serves as an independent implementation for validation purposes.
        The metric components are:
        g_tt = -(1 - 2M/r)
        g_ij = delta_ij + (2M/r) * (x_i * x_j / r^2)

        Returns:
            A tuple (g4DD, xx), where g4DD is the symbolic 4x4 metric tensor
            and xx is the list of symbolic coordinate variables.
        """
        t, x, y, z = sp.symbols("y[0] y[1] y[2] y[3]", real=True)
        xx = [t, x, y, z]
        r = sp.sqrt(x**2 + y**2 + z**2)
        x_i = [x, y, z]
        g4DD = ixp.zerorank2(dimension=4)
        g4DD[0][0] = -(1 - 2 * M_scale / r)
        for i in range(3):
            for j in range(3):
                delta_ij = sp.sympify(1) if i == j else sp.sympify(0)
                g4DD[i + 1][j + 1] = delta_ij + (2 * M_scale / r) * (
                    x_i[i] * x_i[j] / (r**2)
                )
        return g4DD, xx

    @staticmethod
    def _derivative_g4DD(
        g4DD: List[List[sp.Expr]], xx: List[sp.Symbol]
    ) -> List[List[List[sp.Expr]]]:
        """
        Compute symbolic first derivatives of the metric tensor, g_{munu,alpha}.

        Args:
            g4DD: The symbolic 4x4 metric tensor.
            xx: The list of symbolic coordinate variables.

        Returns:
            A rank-3 tensor g_{munu,alpha}.
        """
        g4DD_dD = ixp.zerorank3(dimension=4)
        for mu in range(4):
            for nu in range(4):
                for alpha in range(4):
                    g4DD_dD[mu][nu][alpha] = sp.diff(g4DD[mu][nu], xx[alpha])
        return g4DD_dD

    @staticmethod
    def _four_connections(
        g4DD: List[List[sp.Expr]], g4DD_dD: List[List[List[sp.Expr]]]
    ) -> List[List[List[sp.Expr]]]:
        """
        Compute and simplify Christoffel symbols of the second kind, Gamma^delta_{mu,nu}.

        The formula is:
        Gamma^delta_{mu,nu} = (1/2) * g^{delta,alpha} * (g_{nu,alpha,mu} + g_{mu,alpha,nu} - g_{mu,nu,alpha})

        Args:
            g4DD: The symbolic 4x4 metric tensor.
            g4DD_dD: The symbolic derivatives of the metric tensor.

        Returns:
            A rank-3 tensor Gamma^delta_{mu,nu}.
        """
        g4UU, _ = ixp.symm_matrix_inverter4x4(g4DD)
        Gamma4UDD = ixp.zerorank3(dimension=4)
        for delta in range(4):
            for mu in range(4):
                for nu in range(4):
                    term = sp.sympify(0)
                    for alpha in range(4):
                        term += sp.Rational(1, 2) * g4UU[delta][alpha] * (
                            g4DD_dD[nu][alpha][mu]
                            + g4DD_dD[mu][alpha][nu]
                            - g4DD_dD[mu][nu][alpha]
                        )
                    Gamma4UDD[delta][mu][nu] = sp.trigsimp(term)
        return Gamma4UDD

    @staticmethod
    def _geodesic_eom_rhs() -> List[sp.Expr]:
        """
        Generate symbolic RHS for the 8 geodesic ODEs.

        The system is decomposed into two first-order systems:
        1. dx^alpha/dtau = u^alpha
        2. du^alpha/dtau = -Gamma^alpha_{mu,nu} u^mu u^nu

        Returns:
            A list of 8 symbolic expressions for the RHS of the ODE system.
        """
        ut, ux, uy, uz = sp.symbols("y[4] y[5] y[6] y[7]", real=True)
        pos_rhs = [ut, ux, uy, uz]
        Gamma4UDD = ixp.declarerank3("conn->Gamma4UDD", dimension=4, sym="sym12")
        uU = [ut, ux, uy, uz]
        vel_rhs = ixp.zerorank1(dimension=4)
        for alpha in range(4):
            sum_term = sp.sympify(0)
            for mu in range(4):
                for nu in range(mu, 4):
                    term = Gamma4UDD[alpha][mu][nu] * uU[mu] * uU[nu]
                    if mu != nu:
                        term *= 2
                    sum_term += term
            vel_rhs[alpha] = -sum_term
        return pos_rhs + vel_rhs

    @staticmethod
    def _ut_massive() -> sp.Expr:
        """
        Symbolically derive u^t for a massive particle.

        This solves the timelike normalization condition g_munu u^mu u^nu = -1
        for u^t, which is a quadratic equation: A(u^t)^2 + B(u^t) + C = 0.

        Returns:
            A symbolic expression for u^t.
        """
        u0, u1, u2, u3 = sp.symbols("u0 u1 u2 u3", real=True)
        uU = [u0, u1, u2, u3]
        g4DD = ixp.declarerank2("metric->g", sym="sym01", dimension=4)
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
        discriminant = (2 * sum_g0i_ui) ** 2 - 4 * g4DD[0][0] * (sum_gij_ui_uj + 1)
        answer = (-2 * sum_g0i_ui - sp.sqrt(discriminant)) / (2 * g4DD[0][0])
        return answer

    @staticmethod
    def _symbolic_ut_uphi_from_r_stable() -> Tuple[sp.Expr, sp.Expr]:
        """
        Derive u^t and u^phi for a circular, equatorial orbit using a stable method.

        Returns:
            A tuple containing symbolic expressions for (u^t, u^phi).
        """
        r = sp.symbols("r_initial", real=True)
        Omega = sp.sqrt(M_scale) / (r ** sp.Rational(3, 2) + a_spin * sp.sqrt(M_scale))
        g_tt = -(1 - 2 * M_scale / r)
        g_tphi = -2 * a_spin * M_scale / r
        g_phiphi = r**2 + a_spin**2 + (2 * M_scale * a_spin**2) / r
        ut_squared_inv_denom = g_tt + 2 * g_tphi * Omega + g_phiphi * Omega**2
        ut_squared = -1 / ut_squared_inv_denom
        ut = sp.sqrt(ut_squared)
        uphi = Omega * ut
        return ut, uphi

    @staticmethod
    def _symbolic_conserved_quantities() -> Tuple[sp.Expr, List[sp.Expr], sp.Expr]:
        """
        Generate symbolic recipes for conserved quantities E, L_i, and Q.

        Returns:
            A tuple (E_expr, L_exprs, Q_expr).
        """
        t, x, y, z = sp.symbols("y[0] y[1] y[2] y[3]", real=True)
        pt, px, py, pz = sp.symbols("y[4] y[5] y[6] y[7]", real=True)
        pU = [pt, px, py, pz]
        g4DD = ixp.declarerank2("metric->g", sym="sym01", dimension=4)
        p_t = sp.sympify(0)
        for mu in range(4):
            p_t += g4DD[0][mu] * pU[mu]
        E_expr = -p_t
        p_down = ixp.zerorank1(dimension=4)
        for k in range(1, 4):
            for mu in range(4):
                p_down[k] += g4DD[k][mu] * pU[mu] if k <= mu else g4DD[mu][k] * pU[mu]
        p_x, p_y, p_z = p_down[1], p_down[2], p_down[3]
        L_x, L_y, L_z = y * p_z - z * p_y, z * p_x - x * p_z, x * p_y - y * p_x
        L_exprs = [L_x, L_y, L_z]
        if a_spin == 0:
            Q_expr = L_x**2 + L_y**2 + L_z**2
        else:
            r_sq, rho_sq = x**2 + y**2 + z**2, x**2 + y**2
            xpx_plus_ypy = x * p_x + y * p_y
            p_theta_sq = (
                (z**2 * xpx_plus_ypy**2 / rho_sq)
                - (2 * z * p_z * xpx_plus_ypy)
                + (rho_sq * p_z**2)
            )
            second_term = (z**2 / r_sq) * (
                a_spin**2 * (1 - E_expr**2) + L_z**2 * (r_sq / rho_sq)
            )
            Q_formula = p_theta_sq + second_term
            Q_expr = sp.Piecewise((0, rho_sq < 1e-12), (Q_formula, True))
        return E_expr, L_exprs, Q_expr


class GeodesicEquations_dict(Dict[str, GeodesicEquations]):
    """A caching dictionary for GeodesicEquations instances."""

    def __getitem__(self, key: str) -> GeodesicEquations:
        """
        Get or create a GeodesicEquations instance for a given configuration.

        Args:
            key: A string key in the format "Spacetime_ParticleType".

        Returns:
            A GeodesicEquations instance for the specified configuration.
        """
        if key not in self:
            # FIX: Robustly parse the key. The last part is the particle type,
            #      and everything before it is the spacetime.
            parts = key.split("_")
            if len(parts) < 2:
                raise ValueError(
                    "Invalid key format. Expected 'Spacetime_ParticleType'."
                )
            particle_type = parts[-1]
            spacetime = "_".join(parts[:-1])

            print(
                f"Setting up GeodesicEquations for spacetime='{spacetime}', particle='{particle_type}'..."
            )
            self[key] = GeodesicEquations(
                spacetime=spacetime,
                particle_type=particle_type,
            )
        return super().__getitem__(key)


Geodesic_Equations = GeodesicEquations_dict()


if __name__ == "__main__":
    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    for config in ["KerrSchild_massive", "Schwarzschild_Cartesian_massive"]:
        geo_eqs = Geodesic_Equations[config]
        results_dict = ve.process_dictionary_of_expressions(
            geo_eqs.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{config}",
            results_dict,
        )

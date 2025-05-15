"""Simulate photon trajectories in Schwarzschild spacetime and map images.

This set of functions provides tools for relativistic ray tracing in a 
Schwarzschild spacetime. It includes functions to:
1. Compute photon trajectories by numerically integrating the geodesic equations,
   handling turning points and various termination conditions.
2. Calculate only the crossing points of trajectories with specified target planes,
   optimized for speed when the full trajectory is not needed.
3. Map the origin of photons on a background source plane to an observer's window plane,
   considering both Cartesian and radial parameterizations of the window.
4. Generate datasets of these mappings for various configurations, with support for
   chunking results to manage memory for large computations.
5. Process these datasets along with a source image to render a gravitationally
   lensed view, as seen by an observer near the Schwarzschild black hole.

The core physics is based on the Schwarzschild metric, and the numerical integration
utilizes an adaptive step-size ODE solver from SciPy. Event detection is used for
accurately finding intersections with target planes and handling turning points.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot* com
"""
import numpy as np
from scipy.integrate import solve_ivp
import time # Optional: for timing checks if needed
import sys # For doctest exit status
from typing import Tuple, List, Optional, Union, Any, Callable

# compute_trajectory: Calculates the trajectory based on initial polar coordinates and psi angle.
# t_end is the max K (affine parameter) value to integrate to, r_max is the max radius before integration stops.
# x_stop: If provided (number), integration stops if x (Cartesian coordinate) decreases below this value.
# x_targets: A list/array of x-values. Record y-values when the trajectory crosses each x-value moving from x > x_target to x < x_target.
#            Uses interpolation for higher accuracy at the crossing point.
# num_interp_points: The number of points to use for the final interpolated output trajectory.
# Returns:
#   x, y, r, phi, K (trajectory arrays, interpolated to num_interp_points)
#   crossings_y_at_x_targets (list of lists for y-values at x_targets, calculated accurately at event times)
def compute_trajectory(r_0: float, phi_0: float, M: float, psi: float, r_max: float, x_stop: Optional[float] = None, x_targets: Optional[List[float]] = None, t_end: float = 10**6, num_interp_points: int = 1000) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, List[List[float]]]:
    """Compute a photon trajectory in Schwarzschild spacetime.

    Integrates the geodesic equations for a photon and returns its path in both
    Cartesian (x, y) and polar (r, phi) coordinates, along with the affine
    parameter K. It can also record y-coordinates at specified x-target crossings.

    DocTests:
    >>> # Test a very short, immediate stop case (r_0 already at r_max)
    >>> x, y, r, p, K, cross = compute_trajectory(r_0=10.0, phi_0=0.0, M=1.0, psi=np.pi/2, r_max=10.0, t_end=1, num_interp_points=2)
    >>> len(x)
    1
    >>> np.isclose(x[0], 10.0)
    True
    >>> len(cross)
    0
    >>> # Test simple radial trajectory (b=0) that should hit event horizon quickly
    >>> # M=1, so event horizon at r=2. Start at r_0=3, psi=np.pi (inward radial).
    >>> x_out, y_out, r_out, phi_out, K_out, _ = compute_trajectory(r_0=3.0, phi_0=0.0, M=1.0, psi=np.pi, r_max=100.0, x_stop=0.0, num_interp_points=10)
    >>> r_out[-1] <= (2*1.0 + 1e-6) # Should end at or very near event horizon
    True
    >>> np.allclose(phi_out, 0.0) # Radial trajectory, phi should not change from phi_0
    True
    """
    # --- Basic Setup ---\n
    # Calculate impact parameter b
    # Ensure r_0 is outside the event horizon for a physical starting point for b calculation
    if r_0 <= 2 * M:
        # Handle cases where starting inside or at event horizon - trajectory is ill-defined or trivial
        # For simplicity, return a single point representing the initial state.
        K_out = np.array([0.0])
        r_out = np.array([r_0])
        phi_out = np.array([phi_0 % (2 * np.pi)])
        x_out = r_out * np.cos(phi_out)
        y_out = r_out * np.sin(phi_out)
        num_x_targets_init = len(x_targets) if x_targets is not None else 0
        crossings_y_at_x_targets_out = [[] for _ in range(num_x_targets_init)]
        return x_out, y_out, r_out, phi_out, K_out, crossings_y_at_x_targets_out

    metric_r0_term = (1 - 2 * M / r_0)
    if metric_r0_term <= 0: # Should be caught by r_0 <= 2M, but as a safeguard for sqrt
        K_out = np.array([0.0])
        r_out = np.array([r_0])
        phi_out = np.array([phi_0 % (2 * np.pi)])
        x_out = r_out * np.cos(phi_out)
        y_out = r_out * np.sin(phi_out)
        num_x_targets_init = len(x_targets) if x_targets is not None else 0
        crossings_y_at_x_targets_out = [[] for _ in range(num_x_targets_init)]
        return x_out, y_out, r_out, phi_out, K_out, crossings_y_at_x_targets_out

    b = np.cos(psi) * r_0 / np.sqrt(metric_r0_term)
    sign = [1]  # sign[0] determines dr/dK direction, mutable list for pass-by-reference behavior
    if psi % (2 * np.pi) > np.pi: # Corresponds to inward radial component if psi is angle with local r-direction
        sign[0] *= -1

    # --- ODE Definition ---\n
    def ODE_Solver(K: float, u: np.ndarray) -> List[float]:
        r, phi = u
        if r <= 2 * M: # Photon captured by black hole
            return [0.0, 0.0]
        metric_term = (1 - 2 * M / r)
        term_b_r_sq = (b / r) ** 2
        # Effective radial potential related term for dr/dK
        fr = 1 - metric_term * term_b_r_sq 
        # Use abs() to prevent sqrt domain error just before event triggers for fr=0
        # Event handling (event_fr_zero) should manage the sign change correctly at turning points.
        drdK = sign[0] * np.sqrt(abs(fr))
        dphidK = b / r**2
        return [drdK, dphidK]

    # --- Standard Terminal Event Definitions ---\n
    def event_fr_zero(K: float, u: np.ndarray) -> float: # Index 0: Radial turning point
        r, phi = u
        if r <= (2 * M + 1e-12): return 0.0 # Avoid issues very close to or inside horizon
        metric_term = (1 - 2 * M / r)
        term_b_r_sq = (b / r) ** 2
        return 1 - metric_term * term_b_r_sq # R(r) = 0 condition
    event_fr_zero.terminal = True
    event_fr_zero.direction = 0 # Trigger on vanishing from positive or negative

    def event_r_leq_2M(K: float, u: np.ndarray) -> float: # Index 1: Reaches event horizon (capture)
        r, phi = u
        return r - (2 * M + 1e-12) # Add small epsilon to avoid numerical issues at boundary
    event_r_leq_2M.terminal = True
    event_r_leq_2M.direction = -1 # Trigger when r is decreasing towards 2M

    def event_r_max(K: float, u: np.ndarray) -> float: # Index 2: Reaches maximum radial boundary (escape)
        r, phi = u
        return r - r_max
    event_r_max.terminal = True
    event_r_max.direction = 1 # Trigger when r is increasing towards r_max

    # --- Initialize Event List & Indices ---\n
    active_events: List[Callable] = [event_fr_zero, event_r_leq_2M, event_r_max]
    event_x_stop_index = -1
    x_target_event_indices: List[int] = []

    # --- Optional Terminal Event: x_stop ---\n
    event_x_stop_active = isinstance(x_stop, (int, float))
    if event_x_stop_active:
        def event_x_stop_func(K: float, u: np.ndarray) -> float:
            r, phi = u
            if r < 1e-10: return 1.0 # Avoid issues if r is very small
            x = r * np.cos(phi)
            return x - x_stop # type: ignore[operator]
        event_x_stop_func.terminal = True
        event_x_stop_func.direction = -1 # Trigger when x decreases past x_stop
        event_x_stop_index = len(active_events)
        active_events.append(event_x_stop_func)

    # --- Optional Non-Terminal Events: x_targets ---\n
    num_x_targets = 0
    event_x_targets_active = False
    x_targets_list: List[float] = []
    if x_targets is not None:
        try:
            _ = iter(x_targets) # Check if iterable
            x_targets_list = list(x_targets)
            num_x_targets = len(x_targets_list)
            if num_x_targets > 0: event_x_targets_active = True
        except TypeError:
            print("Warning: x_targets provided but not iterable. Ignoring.")
            num_x_targets = 0

    if event_x_targets_active:
        base_idx = len(active_events)
        for i, xt_val in enumerate(x_targets_list):
            # Use a factory function to correctly capture xt_val in the closure
            def make_event_func(target_x: float) -> Callable:
                def event_x_target_dynamic(K: float, u: np.ndarray) -> float:
                    r, phi = u
                    if r < 1e-10: return 1.0 # Avoid issues if r is very small
                    x = r * np.cos(phi)
                    return x - target_x
                event_x_target_dynamic.terminal = False # Non-terminal event
                event_x_target_dynamic.direction = -1 # Trigger when x decreases past target_x
                return event_x_target_dynamic
            event_func = make_event_func(xt_val)
            active_events.append(event_func)
            x_target_event_indices.append(base_idx + i)

    # --- Initialization for Integration Loop ---\n
    segments_data: List[dict] = [] # List to store {'t_start':..., 't_end':..., 'interpolator':...}
    crossings_y_at_x_targets: List[List[float]] = [[] for _ in range(num_x_targets)]
    current_t, current_y = 0.0, [r_0, phi_0]
    final_K = 0.0 # Track the maximum K (affine parameter) reached

    # --- Integration Loop ---\n
    while current_t < t_end:
        sol = solve_ivp(
            ODE_Solver,
            (current_t, t_end),
            current_y,
            events=active_events,
            dense_output=True, # MUST BE TRUE for accurate event interpolation and final trajectory interpolation
            rtol=1e-10, atol=1e-10,
        )

        # Store the dense output function and its K range for this segment
        if sol.t.size > 1: # Need at least two points for a valid interpolation interval
             segments_data.append({
                 't_start': sol.t[0],
                 't_end': sol.t[-1],
                 'interpolator': sol.sol # Store the callable interpolator object
             })
             final_K = sol.t[-1] # Update the final K reached in this integration attempt
        elif sol.t.size == 1: # Handle case where solver stops immediately (e.g., event at t=0 for this segment)
             final_K = sol.t[0]
             # This single point doesn't form an interpolation segment, handled in post-processing if it's the only point

        # --- Process NON-TERMINAL x_target events using INTERPOLATION ---\n
        if event_x_targets_active and sol.sol:
            for i, event_idx in enumerate(x_target_event_indices):
                # Check if the event index is valid and if any events of this type occurred
                if event_idx < len(sol.t_events) and sol.t_events[event_idx].size > 0:
                    target_list_index = i # Corresponds to the index in x_targets_list
                    for t_event in sol.t_events[event_idx]: # Iterate through all occurrences of this event
                        try:
                            # Interpolate state (r, phi) at the precise event time K_event
                            interpolated_u = sol.sol(t_event)
                            r_event, phi_event = interpolated_u
                            # Ensure the event is physical before recording
                            if r_event < 2 * M or r_event < 1e-10: continue
                            y_event = r_event * np.sin(phi_event)
                            crossings_y_at_x_targets[target_list_index].append(y_event)
                        except ValueError:
                            # Handle cases where t_event might be slightly outside sol.t range due to numerics
                            print(f"Warning: Could not interpolate at x_target event K={t_event}. Skipping point.")

        # --- Handle ALL TERMINAL events ---\n
        terminal_event_occurred = False
        event_time = t_end # Default if no event occurs before t_end
        event_state = sol.y[:,-1] if sol.t.size > 0 else np.array(current_y)

        min_terminal_event_time = t_end
        triggering_event_index = -1 # Index in active_events list
        triggering_event_subindex = -1 # Index within sol.y_events[triggering_event_index]

        # Define which events are terminal for this logic block
        terminal_indices = [0, 1, 2] # Indices for event_fr_zero, event_r_leq_2M, event_r_max
        if event_x_stop_active: terminal_indices.append(event_x_stop_index)

        # Find the earliest terminal event that occurred in this segment
        for term_idx in terminal_indices:
             if term_idx < len(sol.t_events) and sol.t_events[term_idx].size > 0:
                 first_event_time_for_this_type = sol.t_events[term_idx][0]
                 if first_event_time_for_this_type <= min_terminal_event_time:
                    # If times are very close, prioritize critical events (fr=0, r=2M)
                    if abs(first_event_time_for_this_type - min_terminal_event_time) < 1e-12:
                         if term_idx == 0: # fr=0 is a high priority turning point
                            min_terminal_event_time = first_event_time_for_this_type
                            triggering_event_index = term_idx
                            triggering_event_subindex = 0
                         elif term_idx == 1 and triggering_event_index != 0: # r=2M, but fr=0 takes precedence if simultaneous
                            min_terminal_event_time = first_event_time_for_this_type
                            triggering_event_index = term_idx
                            triggering_event_subindex = 0
                    else: # This event is strictly earlier
                        min_terminal_event_time = first_event_time_for_this_type
                        triggering_event_index = term_idx
                        triggering_event_subindex = 0

        if triggering_event_index != -1: # A terminal event was identified
            terminal_event_occurred = True
            event_time = min_terminal_event_time
            try:
                # Use dense output to get state at the precise event time K_event
                event_state = sol.sol(event_time)
            except ValueError:
                 # Fallback if interpolation fails (e.g., event exactly at segment boundary)
                 print(f"Warning: Interpolation failed for terminal event {triggering_event_index} at K={event_time}. Using endpoint/event state.")
                 # Use the state recorded by solve_ivp for the event if available
                 if triggering_event_subindex != -1 and \
                    triggering_event_index < len(sol.y_events) and \
                    sol.y_events[triggering_event_index].ndim > 1 and \
                    triggering_event_subindex < sol.y_events[triggering_event_index].shape[0]:
                     event_state = sol.y_events[triggering_event_index][triggering_event_subindex].copy()
                 elif sol.t.size > 0: # Fallback to last point of trajectory segment
                     event_state = sol.y[:,-1].copy()
                 else: # Fallback to current_y if segment was empty
                     event_state = np.array(current_y).copy()

            if triggering_event_index == 0: # fr=0: Radial turning point, reverse sign and continue
                sign[0] *= -1
                current_t = event_time
                current_y = event_state.copy()

                # Nudge slightly past the turning point to avoid immediate re-triggering or numerical issues
                nudge_r_amount = 1e-7 * sign[0] # Nudge in the new direction of dr/dK
                current_y[0] = current_y[0] + nudge_r_amount
                current_t += 1e-7 # Also nudge affine parameter K slightly

                # Check if this nudge immediately caused another termination condition
                if current_y[0] <= 2 * M or current_y[0] >= r_max: break
                if event_x_stop_active and x_stop is not None and (current_y[0] * np.cos(current_y[1]) <= x_stop): break

                continue # Continue to the next integration segment

            else: # Any other terminal event (r_leq_2M, r_max, x_stop): Stop integration
                final_K = event_time # Update final_K to the exact event time
                break # Exit while loop

        # Handle solver status if no specific terminal event caused termination of this segment
        if not terminal_event_occurred:
            if sol.status == 1: # Reached t_end (the end of the current integration interval)
                final_K = t_end # The overall t_end for the whole trajectory was reached
                break
            elif sol.status < 0: # Solver failed
                print(f"Warning: Solver failed with status {sol.status} at K={current_t}, state={current_y}")
                final_K = sol.t[-1] if sol.t.size > 0 else current_t
                break
            elif sol.status == 0 and sol.t[-1] >= t_end: # Normal completion at t_end
                 final_K = sol.t[-1]
                 break
            # If sol.status == 0 and sol.t[-1] < t_end, it means a step was completed but no event
            # and not at t_end yet. This should be rare if events are set up correctly or t_eval is used.
            # The loop condition `current_t < t_end` and updating `current_t = sol.t[-1]` would handle this.
            # However, solve_ivp usually stops at an event or t_end[1] of the interval.

    # --- Post-Integration Interpolation ---\n
    if not segments_data and final_K == 0.0 and r_0 == current_y[0] and phi_0 == current_y[1]:
        # This handles cases like immediate stop (e.g. r_0 = r_max from start)
        # or if the initial conditions themselves were at an event.
        K_out = np.array([0.0])
        r_out = np.array([r_0])
        phi_out = np.array([phi_0 % (2 * np.pi)])
        x_out = r_out * np.cos(phi_out)
        y_out = r_out * np.sin(phi_out)
        return x_out, y_out, r_out, phi_out, K_out, crossings_y_at_x_targets

    # Create the unified K grid for interpolation, ensuring it spans 0 to final_K
    K_interp = np.linspace(0, final_K, num_interp_points)

    # Initialize output arrays for interpolated r and phi
    r_interp = np.zeros_like(K_interp)
    phi_interp = np.zeros_like(K_interp)

    # Perform interpolation across all collected segments
    current_segment_idx = 0
    for i, k_val in enumerate(K_interp):
        # Find the correct segment for this k_val
        # Start search from the last used segment index for efficiency
        found_segment = False
        for seg_idx in range(current_segment_idx, len(segments_data)):
            segment = segments_data[seg_idx]
            # Check if k_val is within the segment's K range (inclusive)
            # Add small tolerance for floating point comparisons at boundaries
            if segment['t_start'] - 1e-12 <= k_val <= segment['t_end'] + 1e-12:
                try:
                    r_interp[i], phi_interp[i] = segment['interpolator'](k_val)
                    current_segment_idx = seg_idx # Update hint for next search
                    found_segment = True
                    break # Move to the next k_val
                except ValueError:
                     # This might happen if k_val is *extremely* close to a boundary
                     # Try using the endpoint value as a fallback if interpolation fails at exact boundary
                     if abs(k_val - segment['t_start']) < 1e-10:
                         r_interp[i], phi_interp[i] = segment['interpolator'](segment['t_start'])
                         current_segment_idx = seg_idx
                         found_segment = True
                         break
                     elif abs(k_val - segment['t_end']) < 1e-10:
                          r_interp[i], phi_interp[i] = segment['interpolator'](segment['t_end'])
                          current_segment_idx = seg_idx
                          found_segment = True
                          break
                     else:
                         print(f"Warning: Interpolation failed for K={k_val} within segment {seg_idx} [{segment['t_start']},{segment['t_end']}]. Setting to NaN.")
                         r_interp[i], phi_interp[i] = np.nan, np.nan
                         found_segment = True # Mark as found to avoid fallback message later
                         break

        if not found_segment:
             # This case should ideally not happen if segments cover 0 to final_K.
             # It might occur if the very first point (k_val=0) wasn't captured in a segment (e.g., if integration started >0).
             if i == 0 and abs(k_val - 0.0) < 1e-12:
                 r_interp[i], phi_interp[i] = r_0, phi_0 # Use initial condition for K=0
             elif i == num_interp_points -1 and abs(k_val - final_K) < 1e-12 and segments_data:
                 # If it's the last point and very close to final_K, use the end of the last segment
                 r_interp[i], phi_interp[i] = segments_data[-1]['interpolator'](segments_data[-1]['t_end'])
             else:
                 print(f"Warning: Could not find segment for K={k_val} (final_K={final_K}). Setting state to NaN.")
                 r_interp[i], phi_interp[i] = np.nan, np.nan

    # Final post-processing of interpolated trajectory
    phi_interp = phi_interp % (2 * np.pi) # Ensure phi is in [0, 2*pi)
    x_interp = r_interp * np.cos(phi_interp)
    y_interp = r_interp * np.sin(phi_interp)

    # Remove NaN values if any occurred (indicates a problem in segment finding or interpolation)
    nan_mask = np.isnan(r_interp) | np.isnan(phi_interp)
    if np.any(nan_mask):
        print(f"Warning: Removing {np.sum(nan_mask)} NaN values from interpolated output.")
        K_out = K_interp[~nan_mask]
        x_out = x_interp[~nan_mask]
        y_out = y_interp[~nan_mask]
        r_out = r_interp[~nan_mask]
        phi_out = phi_interp[~nan_mask]
    else:
        K_out, x_out, y_out, r_out, phi_out = K_interp, x_interp, y_interp, r_interp, phi_interp

    return x_out, y_out, r_out, phi_out, K_out, crossings_y_at_x_targets

# compute_trajectory_crossings_only: A slimmed-down version of compute_trajectory,
# this function only returns the y-crossings at specified x_targets.
from typing import List, Optional, Tuple, Callable # Already imported, but good for explicitness here

def compute_trajectory_crossings_only(r_0: float, phi_0: float, M: float, psi: float, r_max: float, x_stop: Optional[float] = None, x_targets: Optional[List[float]] = None, t_end: float = 10**6) -> List[List[float]]:
    """Calculate trajectory intersections with x_targets accurately.

    Based on the robust logic of `compute_trajectory`, but optimized to return
    *only* the `crossings_y_at_x_targets` list for improved speed when the full
    trajectory is not needed. Uses interpolation for accuracy at crossing points.

    :param r_0: Initial radial coordinate.
    :param phi_0: Initial azimuthal coordinate (angle).
    :param M: Mass of the central object (e.g., black hole).
    :param psi: Initial angle related to the photon's direction, affecting impact parameter.
    :param r_max: Maximum radial distance before integration stops (escape boundary).
    :param x_stop: If provided (number), integration stops if the Cartesian x-coordinate
                   decreases below this value.
    :param x_targets: A list/array of Cartesian x-values. Records y-values when the trajectory
                      crosses each x-value (specifically, when x is decreasing past the target).
                      Uses interpolation for accuracy.
    :param t_end: Maximum integration parameter (affine parameter K) value.
    :return: A list of lists. Each inner list corresponds to an x_target and
             contains the y-values recorded at each crossing of that target.
             Example: [[y_cross1_xtarget1, y_cross2_xtarget1], [y_cross1_xtarget2], ...]
             Returns empty lists ([[], [], ...]) if x_targets is None or empty,
             or if the trajectory cannot be computed (e.g., r_0 <= 2M).

    DocTests:
    >>> # Test simple radial trajectory (b=0) crossing specific x (r) targets
    >>> # M=1, r_0=5, inward radial (psi=np.pi), x_targets=[3.0, 2.5]
    >>> # Since phi_0=0 and d(phi)/dK = 0 for b=0, y should be 0 at crossings.
    >>> crossings = compute_trajectory_crossings_only(r_0=5.0, phi_0=0.0, M=1.0, psi=np.pi, r_max=10.0, x_targets=[3.0, 2.5], x_stop=1.0)
    >>> len(crossings)
    2
    >>> len(crossings[0]) # Crossings for x_target = 3.0
    1
    >>> np.isclose(crossings[0][0], 0.0)
    True
    >>> len(crossings[1]) # Crossings for x_target = 2.5
    1
    >>> np.isclose(crossings[1][0], 0.0)
    True
    >>> # Test case where no targets are provided
    >>> crossings_none = compute_trajectory_crossings_only(r_0=5.0, phi_0=0.0, M=1.0, psi=0.0, r_max=10.0)
    >>> crossings_none
    []
    >>> # Test case starting inside event horizon
    >>> crossings_inside = compute_trajectory_crossings_only(r_0=1.0, phi_0=0.0, M=1.0, psi=0.0, r_max=10.0, x_targets=[0.5])
    >>> crossings_inside
    [[]]
    """
    # --- Basic Setup ---\n
    # Handle cases where x_targets is None or empty early
    num_x_targets_initial = 0
    x_targets_list_initial: List[float] = []
    if x_targets is not None:
         try:
             x_targets_list_initial = list(x_targets) # Check iterability and copy early
             num_x_targets_initial = len(x_targets_list_initial)
         except TypeError:
             print("Warning: x_targets provided but not iterable. Will return empty lists as if x_targets was empty.")
             # num_x_targets_initial remains 0, x_targets_list_initial remains []

    # Ensure r_0 is far enough for b calculation and initial state is physical
    if r_0 <= 2 * M: return [[] for _ in range(num_x_targets_initial)]
    metric_term_r0 = (1 - 2 * M / r_0)
    if metric_term_r0 <= 0: return [[] for _ in range(num_x_targets_initial)] # Should be caught by r_0 <= 2M
    b = np.cos(psi) * r_0 / np.sqrt(metric_term_r0)
    sign = [1] # Mutable list for pass-by-reference behavior in ODE_Solver
    if psi % (2 * np.pi) > np.pi: sign[0] *= -1

    # --- ODE Definition ---\n
    def ODE_Solver(K: float, u: np.ndarray) -> List[float]:
        r, phi = u
        if r <= 2 * M: return [0.0, 0.0] # Captured
        metric_term = (1 - 2 * M / r)
        term_b_r_sq = (b / r) ** 2
        fr = 1 - metric_term * term_b_r_sq # Related to radial potential
        drdK = sign[0] * np.sqrt(abs(fr)) # Use abs() for robustness near turning points (fr=0)
        if r < 1e-15: dphidK = 0.0 # Avoid division by zero if r is extremely small (though unlikely if r > 2M)
        else: dphidK = b / r**2
        return [drdK, dphidK]

    # --- Standard Terminal Event Definitions ---\n
    def event_fr_zero(K: float, u: np.ndarray) -> float: # Index 0: Radial turning point
        r, phi = u
        if r <= (2 * M + 1e-12): return 0.0 # Robustness near horizon
        metric_term = (1 - 2 * M / r)
        term_b_r_sq = (b / r) ** 2
        return 1 - metric_term * term_b_r_sq
    event_fr_zero.terminal = True
    event_fr_zero.direction = 0

    def event_r_leq_2M(K: float, u: np.ndarray) -> float: # Index 1: Reaches event horizon
        r, phi = u
        return r - (2 * M + 1e-12)
    event_r_leq_2M.terminal = True
    event_r_leq_2M.direction = -1

    def event_r_max(K: float, u: np.ndarray) -> float: # Index 2: Reaches outer boundary
        r, phi = u
        return r - r_max
    event_r_max.terminal = True
    event_r_max.direction = 1

    # --- Initialize Event List & Indices ---\n
    active_events: List[Callable] = [event_fr_zero, event_r_leq_2M, event_r_max]
    event_x_stop_index = -1
    x_target_event_indices: List[int] = []
    event_x_stop_active = isinstance(x_stop, (int, float))
    if event_x_stop_active:
        def event_x_stop_func(K: float, u: np.ndarray) -> float:
            r, phi = u; x = r * np.cos(phi)
            if r < 1e-10: return 1.0 # Avoid issues if r is very small
            return x - x_stop # type: ignore[operator]
        event_x_stop_func.terminal = True; event_x_stop_func.direction = -1
        event_x_stop_index = len(active_events); active_events.append(event_x_stop_func)

    num_x_targets = 0; event_x_targets_active = False; x_targets_list: List[float] = []
    if x_targets is not None and num_x_targets_initial > 0:
        x_targets_list = x_targets_list_initial # Use the validated and copied list
        num_x_targets = num_x_targets_initial; event_x_targets_active = True

    if event_x_targets_active:
        base_idx = len(active_events)
        for i, xt_val in enumerate(x_targets_list):
            def make_event_func(target_x: float) -> Callable:
                def event_x_target_dynamic(K: float, u: np.ndarray) -> float:
                    r, phi = u; x = r * np.cos(phi)
                    if r < 1e-10: return 1.0 # Avoid issues if r is very small
                    return x - target_x
                event_x_target_dynamic.terminal = False; event_x_target_dynamic.direction = -1
                return event_x_target_dynamic
            event_func = make_event_func(xt_val)
            active_events.append(event_func); x_target_event_indices.append(base_idx + i)

    # --- Initialization for Integration Loop ---\n
    crossings_y_at_x_targets: List[List[float]] = [[] for _ in range(num_x_targets)]
    current_t, current_y = 0.0, [r_0, phi_0]

    # --- Integration Loop ---\n
    while current_t < t_end:
        sol = solve_ivp(
            ODE_Solver,
            (current_t, t_end),
            current_y,
            events=active_events,
            dense_output=True, # KEEP True for accurate event interpolation
            rtol=1e-8, atol=1e-10, # Tolerances can be adjusted; may differ from full trajectory func
        )

        # --- Process NON-TERMINAL x_target events using INTERPOLATION ---\n
        if event_x_targets_active and sol.sol is not None:
            for i, event_idx in enumerate(x_target_event_indices):
                if event_idx < len(sol.t_events) and sol.t_events[event_idx].size > 0:
                    target_list_index = i
                    for t_event in sol.t_events[event_idx]:
                        # Check event time is within the solved segment for robust interpolation
                        if not (sol.t[0] <= t_event <= sol.t[-1]): continue
                        try:
                            interpolated_u = sol.sol(t_event); r_event, phi_event = interpolated_u
                            if r_event < 2 * M or r_event < 1e-10: continue # Physical check
                            y_event = r_event * np.sin(phi_event)
                            crossings_y_at_x_targets[target_list_index].append(y_event)
                        except ValueError:
                            print(f"Warning: Could not interpolate state at x_target event K={t_event}. Skipping.")

        # --- Handle ALL TERMINAL events ---\n
        terminal_event_occurred = False
        event_time = t_end # Default if no event before t_end of segment
        event_state = sol.y[:,-1] if sol.t.size > 0 else np.array(current_y) # Default state

        min_terminal_event_time = t_end
        triggering_event_index = -1
        triggering_event_subindex = -1 # For potential fallback state from sol.y_events

        terminal_indices = [0, 1, 2] # Indices for event_fr_zero, event_r_leq_2M, event_r_max
        if event_x_stop_active: terminal_indices.append(event_x_stop_index)

        # Find the earliest terminal event in this segment
        for term_idx in terminal_indices:
             if term_idx < len(sol.t_events) and sol.t_events[term_idx].size > 0:
                 first_event_time_for_this_type = sol.t_events[term_idx][0]
                 # Check if this event is strictly earlier (with tolerance)
                 if first_event_time_for_this_type < min_terminal_event_time - 1e-12:
                     min_terminal_event_time = first_event_time_for_this_type
                     triggering_event_index = term_idx
                     triggering_event_subindex = 0
                 elif abs(first_event_time_for_this_type - min_terminal_event_time) < 1e-12:
                     # If times are near-identical, prioritize critical events (fr=0, then r=2M)
                     if term_idx == 0 and triggering_event_index not in [0]: # fr=0 highest priority
                         min_terminal_event_time = first_event_time_for_this_type
                         triggering_event_index = term_idx
                         triggering_event_subindex = 0
                     elif term_idx == 1 and triggering_event_index not in [0, 1]: # r=2M next priority
                         min_terminal_event_time = first_event_time_for_this_type
                         triggering_event_index = term_idx
                         triggering_event_subindex = 0
                     # else: keep the already found event if it was of higher or equal priority

        # Process the first terminal event that occurred, if any, within this segment's K range
        if triggering_event_index != -1 and min_terminal_event_time <= sol.t[-1] + 1e-12 :
            terminal_event_occurred = True
            event_time = min_terminal_event_time
            try:
                # Use dense output for accuracy, ensuring K is within segment bounds for interpolation
                t_interp = np.clip(event_time, sol.t[0], sol.t[-1])
                event_state = sol.sol(t_interp)
            except ValueError:
                 print(f"Warning: Interpolation failed for terminal event {triggering_event_index} at K={event_time}. Using endpoint/event state.")
                 # Fallback: try using y_events if available and index is valid
                 if triggering_event_subindex != -1 and \
                    triggering_event_index < len(sol.y_events) and \
                    sol.y_events[triggering_event_index].ndim > 1 and \
                    triggering_event_subindex < sol.y_events[triggering_event_index].shape[0]:
                     event_state = sol.y_events[triggering_event_index][triggering_event_subindex].copy()
                 elif sol.t.size > 0: # Fallback to last point of segment
                     event_state = sol.y[:,-1].copy()
                 else: # Fallback to current_y if segment was empty
                     event_state = np.array(current_y).copy()

            if triggering_event_index == 0: # fr=0: Radial turning point, reverse sign and continue
                sign[0] *= -1
                current_t = event_time
                current_y = event_state.copy()
                # Nudge slightly past turning point
                nudge_r_amount = 1e-7 * sign[0]
                current_y[0] = current_y[0] + nudge_r_amount
                current_t += 1e-7 # Also nudge K slightly
                # Check if nudge immediately caused another termination
                if current_y[0] <= 2 * M or current_y[0] >= r_max: break
                if event_x_stop_active and x_stop is not None and (current_y[0] * np.cos(current_y[1]) <= x_stop): break
                continue # Continue to next integration segment

            else: # Any other terminal event: Stop integration
                break # Exit while loop

        # --- Handle solver status if no specific terminal event caused termination ---\n
        if not terminal_event_occurred:
            if sol.status == 1: # Reached t_end for this segment (sol.t[-1])
                 if sol.t[-1] >= t_end - 1e-12: break # Reached overall t_end
                 else: # Should not happen if interval is (current_t, t_end) unless solver stopped early
                      current_t = sol.t[-1]; current_y = sol.y[:, -1].copy()
            elif sol.status < 0: # Solver failed
                print(f"Warning: Solver failed with status {sol.status} at K={current_t}, state={current_y}")
                break
            elif sol.status == 0: # Step completed successfully
                 if sol.t[-1] >= t_end - 1e-12: break # Reached overall t_end exactly
                 else: # Step finished before t_end, no event? Update and continue.
                      current_t = sol.t[-1]; current_y = sol.y[:, -1].copy()

    # --- Return ONLY the accumulated y-values at x-target crossings ---\n
    # No final trajectory interpolation is needed for this specialized function.
    return crossings_y_at_x_targets

# Image_map: Cartesian input for initial photon direction parameters.
from typing import Union, Tuple, List # Already imported

def Image_map(Y: float, Z: float, x_0: float, x_1: float, x_2: float, a: float) -> Union[Tuple[float, float], str]:
    """Map a point on an observer's window to intersection points on source/window planes.

    Takes initial parameters (Y,Z) defining a point on a plane (related to the
    observer's window at x=x_1) and traces a photon backwards from an observer
    at x=x_2. It calculates where this photon path intersects a source plane at x=x_0
    and the window plane at x=x_1. Coordinates (Y,Z) are Cartesian-like on the
    plane perpendicular to the x-axis.

    :param Y: Y-coordinate on the initial parameterization plane (related to window).
    :param Z: Z-coordinate on the initial parameterization plane (related to window).
    :param x_0: X-coordinate of the background source plane.
    :param x_1: X-coordinate of the observer's window plane.
    :param x_2: X-coordinate of the observer (focal point).
    :param a: A scaling factor, potentially related to the size of the window or scene.
    :return: A tuple (y_window, y_image) representing the y-coordinates (in the
             Schwarzschild 2D plane of calculation) where the photon crosses the window
             (at x_1) and image/source plane (at x_0) respectively. Returns "Miss" if
             either crossing is not found.

    DocTests:
    >>> # Test a radial shot (Y=0, Z=0 so b=0) from x_2=100 towards origin (M=1).
    >>> # Observer at x_2=20, window at x_1=10, source at x_0=3. 'a' is scale.
    >>> # For radial, y_window and y_image should be 0.
    >>> result = Image_map(Y=0.0, Z=0.0, x_0=3.0, x_1=10.0, x_2=20.0, a=1.0)
    >>> isinstance(result, tuple) and len(result) == 2 and np.isclose(result[0], 0.0) and np.isclose(result[1], 0.0)
    True
    >>> # Test a case expected to miss (e.g., starting too far inside event horizon effectively for psi calc)
    >>> # (Actual miss depends on compute_trajectory_crossings_only logic for invalid r_0 for b)
    >>> # If x_2 (r_0 for trajectory) is inside 2M, compute_trajectory_crossings_only returns empty.
    >>> result_miss = Image_map(Y=0.0, Z=0.0, x_0=-10.0, x_1=5.0, x_2=1.0, a=1.0) # x_2=1, M=1
    >>> result_miss
    'Miss'
    """
    # Effective size parameter for r_max calculation, related to 'a' and observer/source distances.
    a_prime = (x_2 - x_0) / (x_2 - x_1) * a if (x_2 - x_1) != 0 else a 
    # Maximum radius for integration, heuristically set based on scene geometry.
    r_max = np.sqrt(2 * a_prime**2 + x_0**2) 

    # Calculate initial psi angle for the photon trajectory.
    # This psi is the angle relative to the local radial direction in the Schwarzschild (r,phi) plane.
    # It's derived from the geometry of the observer at x_2 looking towards (Y,Z) on the plane x_1.
    rho_sq = Y**2 + Z**2
    observer_to_window_dist_sq = (x_2 - x_1)**2 + rho_sq
    if observer_to_window_dist_sq <= 0: # Avoid sqrt domain error or division by zero
        return "Miss" # Geometrically problematic setup
    
    # cos_angle_with_normal = (x_2-x_1) / sqrt(observer_to_window_dist_sq)
    # sin_angle_with_normal = sqrt(rho_sq) / sqrt(observer_to_window_dist_sq)
    # psi is angle with local r-direction. If Y,Z are small, ray is mostly along -x (psi ~ pi or 0 depending on convention)
    # The formula used here defines psi relative to the direction perpendicular to the line of sight to the origin.
    psi = -np.arccos(np.sqrt(rho_sq) / np.sqrt(observer_to_window_dist_sq))

    # Initial conditions for compute_trajectory: r_0 is observer's x-position, phi_0 is 0 by convention.
    # M=1 sets the scale (mass of black hole).
    # x_targets are the x-planes for which we want y-crossings.
    # x_stop is a safety boundary. 
    y_list = compute_trajectory_crossings_only(r_0=x_2, phi_0=0.0, M=1.0, psi=psi, r_max=r_max, x_targets=[x_0, x_1], x_stop=x_0 - 1)
    
    # y_list[0] corresponds to crossings of x_targets[0] (i.e., x_0, the source plane)
    # y_list[1] corresponds to crossings of x_targets[1] (i.e., x_1, the window plane)
    if not y_list or len(y_list) < 2: # Ensure y_list has expected structure
        return "Miss"

    y_image_crossings = y_list[0]
    y_window_crossings = y_list[1]
    
    if len(y_image_crossings) == 0 or len(y_window_crossings) == 0:
        return "Miss"
    
    # Return the first recorded crossing for each target plane.
    return y_window_crossings[0], y_image_crossings[0]
    
# Image_map_radial: Radial input 'r' for initial photon direction parameters.
from typing import Union, Tuple, List # Already imported

def Image_map_radial(r: float, x_0: float, x_1: float, x_2: float, a: float) -> Union[Tuple[float, float], str]:
    """Map a radial distance on an observer's window to intersection points.

    Similar to Image_map, but takes an initial radial parameter 'r' (rho)
    on a plane (related to the observer's window at x=x_1) to define the
    photon's initial direction from an observer at x=x_2. Traces this photon
    backwards to find intersections with a source plane at x=x_0 and the
    window plane at x=x_1.

    :param r: Radial distance on the initial parameterization plane (related to window).
    :param x_0: X-coordinate of the background source plane.
    :param x_1: X-coordinate of the observer's window plane.
    :param x_2: X-coordinate of the observer (focal point).
    :param a: A scaling factor, potentially related to the size of the window or scene.
    :return: A tuple (y_window, y_image) representing the y-coordinates (in the
             Schwarzschild 2D plane of calculation) where the photon crosses the window
             (at x_1) and image/source plane (at x_0) respectively. Returns "Miss" if
             either crossing is not found.
    DocTests:
    >>> # Test a radial shot (r=0 so b=0) from x_2=100 towards origin (M=1).
    >>> # Observer at x_2=20, window at x_1=10, source at x_0=3. 'a' is scale.
    >>> # For radial, y_window and y_image should be 0.
    >>> result = Image_map_radial(r=0.0, x_0=3.0, x_1=10.0, x_2=20.0, a=1.0)
    >>> isinstance(result, tuple) and len(result) == 2 and np.isclose(result[0], 0.0) and np.isclose(result[1], 0.0)
    True
    >>> # Test a case expected to miss (e.g., x_2 (r_0 for trajectory) inside 2M)
    >>> result_miss = Image_map_radial(r=0.0, x_0=-10.0, x_1=5.0, x_2=1.0, a=1.0) # x_2=1, M=1
    >>> result_miss
    'Miss'
    """
    # Effective size parameter for r_max calculation.
    a_prime = (x_2 - x_0) / (x_2 - x_1) * a if (x_2 - x_1) != 0 else a
    # Maximum radius for integration.
    r_max = np.sqrt(2 * a_prime**2 + x_0**2)

    # Calculate initial psi angle for the photon trajectory using radial distance 'r'.
    # 'r' here is equivalent to rho = sqrt(Y^2+Z^2) in Image_map.
    r_sq = r**2
    observer_to_window_dist_sq = (x_2 - x_1)**2 + r_sq
    if observer_to_window_dist_sq <= 0:
        return "Miss"
    
    psi = -np.arccos(np.sqrt(r_sq) / np.sqrt(observer_to_window_dist_sq))
    
    # Call compute_trajectory_crossings_only with M=1 by convention for these mapping functions.
    y_list = compute_trajectory_crossings_only(r_0=x_2, phi_0=0.0, M=1.0, psi=psi, r_max=r_max, x_targets=[x_0, x_1], x_stop=x_0 - 1)
    
    if not y_list or len(y_list) < 2:
        return "Miss"

    y_image_crossings = y_list[0]
    y_window_crossings = y_list[1]
    
    if len(y_image_crossings) == 0 or len(y_window_crossings) == 0:
        return "Miss"
    
    return y_window_crossings[0], y_image_crossings[0]

# Results_Cartesian: Generates and saves photon mapping data using Cartesian grid sampling.
import numpy as np
import time
import os # Added for path operations
from typing import List, Tuple # Already imported, but good for explicitness

# NOTE: This function assumes that a function named 'Image_map'
#       is defined or imported in the scope where Results_Cartesian is called.

def Results_Cartesian(x_0: float, x_1: float, x_2: float, a: float, n: int, chunk_size: float = 1e7) -> List[str]:
    """Generate and save photon mapping results by sampling a Cartesian grid.

    Iterates over a grid defined by (Y,Z) coordinates (up to 'a' with 'n' steps),
    calls 'Image_map' for each point, and stores the results (transformed to
    Cartesian coordinates using multiple phi_prime angles) in chunks as .npy files.

    :param x_0: X-coordinate of the background source plane.
    :param x_1: X-coordinate of the observer's window plane.
    :param x_2: X-coordinate of the observer (focal point).
    :param a: Maximum value for Y and Z coordinates in the sampling grid.
    :param n: Number of divisions for Y and Z up to 'a', defining grid resolution.
    :param chunk_size: Maximum number of results to store in memory before saving to a file.
    :return: A list of file paths to the saved .npy chunk files.
    :raises IOError: If saving a chunk to a file fails.

    DocTests:
    >>> # Test with minimal n to check file creation and basic loop structure
    >>> # This will be slow if Image_map is actually running full computations.
    >>> # For a real doctest, Image_map might need to be mocked or a very fast scenario used.
    >>> # As a placeholder, we'll assume it runs and check return type. 
    >>> # Note: Actual file creation is a side effect not easily tested in doctest without cleanup.
    >>> # saved_files = Results_Cartesian(x_0=-10.0, x_1=5.0, x_2=20.0, a=0.1, n=1, chunk_size=10)
    >>> # isinstance(saved_files, list) # Expected True
    >>> # This is too slow for a doctest. Example of a call:
    >>> # Results_Cartesian(x_0=-20, x_1=10, x_2=100, a=1, n=2, chunk_size=5)
    """

    start_time = time.time()

    results_chunk: List[Tuple[float, float, float, float]] = [] 
    saved_files: List[str] = []
    chunk_index = 0
    chunk_size = int(chunk_size) # Ensure chunk_size is an integer

    print("Starting Results_Cartesian generation...")

    # --- Define base filename for saved chunks ---\n
    x0_str = str(x_0).replace('.','p').replace('-','neg')
    x1_str = str(x_1).replace('.','p').replace('-','neg')
    x2_str = str(x_2).replace('.','p').replace('-','neg')
    n_str = str(n)
    base_save_path = f'Results_Light_Ring_Cartesian_{x0_str}_{x1_str}_{x2_str}_{n_str}'
    output_directory = "." # Save in current directory
    os.makedirs(output_directory, exist_ok=True)

    # --- Helper function to save a chunk of results ---\n
    def save_chunk(current_chunk_list: List[Tuple[float, float, float, float]], current_chunk_index: int) -> Optional[str]:
        if not current_chunk_list: return None
        chunk_save_path = os.path.join(output_directory, f"{base_save_path}_chunk_{current_chunk_index:03d}.npy")
        print(f"\nSaving chunk {current_chunk_index} ({len(current_chunk_list)} results) to {chunk_save_path}...")
        try:
            np.save(chunk_save_path, np.array(current_chunk_list, dtype=np.float64))
            print(f"Chunk {current_chunk_index} saved.")
            return chunk_save_path
        except Exception as e:
            print(f"\nERROR saving chunk {current_chunk_index} to {chunk_save_path}: {e}")
            return None

    # --- Main Loops for Cartesian grid sampling ---\n
    try:
        for i in range(n + 1):
            y = a * i / n if n > 0 else (a if i > 0 else 0) # Y-coordinate on sampling grid

            for j in range(i + 1): # Exploit symmetry: Z <= Y
                z = a * j / n if n > 0 else (a if j > 0 else 0) # Z-coordinate on sampling grid

                if j % 10**4 == 0: # Progress indicator
                    print(f"y={y}, z={z}")
                # Calculate effective phi_prime for this (y,z) point
                phi_prime_base = np.arctan2(z, y) if not (y == 0 and z == 0) else 0.0

                # Call Image_map to get window and image plane y-crossings
                output = Image_map(Y=y, Z=z, x_0=x_0, x_1=x_1, x_2=x_2, a=a)
                if output == 'Miss':
                    continue
                else:
                    y_window, y_image = output # type: ignore
                
                # Define symmetry angles for phi_prime
                if y == z or z == 0: # Higher symmetry cases (e.g., on diagonal or axis)
                    phi_values = [phi_prime_base, phi_prime_base + np.pi/2, phi_prime_base + np.pi, phi_prime_base + 3/2 * np.pi]
                else: # General case with 8-fold symmetry
                    phi_values = [
                        phi_prime_base, phi_prime_base + np.pi/2, phi_prime_base + np.pi, phi_prime_base + 3/2 * np.pi,
                        np.pi/2 - phi_prime_base, np.pi - phi_prime_base, 3*np.pi/2 - phi_prime_base, 2*np.pi - phi_prime_base]
                
                for phi_prime_val in phi_values:
                    results_chunk.append((
                        y_window * np.cos(phi_prime_val), y_window * np.sin(phi_prime_val),
                        y_image * np.cos(phi_prime_val), y_image * np.sin(phi_prime_val)))

                # --- Check chunk size and save if needed ---\n
                if len(results_chunk) >= chunk_size:
                    saved_path = save_chunk(results_chunk, chunk_index)
                    if saved_path:
                        saved_files.append(saved_path)
                    else:
                        raise IOError(f"Failed to save chunk {chunk_index}")
                    results_chunk = [] # Reset chunk for next batch of results
                    chunk_index += 1

    except KeyboardInterrupt:
         print("\n--- Process Interrupted by User ---")
    except Exception as e:
         print(f"\n--- Error during processing: {e} ---")
         import traceback
         traceback.print_exc()
    finally:
        # --- Save the final (potentially incomplete) chunk ---\n
        print("\nAttempting to save final chunk...")
        final_saved_path = save_chunk(results_chunk, chunk_index)
        if final_saved_path:
            saved_files.append(final_saved_path)

        end_time = time.time()
        duration = end_time - start_time

        total_results_saved = 0
        for f_path in saved_files:
            if os.path.exists(f_path):
                try:
                    total_results_saved += len(np.load(f_path))
                except Exception as load_err:
                    print(f"Warning: Could not load file {f_path} to count results: {load_err}")

        print(f"\nFinished Results_Cartesian generation attempt.")
        print(f"Total results saved: {total_results_saved}")
        print(f"Number of chunk files created: {len(saved_files)}")
        print(f"Total execution time: {duration:.2f} seconds")

        return saved_files

# Results_Radial: Generates and saves photon mapping data using radial grid sampling.
import numpy as np
import time
import os # Added for path operations
from typing import List, Tuple, Optional # Already imported

# NOTE: This function assumes that a function named 'Image_map_radial'
#       is defined or imported in the scope where Results_Radial is called.

def Results_Radial(x_0: float, x_1: float, x_2: float, R_max: float, n: int, chunk_size: float = 1e7) -> List[str]:
    """Generate and save photon mapping results by sampling a radial grid.

    Iterates over radial distances 'r' (up to 'R_max' with 'n' steps) and angles
    'phi_prime'. For each 'r', it calls 'Image_map_radial' once, then generates
    multiple data points by rotating the result by 'phi_prime'. Stores results
    in chunks as .npy files.

    :param x_0: X-coordinate of the background source plane.
    :param x_1: X-coordinate of the observer's window plane.
    :param x_2: X-coordinate of the observer (focal point).
    :param R_max: Maximum radial distance 'r' for sampling on the window plane.
    :param n: Number of divisions for 'r' up to 'R_max', defining radial resolution.
    :param chunk_size: Maximum number of results to store in memory before saving to a file.
    :return: A list of file paths to the saved .npy chunk files.
    :raises IOError: If saving a chunk to a file fails.

    DocTests:
    >>> # Similar to Results_Cartesian, a full doctest is slow and involves side effects.
    >>> # Placeholder for call structure:
    >>> # Results_Radial(x_0=-20, x_1=10, x_2=100, R_max=0.5, n=2, chunk_size=10)
    """

    start_time = time.time()

    a = 1.0 * (x_2 - x_1) # Scaling factor 'a' for Image_map_radial, related to observer-window distance
    results_chunk: List[Tuple[float, float, float, float]] = []
    saved_files: List[str] = []
    chunk_index = 0
    chunk_size = int(chunk_size) # Ensure chunk_size is an integer

    print("Starting Results_Radial generation...")

    # --- Define base filename for saved chunks ---\n
    x0_str = str(x_0).replace('.','p').replace('-','neg')
    x1_str = str(x_1).replace('.','p').replace('-','neg')
    x2_str = str(x_2).replace('.','p').replace('-','neg')
    n_str = str(n)
    # Filename indicates this might be for a specific visualization (rainbow_thesis)
    base_save_path = f'Results_rainbow_thesis_{x0_str}_{x1_str}_{x2_str}_{n_str}'
    output_directory = "." # Save in current directory
    os.makedirs(output_directory, exist_ok=True)

    # --- Helper function to save a chunk of results ---\n
    def save_chunk(current_chunk_list: List[Tuple[float, float, float, float]], current_chunk_index: int) -> Optional[str]:
        if not current_chunk_list: return None
        chunk_save_path = os.path.join(output_directory, f"{base_save_path}_chunk_{current_chunk_index:03d}.npy")
        print(f"\nSaving chunk {current_chunk_index} ({len(current_chunk_list)} results) to {chunk_save_path}...")
        try:
            np.save(chunk_save_path, np.array(current_chunk_list, dtype=np.float64))
            print(f"Chunk {current_chunk_index} saved.")
            return chunk_save_path
        except Exception as e:
            print(f"\nERROR saving chunk {current_chunk_index} to {chunk_save_path}: {e}")
            return None

    # --- Main Loops for radial grid sampling ---\n
    try:
        for i in range(n + 1):
            r_sample = R_max * i / n if n > 0 else (R_max if i > 0 else 0) # Radial sample point

            # Number of angular steps 'k', increases with r_sample for denser sampling at larger radii
            k = max(int(5 * 10**4 * r_sample), int(1))

            # Call Image_map_radial once per r_sample, as the y_window/y_image depends only on magnitude of r_sample
            output = Image_map_radial(r_sample, x_0, x_1, x_2, a)
            valid_output = (output != "Miss")
            y_window: float; y_image: float # Declare types for mypy
            if valid_output:
                 y_window, y_image = output # type: ignore
            # If not valid_output, y_window and y_image are not assigned here, 
            # but the inner loop conditional `if valid_output:` handles this.

            for j in range(k + 1):
                # Calculate phi_prime for rotating the (y_window, y_image) result
                phi_prime_val = (j / k * 2 * np.pi) if k > 0 else 0.0

                if j % 10**4 == 0: # Progress indicator
                    print(f"r_sample={r_sample}, phi_prime_val={phi_prime_val}")

                if valid_output: # Only append if Image_map_radial returned a valid hit
                    cos_phi = np.cos(phi_prime_val)
                    sin_phi = np.sin(phi_prime_val)
                    results_chunk.append((y_window * cos_phi, y_window * sin_phi, y_image * cos_phi, y_image * sin_phi))

                # --- Check chunk size and save if needed ---\n
                if len(results_chunk) >= chunk_size:
                    saved_path = save_chunk(results_chunk, chunk_index)
                    if saved_path:
                        saved_files.append(saved_path)
                    else:
                        raise IOError(f"Failed to save chunk {chunk_index}")
                    results_chunk = [] # Reset chunk
                    chunk_index += 1

    except KeyboardInterrupt:
         print("\n--- Process Interrupted by User ---")
    except Exception as e:
         print(f"\n--- Error during processing: {e} ---")
         import traceback
         traceback.print_exc()
    finally:
        # --- Save the final (potentially incomplete) chunk ---\n
        print("\nAttempting to save final chunk...")
        final_saved_path = save_chunk(results_chunk, chunk_index)
        if final_saved_path:
            saved_files.append(final_saved_path)

        end_time = time.time()
        duration = end_time - start_time

        total_results_saved = 0
        for f_path in saved_files:
            if os.path.exists(f_path):
                try:
                    total_results_saved += len(np.load(f_path))
                except Exception as load_err:
                    print(f"Warning: Could not load file {f_path} to count results: {load_err}")

        print(f"\nFinished Results_Radial generation attempt.")
        print(f"Total results saved: {total_results_saved}")
        print(f"Number of chunk files created: {len(saved_files)}")
        print(f"Total execution time: {duration:.2f} seconds")

        return saved_files

# Results_Radial_Light_Ring: Generates data for light ring phenomena with specific radial/angular sampling.
import numpy as np
import time
import os # Added for path operations
from typing import List, Tuple, Optional # Already imported

# NOTE: This function assumes that a function named 'Image_map_radial'
#       is defined or imported in the scope where Results_Radial_Light_Ring is called.

def Results_Radial_Light_Ring(x_0: float, x_1: float, x_2: float, n: int, chunk_size: float = 1e7) -> Tuple[List[str], List[float], List[float]]:
    """Generate and save photon mapping results focused on light ring phenomena.

    This function uses specific, narrow ranges for radial ('r_sample') and angular
    ('phi_prime_val') sampling, likely to finely probe regions where light ring
    effects are prominent. It calls 'Image_map_radial' and stores results in chunks.
    It also attempts to track 'hits' and 'misses' based on 'Image_map_radial' output.

    :param x_0: X-coordinate of the background source plane.
    :param x_1: X-coordinate of the observer's window plane.
    :param x_2: X-coordinate of the observer (focal point).
    :param n: Number of divisions for the 'r_sample' range, defining radial resolution.
    :param chunk_size: Maximum number of results to store in memory before saving to a file.
    :return: A tuple containing: 
             - A list of file paths to the saved .npy chunk files.
             - A list of 'r_sample' values for which 'Image_map_radial' resulted in a 'Miss'.
             - A list of 'r_sample' values for which 'Image_map_radial' resulted in a 'Hit'.
    :raises IOError: If saving a chunk to a file fails.

    DocTests:
    >>> # Full doctest is slow. Example of a call structure:
    >>> # files, misses, hits = Results_Radial_Light_Ring(x_0=-20, x_1=10, x_2=100, n=2, chunk_size=10)
    """

    start_time = time.time()

    a = 1.0 * (x_2 - x_1) # Scaling factor 'a' for Image_map_radial
    results_chunk: List[Tuple[float, float, float, float]] = []
    saved_files: List[str] = []
    chunk_index = 0
    chunk_size = int(chunk_size)
    misses: List[float] = []
    hits: List[float] = []

    print("Starting Results_Radial_Light_Ring generation...")

    # --- Define base filename for saved chunks ---\n
    x0_str = str(x_0).replace('.','p').replace('-','neg')
    x1_str = str(x_1).replace('.','p').replace('-','neg')
    x2_str = str(x_2).replace('.','p').replace('-','neg')
    n_str = str(n)
    # Filename suggests focus on a light ring segment (0th order, large view?)
    base_save_path = f'Light_Ring_Segement_0th_large_{x0_str}_{x1_str}_{x2_str}_{n_str}'
    output_directory = "." # Save in current directory
    os.makedirs(output_directory, exist_ok=True)

    # --- Helper function to save a chunk of results ---\n
    def save_chunk(current_chunk_list: List[Tuple[float, float, float, float]], current_chunk_index: int) -> Optional[str]:
        if not current_chunk_list: return None
        chunk_save_path = os.path.join(output_directory, f"{base_save_path}_chunk_{current_chunk_index:03d}.npy")
        print(f"\nSaving chunk {current_chunk_index} ({len(current_chunk_list)} results) to {chunk_save_path}...")
        try:
            np.save(chunk_save_path, np.array(current_chunk_list, dtype=np.float64))
            print(f"Chunk {current_chunk_index} saved.")
            return chunk_save_path
        except Exception as e:
            print(f"\nERROR saving chunk {current_chunk_index} to {chunk_save_path}: {e}")
            return None

    # --- Main Loops for specialized radial and angular sampling ---\n
    try:
        for i in range(n + 1):
            # Specific radial sampling range, likely targeting light ring phenomena
            r_start = 0.2543706590950274
            r_delta = 0.0008610946154524735
            r_sample = r_start + (i / n * r_delta) if n > 0 else (r_start if i==0 else r_start + r_delta) 
            
            # Number of angular steps 'k'
            k = int(6 * 10**4 * r_sample)
            # k=0 # This was a commented-out test value
            # print(f"{k=}") # This was a commented-out debug print

            # Call Image_map_radial once per r_sample
            output = Image_map_radial(r_sample, x_0, x_1, x_2, a)
            valid_output = (output != "Miss")
            y_window: float; y_image: float # Declare types for mypy
            if valid_output:
                 y_window, y_image = output # type: ignore
                 hits.append(r_sample)
            else:
                misses.append(r_sample)

            for j in range(k + 1):
                # Specific angular sampling range for phi_prime, e.g. a narrow slice
                phi_prime_start_angle = 3/2 * np.pi
                phi_prime_angular_width = np.pi / 180 / 10 # Small angular width (0.1 degrees)
                phi_prime_val = phi_prime_start_angle + (j / k * phi_prime_angular_width) if k > 0 else phi_prime_start_angle

                if j % 10**4 == 0: # Progress indicator
                    print(f"r_sample={r_sample}, phi_prime_val={phi_prime_val}")

                if valid_output:
                    cos_phi = np.cos(phi_prime_val)
                    sin_phi = np.sin(phi_prime_val)
                    results_chunk.append((y_window * cos_phi, y_window * sin_phi, y_image * cos_phi, y_image * sin_phi))

                # --- Check chunk size and save if needed ---\n
                if len(results_chunk) >= chunk_size:
                    saved_path = save_chunk(results_chunk, chunk_index)
                    if saved_path:
                        saved_files.append(saved_path)
                    else:
                        raise IOError(f"Failed to save chunk {chunk_index}")
                    results_chunk = [] # Reset chunk
                    chunk_index += 1

    except KeyboardInterrupt:
         print("\n--- Process Interrupted by User ---")
    except Exception as e:
         print(f"\n--- Error during processing: {e} ---")
         import traceback
         traceback.print_exc()
    finally:
        # --- Save the final (potentially incomplete) chunk ---\n
        print("\nAttempting to save final chunk...")
        final_saved_path = save_chunk(results_chunk, chunk_index)
        if final_saved_path:
            saved_files.append(final_saved_path)

        end_time = time.time()
        duration = end_time - start_time

        total_results_saved = 0
        for f_path in saved_files:
            if os.path.exists(f_path):
                try:
                    total_results_saved += len(np.load(f_path))
                except Exception as load_err:
                    print(f"Warning: Could not load file {f_path} to count results: {load_err}")

        print(f"\nFinished Results_Radial_Light_Ring generation attempt.")
        print(f"Total results saved: {total_results_saved}")
        print(f"Number of chunk files created: {len(saved_files)}")
        print(f"Total execution time: {duration:.2f} seconds")

        return saved_files, misses, hits

# map_photons: Creates a lensed image from photon data and a source image.
import numpy as np
import time
import os
import sys # Already imported, but good for clarity
from PIL import Image, UnidentifiedImageError
import traceback
import math # Import math for ceiling function
from typing import Union, Tuple, List, Optional, Any # Already imported

def _save_image(image_array: np.ndarray, save_path: str) -> bool:
    """Helper function to save a NumPy array as an RGB image using PIL.

    Handles type conversion to uint8 and ensures the array is a valid RGB image.

    :param image_array: NumPy array representing the image data.
    :param save_path: File path where the image will be saved.
    :return: True if the image was saved successfully, False otherwise.
    """
    try:
        parent_dir = os.path.dirname(save_path)
        if parent_dir and not os.path.exists(parent_dir):
            os.makedirs(parent_dir, exist_ok=True)
        
        # Convert to uint8 if not already, handling different input types and ranges
        if image_array.dtype != np.uint8:
             if image_array.dtype in [np.float32, np.float64]:
                 # Ensure conversion handles potential NaN/Inf before clipping/casting
                 image_array_conv = np.nan_to_num(image_array, nan=0.0, posinf=255.0, neginf=0.0)
                 image_array_conv = np.clip(image_array_conv, 0, 255).astype(np.uint8)
             else:
                 try:
                     # Attempt direct cast for integer types after range check
                     if np.issubdtype(image_array.dtype, np.integer):
                         min_val, max_val = np.min(image_array), np.max(image_array)
                         if min_val < 0 or max_val > 255:
                              print(f"Warning: Integer data out of uint8 range [{min_val}, {max_val}]. Clipping.")
                              image_array_conv = np.clip(image_array, 0, 255).astype(np.uint8)
                         else:
                              image_array_conv = image_array.astype(np.uint8)
                     else: # For other types, attempt standard conversion after ensuring finite
                          image_array_finite = np.nan_to_num(image_array, nan=0.0, posinf=255.0, neginf=0.0)
                          image_array_conv = np.clip(image_array_finite,0,255).astype(np.uint8)
                 except (ValueError, OverflowError, TypeError) as e:
                      print(f"Warning: Could not convert final image dtype {image_array.dtype} to uint8. Error: {e}");
                      return False
             image_array = image_array_conv # Use the converted array
        
        # Validate image shape for RGB
        if image_array.ndim != 3 or image_array.shape[0] == 0 or image_array.shape[1] == 0 or image_array.shape[2] != 3:
             print(f"Error: Final image is not valid RGB (shape: {image_array.shape}). Cannot save.");
             return False

        # Double-check for NaN/Inf again after potential type conversions (should be caught by nan_to_num)
        if np.any(np.isnan(image_array)) or np.any(np.isinf(image_array)):
             print(f"Error: Final image array contains NaN or Inf values after conversion attempts. Cannot save.")
             return False

        pil_image = Image.fromarray(image_array, mode='RGB')
        pil_image.save(save_path)
        print(f"Successfully saved image to {save_path}")
        return True
    except Exception as e:
        print(f"Error saving image to {save_path}: {e}");
        traceback.print_exc(limit=2);
        return False


def map_photons(
    image_source: Union[str, np.ndarray],
    photon_chunk_files: List[str],
    save_path: str,
    dest_logical_bounds: Optional[Union[List[float], Tuple[float, ...]]] = None,      # Format: [y_min(horiz), y_max(horiz), z_min(vert), z_max(vert)]
    pixels_per_logical_unit: Optional[Union[float, int]] = None, # Resolution scale parameter
    output_shape: Optional[Tuple[int, int]] = None,           # Optional, (height, width), only used if dest_logical_bounds is None
    default_color: Union[int, float, np.number, Tuple[Any, ...]] = (0, 0, 0),
    epsilon: float = 1e-9,
    flip_z_axis_render: bool = True # Parameter to control Z-axis rendering direction
) -> bool:
    """Map photon data from chunks to create a lensed image from a source image.

    Processes photon trajectory data (from .npy chunk files) which links points on
    an observer's window (destination) to points on a background source image.
    It then samples colors from the source image and accumulates them onto a
    final rendered image, applying specified coordinate mappings and resolution.
    Supports two main modes for defining the output image dimensions and mapping:
    1. Manual Bounds + Scale: Uses 'dest_logical_bounds' and 'pixels_per_logical_unit'.
    2. Auto Bounds + Manual Shape: Uses 'output_shape'; bounds are auto-calculated.

    The Z-axis rendering can be flipped to achieve a "bottom-up" view.

    :param image_source: Path to the source image file or a NumPy array (H,W[,C]) of image data.
    :param photon_chunk_files: List of file paths to .npy chunk files containing photon data.
                               Each row in a chunk is expected to be (y0, z0, y1, z1), where
                               (y0,z0) are destination coords and (y1,z1) are source coords.
    :param save_path: Path where the output lensed image will be saved.
    :param dest_logical_bounds: Optional. Defines the destination mapping window in logical
                                coordinates [y_min(horiz), y_max(horiz), z_min(vert), z_max(vert)].
                                If set, 'pixels_per_logical_unit' MUST also be set.
    :param pixels_per_logical_unit: Optional. Resolution scale (pixels per logical unit).
                                      Required if 'dest_logical_bounds' is set.
    :param output_shape: Optional. Output image (height, width) in pixels.
                         Required if 'dest_logical_bounds' is None.
    :param default_color: Background RGB color for pixels not hit by any photons.
    :param epsilon: Small value for division-by-zero checks and float comparisons.
    :param flip_z_axis_render: If True, renders Z-axis so dest_map_bound_min_z maps to image bottom.
                               If False, min_z maps to image top.
    :return: True if the image was processed and saved successfully, False otherwise.
    :raises ValueError: If required parameters for a mode are missing or invalid.
    :raises TypeError: If 'photon_chunk_files' or 'image_source' (if path) is invalid type.
    :raises FileNotFoundError: If 'image_source' (if path) is not found.

    DocTests:
    >>> # Test with empty chunk list - should create a default background image
    >>> # Needs a valid (but can be dummy) save_path for the doctest framework.
    >>> import tempfile
    >>> with tempfile.NamedTemporaryFile(suffix='.png') as tmpfile:
    ...     result = map_photons(image_source=np.zeros((10,10,3), dtype=np.uint8), photon_chunk_files=[], save_path=tmpfile.name, output_shape=(5,5))
    >>> result
    True
    """
    start_time_total = time.time()
    mapped_image: Optional[np.ndarray] = None
    mapped_h, mapped_w = 1, 1 # Will be overwritten by logic below
    using_manual_dest_bounds = False # Flag to track operating mode

    # --- DEBUG COUNTERS ---\n
    total_photons_loaded_from_chunks = 0
    photons_lost_nan = 0
    photons_lost_source_bounds = 0
    photons_lost_dest_logical_filter = 0
    photons_lost_dest_pixel_map = 0
    photons_lost_invalid_source_idx = 0
    photons_accumulated = 0
    # --- END DEBUG COUNTERS ---\n

    # --- 0. Validate save_path & chunk list ---\n
    if not save_path or not isinstance(save_path, str): 
        raise ValueError("A valid string save_path is required.")
    if not isinstance(photon_chunk_files, list) or not all(isinstance(f, str) for f in photon_chunk_files):
        raise TypeError("photon_chunk_files must be a list of file path strings.")

    # --- Determine Mode and Validate Inputs / Set Output Shape ---\n
    dest_map_bound_min_y, dest_map_bound_max_y = 0.0, 0.0 # y now maps to horizontal axis
    dest_map_bound_min_z, dest_map_bound_max_z = 0.0, 0.0 # z now maps to vertical axis

    if dest_logical_bounds is not None:
        # --- Mode 1: Manual Bounds + Scale ---\n
        using_manual_dest_bounds = True
        print("Using MANUAL destination bounds + SCALE mode.")
        if flip_z_axis_render:
            print("  Z-axis rendering: FLIPPED (min_z to image bottom, max_z to image top).")
        else:
            print("  Z-axis rendering: NORMAL (min_z to image top, max_z to image bottom).")

        # Check required inputs for this mode
        if pixels_per_logical_unit is None:
            raise ValueError("Must provide 'pixels_per_logical_unit' when 'dest_logical_bounds' is set.")
        # Check for conflicting inputs
        if output_shape is not None:
            print("    Warning: 'output_shape' parameter provided but ignored because 'dest_logical_bounds' is set.")

        # Validate dest_logical_bounds
        if not isinstance(dest_logical_bounds, (list, tuple)) or len(dest_logical_bounds) != 4:
            raise ValueError("dest_logical_bounds must be a list or tuple of four numbers: [y_min(horiz), y_max(horiz), z_min(vert), z_max(vert)]")
        try:
            # y_min, y_max, z_min, z_max for logical bounds
            dlb_y_min, dlb_y_max, dlb_z_min, dlb_z_max = [float(x) for x in dest_logical_bounds]
            if dlb_y_min > dlb_y_max or dlb_z_min > dlb_z_max:
                raise ValueError("Destination logical bounds must have min <= max for both y ([a,b], horizontal) and z ([c,d], vertical).")
            dest_map_bound_min_y = dlb_y_min
            dest_map_bound_max_y = dlb_y_max
            dest_map_bound_min_z = dlb_z_min # This is the logical minimum Z for the range
            dest_map_bound_max_z = dlb_z_max # This is the logical maximum Z for the range
            print(f"  Logical Bounds: y(horiz)=[{dlb_y_min:.3f}, {dlb_y_max:.3f}], z(vert)=[{dlb_z_min:.3f}, {dlb_z_max:.3f}]")
        except (ValueError, TypeError) as e:
             raise ValueError(f"Invalid dest_logical_bounds format or values: {e}")

        # Validate pixels_per_logical_unit
        try:
            scale = float(pixels_per_logical_unit)
            if scale <= 0:
                raise ValueError("pixels_per_logical_unit must be a positive number.")
            print(f"  Scale: {scale:.3f} pixels per logical unit")
        except (ValueError, TypeError):
            raise ValueError("pixels_per_logical_unit must be a positive number.")

        # Calculate output dimensions based on logical range and scale
        logical_width = dest_map_bound_max_y - dest_map_bound_min_y
        logical_height = dest_map_bound_max_z - dest_map_bound_min_z # This is the extent of the Z range

        # Use ceiling to ensure the pixel grid covers the entire logical range
        # Ensure minimum dimension is 1 pixel
        mapped_w = max(1, int(math.ceil(logical_width * scale)))
        mapped_h = max(1, int(math.ceil(logical_height * scale)))

        print(f"--> Calculated Output Shape (Preserving Aspect Ratio): (Height={mapped_h}, Width={mapped_w})")

    else:
        # --- Mode 2: Auto Bounds + Manual Shape ---\n
        print("Using AUTO destination bounds + SHAPE mode.")
        if flip_z_axis_render:
            print("  Z-axis rendering: FLIPPED (min_z to image bottom, max_z to image top).")
        else:
            print("  Z-axis rendering: NORMAL (min_z to image top, max_z to image bottom).")

        # Check required input for this mode
        if output_shape is None:
             raise ValueError("Must provide 'output_shape' when 'dest_logical_bounds' is not set.")
        # Check for conflicting inputs
        if pixels_per_logical_unit is not None:
             print("    Warning: 'pixels_per_logical_unit' parameter provided but ignored because 'dest_logical_bounds' is not set.")

        # Validate output_shape
        if not isinstance(output_shape, tuple) or len(output_shape) != 2:
            raise ValueError("output_shape must be a tuple of (height, width).")
        try:
            mapped_h, mapped_w = int(output_shape[0]), int(output_shape[1])
            if mapped_h <= 0 or mapped_w <= 0:
                raise ValueError("output_shape dimensions must be positive integers.")
        except (ValueError, TypeError):
            raise ValueError("output_shape dimensions must be positive integers.")
        print(f"  Using provided Output Shape: (Height={mapped_h}, Width={mapped_w})")


    # --- Handle Empty Chunk List (Uses calculated/provided mapped_h/w) ---\n
    if not photon_chunk_files:
        print("Warning: photon_chunk_files list is empty. Nothing to process. Creating default background image.")
        try:
            default_color_rgb_tuple: Tuple[int, int, int]
            if isinstance(default_color, (int,float,np.number)): default_color_rgb_tuple = (int(default_color),)*3 # type: ignore
            elif isinstance(default_color, (tuple, list)) and len(default_color) == 3: default_color_rgb_tuple = tuple(int(c) for c in default_color) # type: ignore
            else: default_color_rgb_tuple = (0, 0, 0); print("Warning: Invalid default_color. Using (0,0,0).")
            default_color_rgb_tuple = tuple(np.clip(c, 0, 255) for c in default_color_rgb_tuple) # Ensure values are valid uint8
            mapped_image_arr = np.full((mapped_h, mapped_w, 3), default_color_rgb_tuple, dtype=np.uint8)
            return _save_image(mapped_image_arr, save_path)
        except Exception as e: print(f"Error setting up background image dimensions: {e}"); return False


    # --- Continue Processing ---\n
    try:
        # --- 1. Calculate Output Aspect Ratio ---\n
        output_pixel_aspect = mapped_w / mapped_h if mapped_h > 0 else 1.0

        # --- 2. Load Source Image ---\n
        original_image_loaded: np.ndarray
        if isinstance(image_source, str):
            try:
                with Image.open(image_source) as img:
                    # img_mode = img.mode; # Original mode, not used further here
                    if img.mode != 'RGB': img = img.convert('RGB')
                    original_image_loaded = np.array(img)
            except FileNotFoundError: raise FileNotFoundError(f"Source image file not found: {image_source}")
            except UnidentifiedImageError: raise ValueError(f"Could not identify or open image file: {image_source}") # Changed to ValueError
            except Exception as e: raise ValueError(f"Error loading image file {image_source}: {e}") # Changed to ValueError
        elif isinstance(image_source, np.ndarray):
            original_image_loaded = image_source.copy()
        else:
            raise TypeError("image_source must be a file path (string) or a NumPy array.")

        # --- Force RGB, Get Source Dimensions, Aspect ---\n
        original_image_rgb: np.ndarray
        if original_image_loaded.ndim == 2: original_image_rgb = np.stack([original_image_loaded] * 3, axis=-1)
        elif original_image_loaded.ndim == 3 and original_image_loaded.shape[2] == 1: original_image_rgb = np.repeat(original_image_loaded, 3, axis=-1)
        elif original_image_loaded.ndim == 3 and original_image_loaded.shape[2] == 3: original_image_rgb = original_image_loaded
        elif original_image_loaded.ndim == 3 and original_image_loaded.shape[2] == 4: original_image_rgb = original_image_loaded[:, :, :3] # Drop alpha channel
        else: raise ValueError(f"Unsupported source image format (shape: {original_image_loaded.shape})")

        orig_h, orig_w, _ = original_image_rgb.shape
        if orig_h == 0 or orig_w == 0: raise ValueError("Source image dimensions cannot be zero.")
        source_pixel_aspect = orig_w / orig_h

        # --- Rescale Logic for low-contrast or non-uint8 source images ---\n
        final_image_for_sampling = original_image_rgb
        try:
            source_max_val = np.max(final_image_for_sampling); source_min_val = np.min(final_image_for_sampling)
            needs_rescaling = False; rescale_threshold = 32 # Rescale if max value is low (e.g., < 32 for uint8)
            if final_image_for_sampling.dtype == np.uint8 and source_max_val > source_min_val and source_max_val < rescale_threshold: needs_rescaling = True

            if needs_rescaling:
                 print("    Rescaling low-contrast source image to full range [0, 255].")
                 if source_max_val > source_min_val:
                     scaled_image_float = ((final_image_for_sampling.astype(np.float64) - source_min_val) / (source_max_val - source_min_val) * 255.0)
                 else: # Handle flat image case (all pixels same color)
                     scaled_image_float = np.full_like(final_image_for_sampling, 128.0) # Assign mid-gray
                 final_image_for_sampling = np.clip(scaled_image_float, 0, 255).astype(np.uint8)
            elif final_image_for_sampling.dtype != np.uint8: # Ensure uint8 otherwise
                 print(f"    Converting source image from {final_image_for_sampling.dtype} to uint8, clipping to [0, 255].")
                 # Handle potential non-numeric data if not float/int before clipping
                 if np.issubdtype(final_image_for_sampling.dtype, np.floating) or np.issubdtype(final_image_for_sampling.dtype, np.integer):
                    final_image_for_sampling = np.clip(final_image_for_sampling, 0, 255).astype(np.uint8)
                 else: # Fallback for other dtypes, try direct conversion (might fail)
                    final_image_for_sampling = final_image_for_sampling.astype(np.uint8)
        except Exception as scale_e: print(f"  Warning: Could not perform source image range check/rescaling: {scale_e}.")

        original_image_rgb = final_image_for_sampling # Use potentially rescaled image

        print(f"Source image loaded: shape={original_image_rgb.shape}, aspect={source_pixel_aspect:.3f}")
        print(f"Output shape set to: {(mapped_h, mapped_w)}, aspect={output_pixel_aspect:.3f}")


        # --- 4. Prepare Default Color (as tuple of ints) ---\n
        default_color_rgb_final: Tuple[int, int, int]
        if isinstance(default_color, (int, float, np.number)): default_color_rgb_final = (int(default_color),) * 3 # type: ignore
        elif isinstance(default_color, (tuple, list)) and len(default_color) == 3:
            try: default_color_rgb_final = tuple(int(c) for c in default_color) # type: ignore
            except (ValueError, TypeError): default_color_rgb_final = (0,0,0); print("Warning: Could not parse default_color components. Using (0,0,0).")
        else: default_color_rgb_final = (0, 0, 0); print("Warning: Invalid default_color format. Using (0,0,0).")
        try:
            default_color_rgb_final = tuple(np.clip(c, 0, 255) for c in default_color_rgb_final)
        except (ValueError, TypeError):
            default_color_rgb_final = (0, 0, 0); print("Warning: Could not parse/clip default_color. Using (0,0,0).")


        # --- 5. Initialize Accumulation Buffers (Uses determined mapped_h/w) ---\n
        print("Initializing accumulation buffers...")
        sum_array = np.zeros((mapped_h, mapped_w, 3), dtype=np.float64)
        count_array = np.zeros((mapped_h, mapped_w), dtype=np.int64)
        # Initialize mapped_image with the final default color
        mapped_image_output = np.full((mapped_h, mapped_w, 3), default_color_rgb_final, dtype=np.uint8)

        # --- Determine Source & Auto Destination Bounds (Scan needed for source + auto dest) ---\n
        print("Scanning chunks to determine coordinate bounds...")
        global_max_abs_y0=epsilon; global_max_abs_z0=epsilon # Needed only for Auto destination mode
        global_max_abs_y1=epsilon; global_max_abs_z1=epsilon # Needed for Source bounds (y1,z1 are source coords)
        scan_start_time = time.time()
        photons_found_in_scan = False
        for chunk_idx, chunk_file in enumerate(photon_chunk_files):
             try:
                photon_list_chunk = np.load(chunk_file)
                if photon_list_chunk.ndim != 2 or photon_list_chunk.shape[1] != 4 or photon_list_chunk.shape[0] == 0: continue
                # Photon chunk data: y0(dest_horiz), z0(dest_vert), y1(src_horiz), z1(src_vert)
                y0c, z0c, y1c, z1c = photon_list_chunk[:,0], photon_list_chunk[:,1], photon_list_chunk[:,2], photon_list_chunk[:,3]
                nan_mask_chunk = np.isnan(y0c)|np.isnan(z0c)|np.isnan(y1c)|np.isnan(z1c)
                valid_indices = ~nan_mask_chunk
                if not np.any(valid_indices): continue
                photons_found_in_scan = True
                if not using_manual_dest_bounds: # Auto destination bounds mode
                    global_max_abs_y0=max(global_max_abs_y0, np.max(np.abs(y0c[valid_indices])))
                    global_max_abs_z0=max(global_max_abs_z0, np.max(np.abs(z0c[valid_indices])))
                # Source bounds always needed
                global_max_abs_y1=max(global_max_abs_y1, np.max(np.abs(y1c[valid_indices])))
                global_max_abs_z1=max(global_max_abs_z1, np.max(np.abs(z1c[valid_indices])))
                del photon_list_chunk # Free memory
             except FileNotFoundError: print(f"\nWarning: File not found during scan: {chunk_file}. Skipping."); continue
             except Exception as e: print(f"\nWarning: Error scanning {chunk_file} for bounds: {e}")
        print(f"\nBounds scan complete. Time: {time.time() - scan_start_time:.2f}s")

        if not photons_found_in_scan:
             print("\nError: No valid (non-NaN) photons found in any chunk during scan. Cannot proceed.")
             return _save_image(mapped_image_output, save_path) # Save default background


        # --- Calculate SOURCE Mapping Bounds (y1, z1 mapping to source image pixels) ---\n
        # These bounds define the extent of the source image sampled by photons.
        source_bound_y1_extent = max(global_max_abs_y1, epsilon)
        source_bound_z1_extent = max(global_max_abs_z1, epsilon)
        # Adjust bounds to match source image aspect ratio to avoid distortion when mapping
        if source_pixel_aspect >= 1.0: # Source image is wider or square
             temp_source_map_bound_y1 = source_bound_y1_extent
             temp_source_map_bound_z1 = max(source_bound_z1_extent, source_bound_y1_extent / source_pixel_aspect)
        else: # Source image is taller
             temp_source_map_bound_z1 = source_bound_z1_extent
             temp_source_map_bound_y1 = max(source_bound_y1_extent, source_bound_z1_extent * source_pixel_aspect)
        source_map_bound_y1 = temp_source_map_bound_y1
        source_map_bound_z1 = temp_source_map_bound_z1
        # Denominators for normalizing source coordinates to [0,1] range (actually [-bound, +bound] to [0, range])
        source_denom_y1 = max(2.0 * source_map_bound_y1, epsilon) 
        source_denom_z1 = max(2.0 * source_map_bound_z1, epsilon)
        print(f"Source mapping bounds (y1,z1) (Aspect Corrected): y1(horiz):[+/-{source_map_bound_y1:.3f}], z1(vert):[+/-{source_map_bound_z1:.3f}]")


        # --- Calculate DESTINATION Mapping Denominators (y0, z0 mapping to output image pixels) ---\n
        dest_denom_y0, dest_denom_z0 = epsilon, epsilon

        if using_manual_dest_bounds:
            # Mode 1: dest_map_bound_min_z/max_z, etc. already set from input 'dest_logical_bounds'
            dest_range_y0 = dest_map_bound_max_y - dest_map_bound_min_y
            dest_range_z0 = dest_map_bound_max_z - dest_map_bound_min_z # Logical range for Z
            dest_denom_y0 = max(dest_range_y0, epsilon)
            dest_denom_z0 = max(dest_range_z0, epsilon) # Denominator for Z mapping
        else: # Mode 2 (Auto Bounds from photon data, adjusted for output aspect ratio)
            logical_dest_bound_y0_extent = max(global_max_abs_y0, epsilon)
            logical_dest_bound_z0_extent = max(global_max_abs_z0, epsilon)
            print(f"Auto-calculating destination bounds from max extent: y0(horiz)=+/-{logical_dest_bound_y0_extent:.3f}, z0(vert)=+/-{logical_dest_bound_z0_extent:.3f}")
            temp_dest_map_bound_y0, temp_dest_map_bound_z0 = 0.0, 0.0
            if output_pixel_aspect >= 1.0: # Output image is wider or square
                 temp_dest_map_bound_y0 = logical_dest_bound_y0_extent
                 temp_dest_map_bound_z0 = max(logical_dest_bound_z0_extent, logical_dest_bound_y0_extent / output_pixel_aspect)
            else: # Output image is taller
                 temp_dest_map_bound_z0 = logical_dest_bound_z0_extent
                 temp_dest_map_bound_y0 = max(logical_dest_bound_y0_extent, logical_dest_bound_z0_extent * output_pixel_aspect)
            # For auto mode, bounds are symmetric around 0
            dest_map_bound_min_y = -temp_dest_map_bound_y0
            dest_map_bound_max_y = temp_dest_map_bound_y0
            dest_map_bound_min_z = -temp_dest_map_bound_z0 # This is logical_min_z for auto mode
            dest_map_bound_max_z = temp_dest_map_bound_z0 # This is logical_max_z for auto mode
            dest_denom_y0 = max(2.0 * temp_dest_map_bound_y0, epsilon)
            dest_denom_z0 = max(2.0 * temp_dest_map_bound_z0, epsilon) # Denominator for Z mapping
            print(f"Destination mapping bounds (y0,z0) (Auto, Aspect Corrected): y(horiz)=[{dest_map_bound_min_y:.3f},{dest_map_bound_max_y:.3f}], z(vert)=[{dest_map_bound_min_z:.3f},{dest_map_bound_max_z:.3f}]")

        print(f"Final Mapping Denominators: src_y1(horiz)={source_denom_y1:.3f}, src_z1(vert)={source_denom_z1:.3f}, dst_y0(horiz)={dest_denom_y0:.3f}, dst_z0(vert)={dest_denom_z0:.3f}")


        # --- Loop Through Chunks for Processing ---\n
        print("\nProcessing photon data chunk by chunk...")
        start_proc_time = time.time()
        for chunk_idx, chunk_file in enumerate(photon_chunk_files):
            if (chunk_idx + 1) % 10 == 0 or chunk_idx == 0 or chunk_idx == len(photon_chunk_files) - 1:
                 print(f"\rProcessing chunk {chunk_idx + 1}/{len(photon_chunk_files)}: {os.path.basename(chunk_file)}...", end="")
            try:
                photon_list = np.load(chunk_file)
                count_loaded_chunk = 0
                if photon_list.ndim == 2 and photon_list.shape[1] == 4 and photon_list.shape[0] > 0:
                    count_loaded_chunk = photon_list.shape[0]
                    total_photons_loaded_from_chunks += count_loaded_chunk
                else:
                    print(f"\nWarning: Invalid shape in chunk {chunk_file}: {photon_list.shape}. Skipping.")
                    del photon_list; continue
                
                # Extract coordinates: y0 (dest_horiz), z0 (dest_vert), y1 (src_horiz), z1 (src_vert)
                y0_f, z0_f, y1_f, z1_f = photon_list[:, 0], photon_list[:, 1], photon_list[:, 2], photon_list[:, 3]
                
                # Filter out NaN values from photon data
                nan_mask = np.isnan(y0_f) | np.isnan(z0_f) | np.isnan(y1_f) | np.isnan(z1_f)
                valid_indices_nan = np.where(~nan_mask)[0]
                count_after_nan = len(valid_indices_nan)
                photons_lost_nan += (count_loaded_chunk - count_after_nan)
                if count_after_nan == 0: del photon_list; continue
                
                # Filter photons based on source mapping bounds (y1, z1)
                source_bounds_ok = (np.abs(y1_f[valid_indices_nan]) <= source_map_bound_y1 + epsilon) & \
                                   (np.abs(z1_f[valid_indices_nan]) <= source_map_bound_z1 + epsilon)
                valid_indices_src_bounds = valid_indices_nan[source_bounds_ok]
                count_after_src_bounds = len(valid_indices_src_bounds)
                photons_lost_source_bounds += (count_after_nan - count_after_src_bounds)
                if count_after_src_bounds == 0: del photon_list; continue
                
                # Filter photons based on destination logical bounds (y0, z0)
                y0_f_filt = y0_f[valid_indices_src_bounds]; z0_f_filt = z0_f[valid_indices_src_bounds]
                dest_logical_ok = (y0_f_filt >= dest_map_bound_min_y - epsilon) & (y0_f_filt <= dest_map_bound_max_y + epsilon) & \
                                  (z0_f_filt >= dest_map_bound_min_z - epsilon) & (z0_f_filt <= dest_map_bound_max_z + epsilon)
                valid_indices_dest_logical = valid_indices_src_bounds[dest_logical_ok]
                count_after_dest_logical = len(valid_indices_dest_logical)
                photons_lost_dest_logical_filter += (count_after_src_bounds - count_after_dest_logical)
                if count_after_dest_logical == 0: del photon_list; continue
                
                # Valid photons after all filters for this chunk
                y0_f_pv = y0_f[valid_indices_dest_logical]; z0_f_pv = z0_f[valid_indices_dest_logical]
                y1_f_pv = y1_f[valid_indices_dest_logical]; z1_f_pv = z1_f[valid_indices_dest_logical]

                # Map source coordinates (y1, z1) to source image pixel indices (row, col)
                # Normalizing y1, z1 from [-bound, +bound] to [0, image_dim-1]
                # Source image: y1 maps to columns (width), z1 maps to rows (height)
                # Wait, paper implies x = x, y=p cos phi', z = p sin phi'. So y is horizontal, z is vertical.
                # Original code has y1_physical (rows from z1) and z1_physical (cols from y1)
                # This seems consistent: z1 (vertical logical) -> row_indices_f (vertical pixel)
                #                       y1 (horizontal logical) -> col_indices_f (horizontal pixel)
                row_indices_f_src = ((z1_f_pv + source_map_bound_z1) / source_denom_z1) * (orig_h - 1)
                col_indices_f_src = ((y1_f_pv + source_map_bound_y1) / source_denom_y1) * (orig_w - 1)
                # Convert to int and clip to source image dimensions
                source_pixel_rows = np.round(row_indices_f_src).astype(int)
                source_pixel_cols = np.round(col_indices_f_src).astype(int)
                source_pixel_rows = np.clip(source_pixel_rows, 0, orig_h - 1)
                source_pixel_cols = np.clip(source_pixel_cols, 0, orig_w - 1)

                # Map destination coordinates (y0, z0) to output image pixel indices
                # y0 maps to columns (width), z0 maps to rows (height) of the output image.
                y0_idx_f: np.ndarray; z0_idx_f: np.ndarray # Destination pixel coordinates (float)
                if using_manual_dest_bounds: # Mode 1: Manual logical bounds
                    # Map y0 from [min_y, max_y] to [0, mapped_w-1]
                    y0_idx_f = ((y0_f_pv - dest_map_bound_min_y) / dest_denom_y0) * (mapped_w - 1)
                    if flip_z_axis_render:
                        # Maps min_z to bottom row (mapped_h-1), max_z to top row (0)
                        z0_idx_f = ((dest_map_bound_max_z - z0_f_pv) / dest_denom_z0) * (mapped_h - 1)
                    else: # Original mapping: min_z to top row (0), max_z to bottom_row (mapped_h-1)
                        z0_idx_f = ((z0_f_pv - dest_map_bound_min_z) / dest_denom_z0) * (mapped_h - 1)
                else: # Mode 2 (Auto Bounds, symmetric around 0: [-bound, +bound])
                    # dest_map_bound_min_y = -auto_bound_y, dest_map_bound_max_y = +auto_bound_y
                    # dest_denom_y0 = 2 * auto_bound_y
                    auto_bound_y0 = dest_map_bound_max_y # Positive extent for y0
                    auto_bound_z0_positive_extent = dest_map_bound_max_z # Positive extent for z0

                    y0_idx_f = ((y0_f_pv + auto_bound_y0) / dest_denom_y0) * (mapped_w - 1)
                    if flip_z_axis_render:
                        # Maps -auto_bound_z0 (logical bottom) to bottom row (mapped_h-1)
                        #      +auto_bound_z0 (logical top) to top row (0)
                        z0_idx_f = ((auto_bound_z0_positive_extent - z0_f_pv) / dest_denom_z0) * (mapped_h - 1)
                    else: # Original mapping: -auto_bound_z0 to top row, +auto_bound_z0 to bottom row
                        z0_idx_f = ((z0_f_pv + auto_bound_z0_positive_extent) / dest_denom_z0) * (mapped_h - 1)

                # Convert to int and clip to output image dimensions
                dest_pixel_cols_clipped = np.round(y0_idx_f).astype(int)
                dest_pixel_rows_clipped = np.round(z0_idx_f).astype(int)
                dest_pixel_cols_clipped = np.clip(dest_pixel_cols_clipped, 0, mapped_w - 1)
                dest_pixel_rows_clipped = np.clip(dest_pixel_rows_clipped, 0, mapped_h - 1)
                
                # Filter out photons that map outside the pixel grid (should be minimal due to clipping and earlier logical filter)
                pixel_map_ok_mask = (y0_idx_f >= -epsilon) & (y0_idx_f < mapped_w + epsilon) & \
                                    (z0_idx_f >= -epsilon) & (z0_idx_f < mapped_h + epsilon)
                count_after_pixel_map = np.sum(pixel_map_ok_mask)
                photons_lost_dest_pixel_map += (count_after_dest_logical - count_after_pixel_map)
                if count_after_pixel_map == 0: del photon_list; continue
                
                # Final arrays of pixel indices for accumulation
                dest_pixel_cols_final = dest_pixel_cols_clipped[pixel_map_ok_mask]
                dest_pixel_rows_final = dest_pixel_rows_clipped[pixel_map_ok_mask]
                source_pixel_rows_final = source_pixel_rows[pixel_map_ok_mask]
                source_pixel_cols_final = source_pixel_cols[pixel_map_ok_mask]
                
                # This check is redundant if clipping to source_pixel_rows/cols is correct
                # invalid_src_idx_mask = (source_pixel_rows_final >= orig_h) | (source_pixel_cols_final >= orig_w) | \
                #                        (source_pixel_rows_final < 0) | (source_pixel_cols_final < 0)
                # num_invalid_src_idx = np.sum(invalid_src_idx_mask)
                # photons_lost_invalid_source_idx += num_invalid_src_idx
                # if num_invalid_src_idx > 0 :
                #      valid_accumulation_mask = ~invalid_src_idx_mask
                #      if not np.any(valid_accumulation_mask):
                #          del photon_list; continue
                #      dest_pixel_cols_final = dest_pixel_cols_final[valid_accumulation_mask]
                #      dest_pixel_rows_final = dest_pixel_rows_final[valid_accumulation_mask]
                #      source_pixel_rows_final = source_pixel_rows_final[valid_accumulation_mask]
                #      source_pixel_cols_final = source_pixel_cols_final[valid_accumulation_mask]
                
                count_accumulated_chunk = len(dest_pixel_cols_final)
                photons_accumulated += count_accumulated_chunk
                
                if count_accumulated_chunk > 0:
                    # Get colors from source image at (source_pixel_rows_final, source_pixel_cols_final)
                    source_colors = original_image_rgb[source_pixel_rows_final, source_pixel_cols_final].astype(np.float64)
                    # Accumulate colors and counts at (dest_pixel_rows_final, dest_pixel_cols_final)
                    np.add.at(sum_array, (dest_pixel_rows_final, dest_pixel_cols_final), source_colors)
                    np.add.at(count_array, (dest_pixel_rows_final, dest_pixel_cols_final), 1)
                
                # Clean up large arrays from this chunk
                del photon_list, y0_f, z0_f, y1_f, z1_f, valid_indices_nan, valid_indices_src_bounds, valid_indices_dest_logical
                del y0_f_pv, z0_f_pv, y1_f_pv, z1_f_pv, source_pixel_rows, source_pixel_cols
                del dest_pixel_cols_clipped, dest_pixel_rows_clipped, pixel_map_ok_mask
                del dest_pixel_cols_final, dest_pixel_rows_final, source_pixel_rows_final, source_pixel_cols_final
                if count_accumulated_chunk > 0: del source_colors

            except FileNotFoundError: print(f"\nWarning: Chunk file not found during processing: {chunk_file}. Skipping."); continue
            except Exception as e: print(f"\nWarning: Error processing chunk {chunk_file}: {e}. Skipping."); traceback.print_exc(limit=1, file=sys.stdout); continue
        
        print(f"\nFinished processing all chunks. Time: {time.time() - start_proc_time:.2f}s")
        print("\n--- Photon Loss Debug Summary ---")
        print(f"Destination Mode: {'MANUAL Bounds + Scale' if using_manual_dest_bounds else 'AUTO Bounds + Shape'}")
        print(f"  Z-axis Render Flip: {flip_z_axis_render}")
        print(f"Total photons loaded from chunks: {total_photons_loaded_from_chunks}")
        print(f"  Lost due to NaN values:         {photons_lost_nan}")
        count_rem_1 = total_photons_loaded_from_chunks - photons_lost_nan
        print(f"  Remaining after NaN filter:     {count_rem_1}")
        # Corrected f-string for Python 3.7 compatibility
        print(
            f"  Lost due to source bounds filter [y1(horiz):+/-{source_map_bound_y1:.3f}, "
            f"z1(vert):+/-{source_map_bound_z1:.3f}]: {photons_lost_source_bounds}"
        )
        count_rem_2 = count_rem_1 - photons_lost_source_bounds
        print(f"  Remaining after source bounds:  {count_rem_2}")
        # Corrected f-string for Python 3.7 compatibility
        print(
            f"  Lost due to DEST LOGICAL filter [y0(horiz)=[{dest_map_bound_min_y:.3f},{dest_map_bound_max_y:.3f}], "
            f"z0(vert)=[{dest_map_bound_min_z:.3f},{dest_map_bound_max_z:.3f}]]: {photons_lost_dest_logical_filter}"
        )
        count_rem_3 = count_rem_2 - photons_lost_dest_logical_filter
        print(f"  Remaining after dest logical:   {count_rem_3}")
        # Corrected f-string for Python 3.7 compatibility
        print(
            f"  Lost due to dest pixel map/clip (outside vert=[{0}-{mapped_h-1}], "
            f"horiz=[{0}-{mapped_w-1}]): {photons_lost_dest_pixel_map}"
        )
        count_rem_4 = count_rem_3 - photons_lost_dest_pixel_map
        print(f"  Remaining after dest pixel map: {count_rem_4}")
        # photons_lost_invalid_source_idx is not actively incremented if the redundant check is removed
        # print(f"  Lost due to invalid source idx (post-map/clip check): {photons_lost_invalid_source_idx}")
        print(f"Total photons accumulated for rendering: {photons_accumulated}")
        print(f"Calculated total loss (before accumulation): {total_photons_loaded_from_chunks - photons_accumulated}")
        calc_total_processed = photons_lost_nan + photons_lost_source_bounds + photons_lost_dest_logical_filter + photons_lost_dest_pixel_map + photons_lost_invalid_source_idx + photons_accumulated
        if abs(calc_total_processed - total_photons_loaded_from_chunks) > 0:
             print(f"!!! WARNING: Photon count sanity check failed! Calculated Total Processed ({calc_total_processed}) != Loaded ({total_photons_loaded_from_chunks}). Diff: {total_photons_loaded_from_chunks - calc_total_processed}")
        print("----------------------------------")
        print("Calculating averages and finalizing image...")
        
        # Finalize image by averaging accumulated colors
        hit_mask = count_array > 0
        num_hit_pixels = np.sum(hit_mask)
        if num_hit_pixels > 0:
            print(f"  Averaging colors for {num_hit_pixels} hit pixels.")
            valid_counts = count_array[hit_mask]
            average_colors_float = sum_array[hit_mask] / valid_counts[..., np.newaxis]
            mapped_image_output[hit_mask] = np.clip(average_colors_float, 0, 255).astype(np.uint8)
        else:
            print("  No pixels were hit by valid photons. Output will be default background color.")
        
        end_time_total = time.time()
        print(f"Processing finished. Total time: {end_time_total - start_time_total:.2f} seconds.")
        return _save_image(mapped_image_output, save_path)
    
    except Exception as e:
        print(f"\n--- A critical error occurred during mapping setup or processing: {e} ---"); traceback.print_exc()
        # Try to save a default/partial image if a critical error occurs mid-process
        if 'mapped_h' in locals() and 'mapped_w' in locals() and save_path:
             try:
                 final_default_color_rgb_error: Tuple[int, int, int] = (0,0,0)
                 if 'default_color_rgb_final' in locals(): final_default_color_rgb_error = default_color_rgb_final
                 
                 error_image_to_save: Optional[np.ndarray] = None
                 if 'mapped_image_output' in locals() and mapped_image_output is not None and mapped_image_output.shape == (mapped_h, mapped_w, 3):
                     error_image_to_save = mapped_image_output
                 else: # Fallback to creating a new default image
                     print("Initializing empty background for error image due to inconsistent state.")
                     error_image_to_save = np.full((mapped_h, mapped_w, 3), final_default_color_rgb_error, dtype=np.uint8)
                 
                 print("Attempting to save background/partial image due to critical error.")
                 _save_image(error_image_to_save, save_path + "_CRITICAL_ERROR.png")
             except Exception as save_err:
                 print(f"Failed to save error image: {save_err}")
        return False

if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
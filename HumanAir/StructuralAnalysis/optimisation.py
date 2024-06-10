import numpy as np
from scipy.optimize import minimize  # type: ignore[import-untyped]
from scipy.integrate import cumulative_trapezoid  # type: ignore[import-untyped]
import sys
import os
import matplotlib.pyplot as plt  # noqa: F401
import time

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.StructuralAnalysis.WingStructure import WingStructure
from HumanAir.StructuralAnalysis.LoadDistributions import (
    InternalLoads,
    force_distribution,
    moment_distribution,
    weight_distribution,
    interpolate_Cl_Cd_Cm,
    get_deflection,
)
from HumanAir.aircraft_data import aircraft_data, airfoil_shape, Cl_data_wing, Cdi_data_wing, Cm_data_wing

# Caching for SPEEEEEEEEEEEEEEEEEEEED
# force_distributions = None
# force_dist_params = None # Check whether parameters to get the distribution have changed


def get_stringers_at_nodes(stringer_sections, no_stringers, nodes_halfspan, centre_node=True):
    if len(stringer_sections) != len(no_stringers):
        raise ValueError("The number of stringer sections should be equal to the number of stringers per section")

    stringers = np.ones(nodes_halfspan)
    stringer_section_ends = np.rint(np.cumsum(stringer_sections) * len(stringers)).astype(int)
    stringer_section_starts = np.insert(stringer_section_ends[:-1], 0, 0)
    for i in range(len(stringer_sections)):
        start = stringer_section_starts[i]
        end = stringer_section_ends[i]
        # print(start, end)
        stringers[start:end] = no_stringers[i]

    if centre_node:
        stringers = np.append(stringers, no_stringers[-1])

    return stringers


def objective(variables, htot, wtot, rho, len_nodes, A_stringer, stringer_sections_halfspan, n_halfspan, centre_node):
    t_spar_tip, t_spar_root, t_skin, no_stringers = unstack_variables(variables)
    stringers_at_nodes = get_stringers_at_nodes(
        stringer_sections_halfspan, no_stringers, n_halfspan, centre_node=centre_node
    )

    t_spar = np.linspace(t_spar_root, t_spar_tip, n_halfspan)

    W_spar = np.sum(rho * t_spar * htot * len_nodes)
    W_skin = np.sum(rho * t_skin * wtot * len_nodes)
    W_stringers = np.sum(rho * A_stringer * stringers_at_nodes * len_nodes)
    weight = W_skin + W_spar + W_stringers
    return weight





def get_axial_forces(Mx, Vy, MOI, h_max, stringers_at_nodes, A_stringer):
    sigma_bending = Mx * h_max / MOI
    sigma_axial = Vy / (A_stringer * stringers_at_nodes)
    sigma = sigma_bending + sigma_axial

    print(f"{np.max(sigma)/1e6:.2f} MPa, with P = {np.max(sigma) * A_stringer} N", flush=True)
    P = sigma * A_stringer
    return P


def get_axial_stresses(Mx, MOI, h_max):
    sigma = Mx * h_max / MOI
    return sigma


def get_I(t_spar, t_skin, no_stringers, A_stringer, w_top, w_bottom, h_avemax, h_15c, h_50c):
    I_skin = t_skin * (w_top + w_bottom) * h_avemax**2
    I_spar = 1 / 12 * t_spar * (h_15c**3 + h_50c**3)
    I_stringers = A_stringer * no_stringers * h_avemax**2
    MOI = I_skin + I_spar + I_stringers
    return MOI


lengths_vars = None


def stack_variables(t_spar_tip, t_spar_root, t_skin, no_stringers):
    return np.hstack(([t_spar_tip, t_spar_root, t_skin], no_stringers)).flatten()


def unstack_variables(variables):
    t_spar_tip, t_spar_root, t_skin = variables[:3].flatten()
    no_stringers = variables[3:].flatten()
    return t_spar_tip, t_spar_root, t_skin, no_stringers


def MOI_from_variables(
    variables, A_stringer, w_top, w_bottom, h_avemax, h_15c, h_50c, n_halfspan, stringer_sections_halfspan, centre_node
):
    return MOI_and_stringers_from_variables(
        variables,
        A_stringer,
        w_top,
        w_bottom,
        h_avemax,
        h_15c,
        h_50c,
        n_halfspan,
        stringer_sections_halfspan,
        centre_node,
    )[0]


def MOI_and_stringers_from_variables(
    variables, A_stringer, w_top, w_bottom, h_avemax, h_15c, h_50c, n_halfspan, stringer_sections_halfspan, centre_node
):
    t_spar_root, t_spar_tip, t_skin, no_stringers = unstack_variables(variables)
    stringers_at_nodes = get_stringers_at_nodes(
        stringer_sections_halfspan, no_stringers, n_halfspan, centre_node=centre_node
    )
    t_spar = np.linspace(t_spar_root, t_spar_tip, n_halfspan)
    MOI = get_I(t_spar, t_skin, stringers_at_nodes, A_stringer, w_top, w_bottom, h_avemax, h_15c, h_50c)

    return MOI, stringers_at_nodes


def deflection_constraint(variables, y, Mx, max_deflect, E, MOI_args):
    MOI = MOI_from_variables(variables, *MOI_args)
    deflection = get_deflection(MOI, y, Mx, E)
    return max_deflect - np.max(deflection)


def stringer_buckling_constraint(variables, Mx, Vy, A_stringer, h_15c, E, I_stringer, L_eff, MOI_args):
    P_cr = (np.pi) ** 2 * E * I_stringer / L_eff**2  # critical force at which stringer will buckle
    h_max = h_15c / 2
    MOI, stringers = MOI_and_stringers_from_variables(variables, *MOI_args)

    # compressive stress
    P = get_axial_forces(Mx, Vy, MOI, h_max, stringers, A_stringer)
    # print(P_cr - np.max(P), flush=True)
    return P_cr - np.max(P)


def tensile_failure_constraint(variables, Mx, h_15c, sigma_yield, MOI_args):
    h_max = h_15c / 2
    MOI = MOI_from_variables(variables, *MOI_args)
    sigma = get_axial_stresses(Mx, MOI, h_max)
    return sigma_yield - np.max(sigma)


def get_force_distributions(AoA, altitude, Vc, Cl_data, Cdi_data, Cm_data, chord_dist, ac_data=aircraft_data):
    L_cruise, D_cruise = force_distribution(AoA, altitude, Vc, chord_dist, Cl_DATA=Cl_data, Cdi_DATA=Cdi_data)

    M_cruise = moment_distribution(AoA, altitude, Vc, chord_dist, Cm_data, ac_data=ac_data)

    W, W_fuel, idxs, _ = weight_distribution(chord_dist, ac_data=ac_data)

    return L_cruise, D_cruise, M_cruise, W, W_fuel, idxs


last_time = time.time()
start = time.time()
times_called = 0


def print_time():
    global last_time, times_called
    now = time.time()
    print(f"{times_called} - Time taken:", now - last_time, flush=True)
    last_time = now
    times_called += 1


# TODO: Make the AoA taken from aircraft_data, need to implement it in the json
def run_optimiser(
    ac_data=aircraft_data, AoA=-5.0, altitude=None, nlim=None, max_deflect=1.7, nodes=401, maxiter=None, full_return=False
):
    """
    Function to run the optimization of the wing structure

    Parameters
    ----------
    ac_data : dict
        Dictionary containing the aircraft data
    AoA : float, default -5.0
        Angle of attack of the aircraft [deg]
    altitude : float or None, default None
        Altitude of the aircraft [m].
        If None, the altitude is taken from the aircraft data (cruise altitude).
    nlim : float or None, default None
        Load factor of the aircraft
        If None, the load factor is taken from the aircraft data (nult/1.5).
    max_deflect : float, default 1.7
        Maximum deflection of the wing [m]
    nodes : int, default 401
        Number of nodes to discretize the wing structure
    maxiter : int or None, default None
        Maximum number of iterations for the optimization.
        If None, the default value of the optimizer is used.
    full_return : bool, default False
        If True, the function returns additional information about the optimization process

    Returns
    -------
    optimized_t_spar_tip, optimized_t_spar_root, optimized_t_skin : np.array
        Optimized thickness of the spar (tip and root) and skin
    optimized_no_stringers : np.array
        Optimized number of stringers
    weight : float
        Minimum weight of the wing structure (i.e. for the optimised structure)
    """

    altitude = ac_data["Performance"]["Altitude_Cruise_m"] if altitude is None else altitude
    nlim = ac_data["Performance"]["n_ult"] / 1.5 if nlim is None else nlim

    n_halfspan = nodes // 2

    # Initialize Wing Structure Class
    wing_structure = WingStructure(ac_data, airfoil_shape, nodes)
    h_mid, h_s1s2 = wing_structure.h_s1s2()
    l_box_up, l_box_down = wing_structure.d_s1s2()

    h_mid, h_s1s2 = h_mid[(h_mid.shape[0] // 2) :], h_s1s2[(h_s1s2.shape[0] // 2) :]
    l_box_up, l_box_down = l_box_up[(l_box_up.shape[0] // 2) :], l_box_down[(l_box_down.shape[0] // 2) :]
    chord_dist = wing_structure.chord_distribution
    y_points = wing_structure.ypts
    chord_halfspan, y_points_halfspan = chord_dist[(chord_dist.shape[0] // 2) :], y_points[(y_points.shape[0] // 2) :]

    # Parameters
    # rho = 2710  # density of aluminium (kg/m^3)

    h_15c = h_s1s2[:, 0]  # height of the spar at 15% of the chord, as a function of y
    h_50c = h_s1s2[:, 1]  # height of the spar at 50% of the chord, as a function of y
    htot = (h_15c + h_50c).flatten()  # total height, just for calculation ease, as a function of y
    h_avemax = htot / 4  # averaged "max" height from the central line, as a function of y
    w_top = l_box_up.flatten()  # width of top "straight" skin, as a function of y
    w_bottom = l_box_down.flatten()  # width of bottom "straight" skin, as a function of y
    wtot = w_bottom + w_top  # total width, just for calculation ease

    # E = 68e9  # Young's Modulus for aluminium (Pa)
    # sigma_yield = 40e6  # Yield strength for aluminium (Pa)

    stringer_sections_halfspan = [0.4, 0.3, 0.3]  # Should add up to 1
    no_string = [50, 30, 20]

    # Initial guess for thickness of spar
    t_spar = wing_structure.t_spar_dist.flatten()
    t_spar_root0, t_spar_tip0 = t_spar[0], t_spar[-1]
    t_skin0 = wing_structure.t_skin
    # t_spar0 = t_spar0[t_spar0.shape[0] // 2 :]
    # Constant skin thickness [m]
    # t_skin0 = wing_structure.t_skin * np.ones(t_spar0.shape)
    # no_stringers0 = get_stringers_at_nodes(
    #     stringer_sections_halfspan, no_string, n_halfspan, centre_node=(nodes % 2 == 1)
    # )
    # Initial guess
    # Spar thickness uniformly decreases from root to tip, skin thickness is constant
    t0 = stack_variables(t_spar_tip0, t_spar_root0, t_skin0, no_string)

    # Define the bounds for each variable
    extra = 1 if nodes % 2 == 1 else 0
    bounds_t_spar = [(0.0005, 0.01)] * 2
    bounds_t_skin = [(0.0005, 0.003)] # [(0.0005, 0.01)]
    bounds_no_stringers = [(1, 60)] * (len(no_string))
    bounds = np.array(bounds_t_spar + bounds_t_skin + bounds_no_stringers)

    Cl_DATA, Cdi_DATA, Cm_DATA = interpolate_Cl_Cd_Cm(Cl_data_wing, Cdi_data_wing, Cm_data_wing, y_points)

    material = ac_data["Geometry"]["wingbox_material"]

    L_cruise, D_cruise, M_cruise, W, W_fuel, idxs = get_force_distributions(
        AoA, altitude, ac_data["Performance"]["Vc_m/s"], Cl_DATA, Cdi_DATA, Cm_DATA, chord_dist, ac_data=ac_data
    )

    Vx, Vy, Vz, Mx, My, Mz = InternalLoads(
        L_cruise, D_cruise, M_cruise, wing_structure, ac_data=ac_data, load_factor=nlim
    )

    centre_node = nodes % 2 == 1

    rib_dists = np.diff(ac_data["Geometry"]["wing_rib_pos"], prepend=0.0, append=1.0)
    max_rib_dist_m = np.max(rib_dists) * wing_structure.b / 2

    # Constraints dictionary
    MOI_args = [
        wing_structure.stringer_area,
        w_top,
        w_bottom,
        h_avemax,
        h_15c,
        h_50c,
        n_halfspan,
        stringer_sections_halfspan,
        centre_node,
    ]

    constraints = [
        {
            "type": "ineq",
            "fun": deflection_constraint,
            "args": (y_points_halfspan, Mx, max_deflect, ac_data["Materials"][material]["E"], MOI_args),
        },
        # Rene says no stringer (column) buckling
        # {
        #     "type": "ineq",
        #     "fun": stringer_buckling_constraint,
        #     "args": (
        #         Mx,
        #         Vy,
        #         wing_structure.stringer_area,
        #         h_15c,
        #         ac_data["Materials"][material]["E"],
        #         ac_data["Geometry"]["wing_stringer_MOI_m4"],
        #         max_rib_dist_m,
        #         MOI_args,
        #     ),
        # },
        {
            "type": "ineq",
            "fun": tensile_failure_constraint,
            "args": (Mx, h_15c, ac_data["Materials"][material]["sigma_y"], MOI_args),
        },
    ]

    # Perform the optimization
    len_node_segment = wing_structure.b / (2 * n_halfspan)

    solver_options = {"maxiter": maxiter} if maxiter is not None else {}

    solution = minimize(
        objective,
        t0,
        args=(
            htot,
            wtot,
            ac_data["Materials"][material]["rho"],
            len_node_segment,
            wing_structure.stringer_area,
            stringer_sections_halfspan,
            n_halfspan,
            centre_node,
        ),
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options=solver_options,
    )

    # Optimized thickness values
    optimized_t_spar_tip, optimized_t_spar_root, optimized_t_skin, optimized_no_stringers = unstack_variables(solution.x)
    weight = solution.fun
    # I = get_I(optimized_t_spar, optimized_t_skin, optimized_no_stringers, wing_structure.stringer_area)
    # deflection = get_deflection(I, y_points_halfspan, Mx, aircraft_data["Materials"][material]["E"])

    if full_return:
        return (
            optimized_t_spar_tip,
            optimized_t_spar_root,
            optimized_t_skin,
            optimized_no_stringers,
            weight,
            y_points_halfspan,
            Mx,
            Vy,
            MOI_args,
            wing_structure.stringer_area,
            h_15c,
            max_deflect,
            max_rib_dist_m,
        )

    return optimized_t_spar_tip, optimized_t_spar_root, optimized_t_skin, optimized_no_stringers, weight


if __name__ == "__main__":

    (
        optimized_t_spar_tip,
        optimized_t_spar_root,
        optimized_t_skin,
        optimized_no_stringers,
        weight,
        y_points_halfspan,
        Mx,
        Vy,
        MOI_args,
        A_stringer,
        h_15c,
        max_deflect,
        max_rib_dist,
    ) = run_optimiser(aircraft_data, nodes=100, maxiter=None, max_deflect=1.7, full_return=True)
    print(max_rib_dist)

    # MOI, y, M, E
    material = aircraft_data["Geometry"]["wingbox_material"]

    MOI, stringer_dist = MOI_and_stringers_from_variables(
        stack_variables(optimized_t_spar_tip, optimized_t_spar_root, optimized_t_skin, optimized_no_stringers), *MOI_args
    )

    deflection = get_deflection(MOI, y_points_halfspan, Mx, aircraft_data["Materials"][material]["E"])
    axial_forces = get_axial_forces(Mx, Vy, MOI, h_15c / 2, stringer_dist, A_stringer)
    axial_stresses = get_axial_stresses(Mx, MOI, h_15c / 2)

    max_force_allowed = (
        (np.pi) ** 2
        * aircraft_data["Materials"][material]["E"]
        * aircraft_data["Geometry"]["wing_stringer_MOI_m4"]
        / max_rib_dist**2
    )
    sigma_yield = aircraft_data["Materials"][material]["sigma_y"]

    plt.figure()
    plt.plot(y_points_halfspan, deflection, "g-", label="Deflection")
    plt.plot(y_points_halfspan[[0, -1]], (max_deflect, max_deflect), "k--", label="$v_{max}$")

    # plt.figure()
    # plt.plot(y_points_halfspan, axial_forces, "r-", label="Axial Forces")
    # # plt.plot(y_points_halfspan, Mx, "y-", label="Bending Moment")
    # plt.plot(y_points_halfspan[[0, -1]], (max_force_allowed, max_force_allowed), "k--", label="$P_{cr}$")
    # plt.legend()

    plt.figure()
    plt.plot(y_points_halfspan, axial_stresses, "b-", label="Axial Stresses")
    plt.plot(y_points_halfspan[[0, -1]], (sigma_yield, sigma_yield), "k--", label="$\\sigma_{y}$")

    plt.show()

    # Print the results
    print("Optimized thickness for spar:", optimized_t_spar_root, " to ", optimized_t_spar_tip)
    print("Optimized thickness for skin:", optimized_t_skin)
    print("Optimized number of stringers:", optimized_no_stringers)
    print("Minimum weight:", weight)

    print("Total time taken:", time.time() - start, "s")

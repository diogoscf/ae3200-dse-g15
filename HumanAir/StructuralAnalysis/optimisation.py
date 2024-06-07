import numpy as np
from scipy.optimize import minimize
from scipy.integrate import cumulative_trapezoid
import sys
import os
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.StructuralAnalysis.WingStructure import WingStructure
from HumanAir.StructuralAnalysis.LoadDistributions import (
    InternalLoads,
    force_distribution,
    moment_distribution,
    weight_distribution,
    interpolate_Cl_Cd_Cm,
)
from HumanAir.aircraft_data import aircraft_data, airfoil_shape, Cl_data_wing, Cdi_data_wing, Cm_data_wing


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


def objective(variables, htot, wtot, rho, len_nodes, A_stringer):
    variables = variables.reshape(3, -1)
    t_spar = variables[0, :].flatten()
    t_skin = variables[1, :].flatten()
    no_stringers = variables[2, :].flatten()

    W_spar = np.sum(rho * t_spar * htot * len_nodes)
    W_skin = np.sum(rho * t_skin * wtot * len_nodes)
    W_stringers = np.sum(rho * A_stringer * no_stringers * len_nodes)
    weight = W_skin + W_spar + W_stringers
    return weight


def get_deflection(I, y, M, E):
    integrand = M / I
    # print(np.sum(I), np.sum(M), np.sum(integrand))
    theta = -1 / E * cumulative_trapezoid(integrand, y, initial=0)
    v = -cumulative_trapezoid(theta, y, initial=0)
    return v


def get_I(t_spar, t_skin, no_stringers, A_stringer, w_top, w_bottom, h_avemax, h_15c, h_50c):
    I_skin = t_skin * (w_top + w_bottom) * h_avemax**2
    I_spar = 1 / 12 * t_spar * (h_15c**3 + h_50c**3)
    I_stringers = A_stringer * no_stringers * h_avemax**2
    I = I_skin + I_spar + I_stringers
    return I


def deflection_constraint(
    variables, y, Mx, A_stringer, max_deflect, E, w_top, w_bottom, h_avemax, h_15c, h_50c, n_halfspan
):
    variables = variables.reshape(3, -1)
    t_spar = variables[0, :].flatten()
    t_skin = variables[1, :].flatten()
    no_stringers = variables[2, :].flatten()
    I = get_I(t_spar, t_skin, no_stringers, A_stringer, w_top, w_bottom, h_avemax, h_15c, h_50c)
    deflection = np.max(get_deflection(I, y, Mx, E))
    return max_deflect - deflection


if __name__ == "__main__":
    ########## Input ###########
    # Sw = 34.56  # [m^2]
    # taper_ratio = 0.4
    AoA = -6  # [deg]
    altitude = 3000  # [m]
    nlim = 3.8
    max_deflect = 1.7  # [m]
    # t1_spar = 0.010  # [m] thickness at the tip
    # t2_spar = 0.025  # [m] thickness at the root
    # t_skin = 0.007  # [m] thickness of skin
    # file_path = "HumanAir/WingBox/airfoil.txt"
    # file_path_y = "HumanAir/WingBox/Cl_DATA.txt"
    # A_str = 8e-5
    # Cr = 2.5  # [m] root chord length
    # b = 19.93
    # x_pos = np.array([0.15, 0.5])
    nodes = 401

    n_halfspan = nodes // 2

    # Initialize Torsional Stiffness Class
    wing_structure = WingStructure(aircraft_data, airfoil_shape, nodes)
    h_mid, h_s1s2 = wing_structure.h_s1s2()
    l_box_up, l_box_down = wing_structure.d_s1s2()

    h_mid, h_s1s2 = h_mid[h_mid.shape[0] // 2 :], h_s1s2[h_s1s2.shape[0] // 2 :]
    l_box_up, l_box_down = l_box_up[l_box_up.shape[0] // 2 :], l_box_down[l_box_down.shape[0] // 2 :]
    chord_dist = wing_structure.chord_distribution
    y_points = wing_structure.ypts
    chord_halfspan, y_points_halfspan = chord_dist[chord_dist.shape[0] // 2 :], y_points[y_points.shape[0] // 2 :]

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
    t_spar0 = wing_structure.t_spar_dist.flatten()
    t_spar0 = t_spar0[t_spar0.shape[0] // 2 :]
    # Constant skin thickness [m]
    t_skin0 = wing_structure.t_skin * np.ones(t_spar0.shape)
    no_stringers0 = get_stringers_at_nodes(
        stringer_sections_halfspan, no_string, n_halfspan, centre_node=(nodes % 2 == 1)
    )
    # Initial guess
    t0 = np.hstack((t_spar0, t_skin0, no_stringers0)).flatten()

    # Define the bounds for each variable
    extra = 1 if nodes % 2 == 1 else 0
    bounds_t_spar = [(0.005, 0.1)] * (n_halfspan + extra)
    bounds_t_skin = [(0.005, 0.1)] * (n_halfspan + extra)
    bounds_no_stringers = [(1, 100)] * (n_halfspan + extra)
    bounds = np.array(bounds_t_spar + bounds_t_skin + bounds_no_stringers)

    Cl_DATA, Cdi_DATA, Cm_DATA = interpolate_Cl_Cd_Cm(Cl_data_wing, Cdi_data_wing, Cm_data_wing, y_points)

    L_cruise, D_cruise = force_distribution(
        AoA, altitude, aircraft_data["Performance"]["Vc_m/s"], Cl_DATA=Cl_DATA, Cdi_DATA=Cdi_DATA
    )
    W_cruise, W_fuel, idxs, c_between_struts = weight_distribution(chord_dist, ac_data=aircraft_data)

    M_cruise = moment_distribution(
        aircraft_data["Performance"]["Vc_m/s"], altitude, Cm_DATA, AoA, ac_data=aircraft_data
    )

    Vx, Vy, Vz, Mx, My, Mz = InternalLoads(
        L_cruise, D_cruise, M_cruise, wing_structure, ac_data=aircraft_data, load_factor=1
    )

    material = aircraft_data["Geometry"]["wingbox_material"]

    # Constraints dictionary
    constraints = [
        {
            "type": "ineq",
            "fun": deflection_constraint,
            "args": (
                y_points_halfspan,
                Mx,
                wing_structure.stringer_area,
                max_deflect,
                aircraft_data["Materials"][material]["E"],
                w_top,
                w_bottom,
                h_avemax,
                h_15c,
                h_50c,
                n_halfspan,
            ),
        }
    ]

    # Perform the optimization
    len_node_segment = wing_structure.b / (2 * n_halfspan)

    solution = minimize(
        objective,
        t0,
        args=(
            htot,
            wtot,
            aircraft_data["Materials"][material]["rho"],
            len_node_segment,
            wing_structure.stringer_area,
        ),
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
    )

    # Optimized thickness values
    optimized_thickness = solution.x.reshape(3, -1)
    optimized_t_spar = optimized_thickness[0, :]
    optimized_t_skin = optimized_thickness[1, :]
    optimized_no_stringers = optimized_thickness[2, :]
    # I = get_I(optimized_t_spar, optimized_t_skin, optimized_no_stringers, wing_structure.stringer_area)
    # deflection = get_deflection(I, y_points_halfspan, Mx, aircraft_data["Materials"][material]["E"])

    # Print the results
    print("Optimized thickness for spar:", optimized_t_spar)
    print("Optimized thickness for skin:", optimized_t_skin)
    print("Optimized number of stringers:", optimized_no_stringers)
    print("Minimum weight:", solution.fun)

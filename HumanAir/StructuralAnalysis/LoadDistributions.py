import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.aircraft_data import aircraft_data, airfoil_shape, Cl_data_wing, Cm_data_wing, Cdi_data_wing
from HumanAir.StructuralAnalysis.WingStructure import WingStructure
from HumanAir.isa import isa
from HumanAir.unit_conversions import G


# Define the forces along half span
# def chord(Sw, taper_ratio, Cl_DATA, AoA, nodes):
#     b = Cl_DATA[AoA]["y_span"][-1] * 2
#     # Generate spanwise coordinate points
#     y = np.linspace(Cl_DATA[AoA]["y_span"][0], Cl_DATA[AoA]["y_span"][-1], nodes)  # n is the number of nodes
#     # Calculate the chord distribution
#     chord_length = 2 * Sw / (1 + taper_ratio) / b * (1 - (1 - taper_ratio) * np.abs(2 * y / b))
#     return chord_length, y


def force_distribution(AoA, altitude, V, Cl_DATA, Cdi_DATA, ac_data=aircraft_data):
    Sw = ac_data["Aero"]["S_Wing"]
    rho = isa(altitude)[2]
    Cl = np.array(Cl_DATA[AoA]["coefficient"])
    Cdi = np.array(Cdi_DATA[AoA]["coefficient"])
    L = Cl * 0.5 * rho * V**2 * Sw  # [N/m]
    D = Cdi * 0.5 * rho * V**2 * Sw  # [N/m]
    # plt.plot(Cl_DATA[AoA]["y_span"], L, label="Lift")
    # plt.ylim(0, np.max(L) * 1.1)
    # plt.show()
    return L, D


def weight_distribution(c, ac_data=aircraft_data):
    structureweight = ac_data["CL2Weight"]["Wing Weight"]
    fuelweight = ac_data["CL2Weight"]["Wfuel_N"]

    # Fuel is located between the struts
    strut_loc = ac_data["Geometry"]["strut_loc_b/2"]
    nodes_half_wing = c.shape[0] // 2
    n_before_strut = np.rint(nodes_half_wing * strut_loc).astype(int)
    extra = 0 if c.shape[0] % 2 == 0 else 1
    c_between_struts = c[nodes_half_wing - n_before_strut : nodes_half_wing + n_before_strut + extra]

    # Assume structure weight is distributed with chord
    # This is not a very good assumption, but it should be conservative
    W_avg_structure = structureweight / c.shape[0]
    W_structure = W_avg_structure * c / (np.mean(c))

    # Assume fuel weight is distributed with airfoil area, i.e. with chord^2
    W_avg_fuel = fuelweight / c_between_struts.shape[0]
    W_fuel = W_avg_fuel * c_between_struts**2 / (np.mean(c_between_struts**2))
    diff = c.shape[0] - c_between_struts.shape[0]
    W_fuel = np.pad(W_fuel, pad_width=diff // 2, mode="constant")
    # print(np.sum(W_fuel))

    return (
        W_structure + W_fuel,
        W_fuel,
        (nodes_half_wing - n_before_strut, nodes_half_wing + n_before_strut + extra),
        c_between_struts,
    )


def moment_distribution(V, altitude, Cm_DATA, AoA, ac_data=aircraft_data):
    Sw = ac_data["Aero"]["S_Wing"]
    MAC = ac_data["Aero"]["MAC_wing"]
    M = np.array(Cm_DATA[AoA]["coefficient"]) * 0.5 * isa(altitude)[2] * V**2 * Sw * MAC  # M_cruise
    return M


# axial forces distribution - during ground - so only the ends (0.6L)
Vy_Strut = 6670 / np.tan(0.546)  # deg in radian


def axial_distribution_ground(nodes, ac_data=aircraft_data):
    strut_loc = ac_data["Geometry"]["strut_loc_b/2"]
    nodes_half_wing = nodes // 2
    n_before_strut = np.rint(nodes_half_wing * strut_loc).astype(int)
    A_ground = np.zeros(nodes_half_wing)
    A_ground[:n_before_strut] = Vy_Strut
    return A_ground


# axial forces distribution - during flight - so only on the inner (inside 0.4L))
def axial_distribution_flight(nodes, ac_data=aircraft_data):
    strut_loc = ac_data["Geometry"]["strut_loc_b/2"]
    nodes_half_wing = nodes // 2
    n_before_strut = np.rint(nodes_half_wing * strut_loc).astype(int)
    A_cruise = np.zeros(nodes_half_wing)
    A_cruise[:n_before_strut] = -Vy_Strut
    return A_cruise


def strut_shear_result():
    return 6670


# def TestForces(Lcruise, W, Cl_DATA, AoA):
#     dy = (Cl_DATA[AoA]["y_span"][-1] - Cl_DATA[AoA]["y_span"][0]) / nodes
#     print("Ltot = ", np.trapz(Lcruise / 9.81, Cl_DATA[AoA]["y_span"]), " kg")
#     print("Wtot = ", np.trapz(W / 9.81, Cl_DATA[AoA]["y_span"]), " kg")


def InternalLoads(L, W, D, M, nodes, y_points, ac_data=aircraft_data, load_factor=1):
    nodes_half_wing = nodes // 2
    b_half = ac_data["Aero"]["b_Wing"] / 2
    sweep = ac_data["Aero"]["QuarterChordSweep_Wing_deg"]
    X_forces = -(D * (load_factor**2)) * b_half / (nodes_half_wing)  # drag and thrust act on the x axis (thrust = 0)
    Vx = integrate.cumulative_trapezoid(np.flip(X_forces[nodes_half_wing:]), y_points[nodes_half_wing:], initial=0)
    Vx = Vx[::-1]
    Z_forces = -(L * (b_half / (nodes_half_wing))) * load_factor + W
    Vz = integrate.cumulative_trapezoid(np.flip(Z_forces[nodes_half_wing:]), y_points[nodes_half_wing:], initial=0)
    Vz = Vz[::-1]

    # Vy - axial force diagram because of the strut
    Ax_total_flight = axial_distribution_flight(
        nodes, ac_data
    )  # axial force in the y span, during flight -- this is what we want now cause of cruise
    Ax_total_ground = axial_distribution_ground(nodes, ac_data)  # axial force in the y span, on the ground
    Vy = np.append(Ax_total_flight, [0])

    # Vz - shear force because of the strut
    strut_loc = ac_data["Geometry"]["strut_loc_b/2"]
    n_before_strut = np.rint(nodes_half_wing * strut_loc).astype(int)
    Vz[:n_before_strut] += strut_shear_result()  # [N]

    # add the moment about x
    # print(b, Vz.shape, n, y_points.shape, M.shape)
    Mx = -integrate.cumulative_trapezoid(
        np.flip(Vz * y_points[nodes_half_wing:]), y_points[nodes_half_wing:], initial=0
    )
    Mx = Mx[::-1]
    Mz = -integrate.cumulative_trapezoid(
        np.flip(Vx * y_points[nodes_half_wing:]), y_points[nodes_half_wing:], initial=0
    )
    Mz = Mz[::-1]

    # plt.plot(y_points[nodes_half_wing:], Z_forces[nodes_half_wing:], label="W-L")
    # plt.plot(y_points[nodes_half_wing:], Vz, label="Vz")
    # plt.plot(y_points[nodes_half_wing:], Mx, label="Mx")
    # plt.legend()
    # plt.show()

    # due to lift and weight moment arm (assume both are at c/4)
    My_lw = integrate.cumulative_trapezoid(
        np.flip(Z_forces[nodes_half_wing:] * np.tan(sweep)), y_points[nodes_half_wing:], initial=0
    )
    # Mw = integrate.cumulative_trapezoid(
    #     np.flip((W * b_half / (nodes_half_wing))[nodes_half_wing:] * np.tan(sweep)), y_points[nodes_half_wing:], initial=0
    # )
    My_m = integrate.cumulative_trapezoid(
        np.flip((M * load_factor * b_half / (nodes_half_wing))[nodes_half_wing:]), y_points[nodes_half_wing:], initial=0
    )

    My = (My_lw + My_m)[::-1]

    return Vx, Vy, Vz, Mx, My, Mz


def interpolate_Cl_Cd_Cm(Cl_data, Cdi_data, Cm_data, y_points):

    nodes_data = len(Cl_data[0]["coefficient"])
    ypts_orig = np.linspace(Cl_data[0]["y_span"][0], Cl_data[0]["y_span"][-1], nodes_data)

    for angle in Cl_data.keys():
        Cl_data[angle]["coefficient"] = np.interp(y_points, ypts_orig, Cl_data[angle]["coefficient"])
        Cm_data[angle]["coefficient"] = np.interp(y_points, ypts_orig, Cm_data[angle]["coefficient"])
        Cdi_data[angle]["coefficient"] = np.interp(y_points, ypts_orig, Cdi_data[angle]["coefficient"])

        Cl_data[angle]["y_span"] = y_points
        Cm_data[angle]["y_span"] = y_points
        Cdi_data[angle]["y_span"] = y_points

    return Cl_data, Cdi_data, Cm_data


if __name__ == "__main__":
    # Data
    AoA = -6
    # Sw = 39  # [m2]
    # taper_ratio = 0.4
    # Vcruise = 60  # [m/s]
    # rho = 0.9  # [kg/m3]
    # structuralmass = 5250 / 9.81
    # batterymass_w = 0
    nl = 3.8  # Load Factor
    nl2 = -1.52
    # sweep = 0.157
    altitude = 3000  # [m]
    nodes = 501

    wing_structure_data = WingStructure(aircraft_data, airfoil_shape, nodes)
    chord_dist = wing_structure_data.chord_distribution
    y_points = wing_structure_data.ypts
    # print(Cl_DATA[-10.0].keys())

    Cl_DATA, Cdi_DATA, Cm_DATA = interpolate_Cl_Cd_Cm(Cl_data_wing, Cdi_data_wing, Cm_data_wing, y_points)

    # Cm_DATA = np.interp(y_points, ypts_orig, Cm_DATA[AoA]["coefficient"])
    # Cdi_DATA = np.interp(y_points, ypts_orig, Cdi_DATA[AoA]["coefficient"])

    L_cruise, D_cruise = force_distribution(
        AoA, altitude, aircraft_data["Performance"]["Vc_m/s"], Cl_DATA=Cl_DATA, Cdi_DATA=Cdi_DATA
    )
    W_cruise, W_fuel, idxs, c_between_struts = weight_distribution(chord_dist, ac_data=aircraft_data)

    M_cruise = moment_distribution(
        aircraft_data["Performance"]["Vc_m/s"], altitude, Cm_DATA, AoA, ac_data=aircraft_data
    )

    # nl = origin
    Vx1, Vy1, Vz1, Mx1, My1, Mz1 = InternalLoads(
        L_cruise, W_cruise, D_cruise, M_cruise, nodes, y_points, ac_data=aircraft_data, load_factor=1
    )

    # nl = 3.8
    Vx, Vy, Vz, Mx, My, Mz = InternalLoads(
        L_cruise, W_cruise, D_cruise, M_cruise, nodes, y_points, ac_data=aircraft_data, load_factor=nl
    )

    # nl = -1
    Vx2, Vy2, Vz2, Mx2, My2, Mz2 = InternalLoads(
        L_cruise, W_cruise, D_cruise, M_cruise, nodes, y_points, ac_data=aircraft_data, load_factor=nl2
    )

    # plt.subplot(2, 2, 1)
    # plt.plot(y_points, nl*L_cruise, label='Lift')
    # plt.plot(y_points, -W_cruise, label = 'Weight')
    # plt.title('Lift and weight distribution along Span')
    # plt.xlabel('Spanwise Position [m]')
    # plt.ylabel('Force [N]')
    # plt.legend()
    # plt.xlim(left=0)
    # plt.tight_layout()
    # plt.grid()

    # plt.subplot(2, 2, 1)
    # plt.plot(y_points[nodes // 2 :], Vx, label="n=3.8")
    # plt.plot(y_points[nodes // 2 :], Vx1, label="cruise")
    # plt.plot(y_points[nodes // 2 :], Vx2, label="n=-1")
    # plt.title("Shear Force Vx along Span")
    # plt.xlabel("Spanwise Position [m]")
    # plt.ylabel("Vx [N]")
    # # plt.xlim(left=0)
    # plt.grid()
    # plt.legend()

    plt.subplot(2, 2, 1)
    plt.plot(y_points[nodes // 2 :], Vy, label="n=3.8")
    plt.plot(y_points[nodes // 2 :], Vy1, label="cruise")
    plt.plot(y_points[nodes // 2 :], Vy2, label="n=-1")
    plt.title("Shear Force Vy along Span")
    plt.xlabel("Spanwise Position [m]")
    plt.ylabel("Vy [N]")
    # plt.xlim(left=0)
    plt.grid()
    plt.legend()

    plt.subplot(2, 2, 2)
    plt.plot(y_points[nodes // 2 :], Vz, label="n=3.8")
    plt.plot(y_points[nodes // 2 :], Vz1, label="cruise")
    plt.plot(y_points[nodes // 2 :], Vz2, label="n=-1")
    plt.title("Shear Force Vz along Span")
    plt.xlabel("Spanwise Position [m]")
    plt.ylabel("Vz [N]")
    # plt.xlim(left=0)
    plt.legend()
    plt.grid()

    plt.subplot(2, 2, 3)
    plt.plot(y_points[nodes // 2 :], Mx, label="n=3.8")
    plt.plot(y_points[nodes // 2 :], Mx1, label="cruise")
    plt.plot(y_points[nodes // 2 :], Mx2, label="n=-1")
    plt.title("Bending Moment Mx along Span")
    plt.xlabel("Spanwise Position [m]")
    plt.ylabel("Mx [Nm]")
    # plt.xlim(left=0)
    plt.legend()
    plt.grid()

    # plt.subplot(2, 2, 4)
    # plt.plot(y_points[nodes // 2 :], Mz, label="n=3.8")
    # plt.plot(y_points[nodes // 2 :], Mz1, label="cruise")
    # plt.plot(y_points[nodes // 2 :], Mz2, label="n=-1")
    # plt.title("Torque Mz along Span")
    # plt.xlabel("Spanwise Position [m]")
    # plt.ylabel("Mz [Nm]")
    # # plt.xlim(left=0)
    # plt.legend()
    # plt.grid()

    plt.subplot(2, 2, 4)
    plt.plot(y_points[y_points.shape[0] // 2 :], My, label="n=3.8")
    plt.plot(y_points[y_points.shape[0] // 2 :], My1, label="cruise")
    plt.plot(y_points[y_points.shape[0] // 2 :], My2, label="n=-1")
    plt.title("Torque My along Span")
    plt.xlabel("Spanwise Position [m]")
    plt.ylabel("My [Nm]")
    plt.legend()
    plt.grid()

    plt.tight_layout()
    plt.show()

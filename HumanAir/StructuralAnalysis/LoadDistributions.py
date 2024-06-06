import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from Functions import import_data2


# Define the forces along half span
def chord(Sw, taper_ratio, Cl_DATA, AoA, n):
    b = Cl_DATA[AoA]["y_span"][-1] * 2
    # Generate spanwise coordinate points
    y = np.linspace(Cl_DATA[AoA]["y_span"][0], Cl_DATA[AoA]["y_span"][-1], n)  # n is the number of nodes
    # Calculate the chord distribution
    chord_length = 2 * Sw / (1 + taper_ratio) / b * (1 - (1 - taper_ratio) * np.abs(2 * y / b))
    return chord_length, y


def force_distribution(Cl_DATA, Cdi_DATA, AoA, Sw, V, rho):
    Cl = np.array(Cl_DATA[AoA]["coefficient"])
    Cdi = np.array(Cdi_DATA[AoA]["coefficient"])
    L = Cl * 0.5 * rho * V**2 * Sw  # [N/m]
    D = Cdi * 0.5 * rho * V**2 * Sw  # [N/m]
    return L, D


def weight_distribution(structuralmass, batterymass_w, Cl_DATA, c, AoA):
    c = np.array(c)
    W_ave = (structuralmass + batterymass_w) * 9.81 / c.shape[0]
    W = W_ave * c / ((c[c.shape[0] // 2] + c[0]) / 2)
    return W


def moment_distribution(c, V, rho, Cm_DATA, AoA):
    c = np.array(c)
    return np.array(Cm_DATA[AoA]["coefficient"]) * 0.5 * rho * V**2 * c  # M_cruise


# axial forces distribution - during ground - so only the ends (0.6L)
Vy_Strut = 6670 / np.tan(0.546)  # deg in radian


def axial_distribution_ground(n):
    A_ground = np.zeros(n)
    A_ground[0:12] = np.ones(12) * Vy_Strut
    A_ground[38 - 12 : 38] = np.ones(12) * Vy_Strut
    return A_ground


# axial forces distribution - during flight - so only on the inner (inside 0.4L))
def axial_distribution_flight(n):
    A_cruise = np.ones(n) * Vy_Strut
    A_cruise[0:12] = np.zeros(12)
    A_cruise[38 - 12 : 38] = np.zeros(12)
    return A_cruise


def TestForces(Lcruise, W, Cl_DATA, AoA):
    dy = (Cl_DATA[AoA]["y_span"][-1] - Cl_DATA[AoA]["y_span"][0]) / n
    print("Ltot = ", np.trapz(Lcruise / 9.81, Cl_DATA[AoA]["y_span"]), " kg")
    print("Wtot = ", np.trapz(W / 9.81, Cl_DATA[AoA]["y_span"]), " kg")


def InternalLoads(L, T, W, D, M, n, y_points, Cl_DATA, AoA, sweep):
    b = Cl_DATA[AoA]["y_span"][-1] * 2
    Dtot = T - D  # drag and thrust act on the x axis
    Vx = integrate.cumulative_trapezoid(
        np.flip(Dtot * b / (2 * n))[Dtot.shape[0] // 2 :], y_points[y_points.shape[0] // 2 :]
    )[::-1]
    Vz = integrate.cumulative_trapezoid(
        np.flip((-L + W) * b / (2 * n))[W.shape[0] // 2 :], y_points[y_points.shape[0] // 2 :]
    )[::-1]
    Vx = np.append(Vx, [0])
    Vz = np.append(Vz, [0])

    # Vy - axial force diagram because of the strut
    Ax_total_flight = axial_distribution_flight(
        n
    )  # axial force in the y span, during flight -- this is what we want now cause of cruise
    Ax_total_ground = axial_distribution_ground(n)  # axial force in the y span, on the ground
    Vy = Ax_total_flight

    # Vz - shear force because of the strut
    Vz[7] = -6670  # [N]
    # Vz[38 - 12] = -6670

    # add the moment about x
    # print(b, Vz.shape, n, y_points.shape, M.shape)
    Mx = -integrate.cumulative_trapezoid(np.flip(Vz * b / (2 * (n / 2))), y_points[y_points.shape[0] // 2 :])[::-1]
    Mx = np.append(Mx, [0])
    Mz = -integrate.cumulative_trapezoid(np.flip(Vx * b / (2 * (n / 2))), y_points[y_points.shape[0] // 2 :])[::-1]
    Mz = np.append(Mz, [0])

    # add the torque function
    Ml = []
    Mw = []
    yp = y_points
    count = 2

    # due to lift and weight moment arm
    for i in yp[1:]:
        Mli = IntegrateTorqueFromLift(count, yp, -L, sweep)
        Mwi = IntegrateTorqueFromLift(count, yp, W, sweep)
        Ml.append(Mli)
        Mw.append(Mwi)
        count += 1

    # due to aerodynamic moment:
    Mym = integrate.cumulative_trapezoid(np.flip(M), np.flip(yp))
    My = (np.array(Ml) + np.array(Mw) + np.array(Mym))[::-1]
    My = np.append(My, [0])

    return Vx, Vy, Vz, Mx, My, Mz


def IntegrateTorqueFromLift(c, axis, data, sweep):
    data = np.flip(data)
    data = data[:c]
    axis = np.flip(axis)
    axis = axis[:c]
    data = data * np.tan(sweep) * (axis - axis[c - 1])

    M = np.trapz(data, axis)
    return M


# Data
AoA = -6
Sw = 39  # [m2]
taper_ratio = 0.4
Vcruise = 60  # [m/s]
rho = 0.9  # [kg/m3]
structuralmass = 5250 / 9.81
batterymass_w = 0
T = 0
nl = 3.8  # Load Factor
nl2 = -1.52
sweep = 0.157

# import files
Cl_DATA = import_data2("Cl_DATA.txt")
Cm_DATA = import_data2("Cm_DATA.txt")
Cdi_DATA = import_data2("Cdi_DATA.txt")


n = len(Cl_DATA[AoA]["coefficient"])
# print("number of elements: ",n)

c, y_points = chord(Sw, taper_ratio, Cl_DATA, AoA, n)
L_cruise, D_cruise = force_distribution(Cl_DATA, Cdi_DATA, AoA, Sw, Vcruise, rho)
W_cruise = weight_distribution(structuralmass, batterymass_w, Cl_DATA, c, AoA)
M_cruise = moment_distribution(c, Vcruise, rho, Cm_DATA, AoA)

# nl = 3.8
Vx, Vy, Vz, Mx, My, Mz = InternalLoads(
    nl * L_cruise, T, W_cruise, abs(nl) * D_cruise, nl * M_cruise, n, y_points, Cl_DATA, AoA, sweep
)

# nl = origin
Vx1, Vy1, Vz1, Mx1, My1, Mz1 = InternalLoads(
    L_cruise, T, W_cruise, D_cruise, M_cruise, n, y_points, Cl_DATA, AoA, sweep
)

# nl = -1
Vx2, Vy2, Vz2, Mx2, My2, Mz2 = InternalLoads(
    nl2 * L_cruise, T, W_cruise, abs(nl2) * D_cruise, nl2 * M_cruise, n, y_points, Cl_DATA, AoA, sweep
)

if __name__ == "__main__":
    """
    plt.subplot(2, 2, 1)
    plt.plot(y_points, nl*L_cruise, label='Lift')
    plt.plot(y_points, -W_cruise, label = 'Weight')
    plt.title('Lift and weight distribution along Span')
    plt.xlabel('Spanwise Position [m]')
    plt.ylabel('Force [N]')
    plt.legend()
    plt.xlim(left=0)
    plt.tight_layout()
    plt.grid()
    #plt.show()


    """
    plt.subplot(2, 2, 1)
    plt.plot(y_points[y_points.shape[0] // 2 :], Vx, label="n=3.8")
    plt.plot(y_points[y_points.shape[0] // 2 :], Vx1, label="cruise")
    plt.plot(y_points[y_points.shape[0] // 2 :], Vx2, label="n=-1")
    plt.title("Shear Force Vx along Span")
    plt.xlabel("Spanwise Position [m]")
    plt.ylabel("Vx [N]")
    # plt.xlim(left=0)
    plt.grid()
    plt.legend()

    plt.subplot(2, 2, 2)
    plt.plot(y_points[y_points.shape[0] // 2 :], Vz, label="n=3.8")
    plt.plot(y_points[y_points.shape[0] // 2 :], Vz1, label="cruise")
    plt.plot(y_points[y_points.shape[0] // 2 :], Vz2, label="n=-1")
    plt.title("Shear Force Vz along Span")
    plt.xlabel("Spanwise Position [m]")
    plt.ylabel("Vz [N]")
    # plt.xlim(left=0)
    plt.legend()
    plt.grid()

    plt.subplot(2, 2, 3)
    plt.plot(y_points[y_points.shape[0] // 2 :], Mx, label="n=3.8")
    plt.plot(y_points[y_points.shape[0] // 2 :], Mx1, label="cruise")
    plt.plot(y_points[y_points.shape[0] // 2 :], Mx2, label="n=-1")
    plt.title("Bending Moment Mx along Span")
    plt.xlabel("Spanwise Position [m]")
    plt.ylabel("Mx [Nm]")
    # plt.xlim(left=0)
    plt.legend()
    plt.grid()

    plt.subplot(2, 2, 4)
    plt.plot(y_points[y_points.shape[0] // 2 :], Mz, label="n=3.8")
    plt.plot(y_points[y_points.shape[0] // 2 :], Mz1, label="cruise")
    plt.plot(y_points[y_points.shape[0] // 2 :], Mz2, label="n=-1")
    plt.title("Torque Mz along Span")
    plt.xlabel("Spanwise Position [m]")
    plt.ylabel("Mz [Nm]")
    # plt.xlim(left=0)
    plt.legend()
    plt.grid()

    # plt.subplot(2, 2, 4)
    # plt.plot(y_points[y_points.shape[0]//2:], My, label="n=3.8")
    # plt.plot(y_points[y_points.shape[0]//2:], My1, label="cruise")
    # plt.plot(y_points[y_points.shape[0]//2:], My2, label="n=-1")
    # plt.title("Torque My along Span")
    # plt.xlabel("Spanwise Position [m]")
    # plt.ylabel("My [Nm]")
    # plt.legend()
    # plt.grid()
    # plt.tight_layout()

    plt.show()

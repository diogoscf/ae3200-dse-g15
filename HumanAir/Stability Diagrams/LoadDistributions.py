import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from Functions import import_data2

#define the forces along half span

def chord(Sw, taper_ratio, Cl_DATA, AoA, n):
    b = Cl_DATA[AoA]['y_span'][-1] * 2  
    # Generate spanwise coordinate points
    y = np.linspace(Cl_DATA[AoA]['y_span'][0], Cl_DATA[AoA]['y_span'][-1], n)  # n is the number of nodes
    # Calculate the chord distribution
    chord_length = 2 * Sw / (1 + taper_ratio) / b * (1 - (1 - taper_ratio) * np.abs(2 * y / b))
    return chord_length, y

def force_distribution(Cl_DATA, Cdi_DATA, AoA, c ,V, rho):
    Cl = np.array(Cl_DATA[AoA]['coefficient'])
    Cdi = np.array(Cdi_DATA[AoA]['coefficient'])
    c = np.array(c)
    L = Cl * 0.5 * rho * V**2 * c  # [N/m]
    D = Cdi * 0.5 * rho * V**2 * c  # [N/m]    
    return L, D

def weight_distribution(structuralmass, batterymass_w, Cl_DATA, c, AoA):
    c = np.array(c)
    W_ave = (structuralmass + batterymass_w) / 2 * 9.81 / Cl_DATA[AoA]['y_span'][-1]  # [N/m]
    W = W_ave * c / ((c[0] + c[-1]) / 2)
    return -W

def moment_distribution(c, V, rho, Cm_DATA, AoA):
    c = np.array(c)
    return np.array(Cm_DATA[AoA]['coefficient']) * 0.5 * rho * V**2 * c

def TestForces(Lcruise, W, Cl_DATA, AoA):
    dy = (Cl_DATA[AoA]['y_span'][-1] - Cl_DATA[AoA]['y_span'][0]) / n
    print('Ltot = ', np.trapz(Lcruise / 9.81, Cl_DATA[AoA]['y_span']), ' kg')
    print('Wtot = ', np.trapz(W / 9.81, Cl_DATA[AoA]['y_span']), ' kg')


def InternalLoads(L, T, W, D, M, n, y_points, Cl_DATA, AoA, sweep):
    b = Cl_DATA[AoA]['y_span'][-1] * 2 
    Dtot = T - D  # drag and thrust act on the x axis
    Vx = integrate.cumtrapz(np.flip(Dtot * b / (2 * n)),y_points)[::-1]
    Vz = -integrate.cumtrapz(np.flip((L - W) * b / (2 * n)), np.abs(y_points))[::-1]
    Vx = np.append(Vx, [0])
    Vz = np.append(Vz, [0])

    #add the moment about x 
    Mx = integrate.cumtrapz(np.flip(Vz * b / (2 * n)), y_points)[::-1]
    Mx = np.append(Mx, [0])

    #add the torque function
    Ml = []
    Mw =[]
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
    Mym = integrate.cumtrapz(np.flip(M), np.flip(yp))
    My = (np.array(Ml) + np.array(Mw) + np.array(Mym))[::-1]
    My = np.append(My, [0])

    return Vx, Vz, Mx, My

def IntegrateTorqueFromLift(c, axis, data, sweep):
    data = np.flip(data)
    data = data[:c]
    axis = np.flip(axis)
    axis = axis[:c]
    data = data * np.tan(sweep) * (axis - axis[c-1])

    M = np.trapz(data, axis)
    return M

# Data
AoA = 6
Sw = 39  # [m2]
taper_ratio = 0.4
Vcruise = 60  # [m/s]
rho = 0.9  # [kg/m3]
structuralmass = 2000
batterymass_w = 0
T = 0
sweep = 0.157

# import files
Cl_DATA = import_data2('HumanAir/FuselageSizing/Cl_DATA.txt')
Cm_DATA = import_data2('HumanAir/FuselageSizing/Cm_DATA.txt')
Cdi_DATA = import_data2('HumanAir/FuselageSizing/Cdi_DATA.txt')

n = len(Cl_DATA[AoA]['coefficient'])

c, y_points = chord(Sw, taper_ratio, Cl_DATA, AoA, n)
L_cruise, D_cruise = force_distribution(Cl_DATA, Cdi_DATA, AoA, c ,Vcruise, rho)
W_cruise = weight_distribution(structuralmass, batterymass_w, Cl_DATA, c, AoA)
M_cruise = moment_distribution(c, Vcruise, rho, Cm_DATA, AoA)

Vx, Vz, Mx, Mz = InternalLoads(L_cruise, T, W_cruise, D_cruise, M_cruise, n, y_points, Cl_DATA, AoA, sweep=0.157)



plt.plot(y_points, L_cruise)
plt.plot(y_points, W_cruise)
plt.title('Lift')
plt.xlim(left=0)
plt.tight_layout()
plt.grid()
plt.show()


plt.subplot(2, 2, 1)
plt.plot(y_points, Vx)
plt.title('Shear Force Vx along span')
plt.xlabel('Spanwise Position [m]')
plt.ylabel('Vx [N]')
plt.xlim(left=0)
plt.grid()


plt.subplot(2, 2, 2)
plt.plot(y_points, Vz)
plt.title('Shear Force Vz along span')
plt.xlabel('Spanwise Position [m]')
plt.ylabel('Vz [N]')
plt.xlim(left=0)
plt.grid()


plt.subplot(2, 2, 3)
plt.plot(y_points, Mx)
plt.title('Bending Moment Mx along span')
plt.xlabel('Spanwise Position [m]')
plt.ylabel('Mx [Nm]')
plt.xlim(left=0)
plt.grid()


plt.subplot(2, 2, 4)
plt.plot(y_points, Mz)
plt.title('Torque Mz along span')
plt.xlabel('Spanwise Position [m]')
plt.ylabel('Mz [Nm]')
plt.xlim(left=0)


plt.tight_layout()
plt.grid()
plt.show()

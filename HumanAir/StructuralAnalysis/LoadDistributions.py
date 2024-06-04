import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from Functions import import_data2
import os
import sys

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(script_dir,'..'))
sys.path.append(project_root)

from aircraft_data import aircraft_data
from HumanAir.Vn_Diagrams.loading_diagram import calc_nmax_nmin_manoeuvre

# Define the forces along half span
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
    return W

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
    Vx = integrate.cumtrapz(np.flip(Dtot * b / (2 * n)), y_points)[::-1]
    Vz = integrate.cumtrapz(np.flip((-L + W) * b / (2 * n)), y_points)[::-1]
    Vx = np.append(Vx, [0])
    Vz = np.append(Vz, [0])

    # add the moment about x 
    Mx = -integrate.cumtrapz(np.flip(Vz * b / (2 * n)), y_points)[::-1]
    Mx = np.append(Mx, [0])
    Mz = -integrate.cumtrapz(np.flip(Vx * b / (2 * n)),y_points)[::-1] 
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
    Mym = integrate.cumtrapz(np.flip(M), np.flip(yp))
    My = (np.array(Ml) + np.array(Mw) + np.array(Mym))[::-1]
    My = np.append(My, [0])

    return Vx, Vz, Mx, My, Mz

def IntegrateTorqueFromLift(c, axis, data, sweep):
    data = np.flip(data)
    data = data[:c]
    axis = np.flip(axis)
    axis = axis[:c]
    data = data * np.tan(sweep) * (axis - axis[c-1])

    M = np.trapz(data, axis)
    return M

def load_distribution_diagram(ac_data = aircraft_data):
    # Data
    AoA = -6
    Sw = ac_data['Aero']['S_Wing']
    taper_ratio = ac_data['Aero']['Taper_Wing']
    Vcruise = ac_data['Performance']['Vc_m/s']  # [m/s]
    rho = ac_data['Performance']['rho_kg/m3'] # [kg/m3]
    structuralmass = ac_data['CL2Weight']['Wing Weight'] / 9.81
    batterymass_w = 0
    T = 0
    nmax, nmin = calc_nmax_nmin_manoeuvre(ac_data['Weights']['MTOW_N'])
    nl = nmax # Maximum Load Factor
    nl2 = nmin # Minimum Load Factor
    sweep = np.deg2rad(ac_data['Aero']['HalfChordSweep_Wing_deg']) # half chord sweep angle of the wing

    # import files
    Cl_DATA = import_data2('Cl_DATA.txt')
    Cm_DATA = import_data2('Cm_DATA.txt')
    Cdi_DATA = import_data2('Cdi_DATA.txt')

    n = len(Cl_DATA[AoA]['coefficient'])

    c, y_points = chord(Sw, taper_ratio, Cl_DATA, AoA, n)
    L_cruise, D_cruise = force_distribution(Cl_DATA, Cdi_DATA, AoA, c ,Vcruise, rho)
    W_cruise = weight_distribution(structuralmass, batterymass_w, Cl_DATA, c, AoA)
    M_cruise = moment_distribution(c, Vcruise, rho, Cm_DATA, AoA)

    #print(c)
    # nl = 3.8
    Vx, Vz, Mx, My, Mz = InternalLoads(nl*L_cruise, T, W_cruise, abs(nl)*D_cruise, nl*M_cruise, n, y_points, Cl_DATA, AoA, sweep)

    # nl = origin
    Vx1, Vz1, Mx1, My1, Mz1 = InternalLoads(L_cruise, T, W_cruise, D_cruise, M_cruise, n, y_points, Cl_DATA, AoA, sweep)

    # nl = -1
    Vx2, Vz2, Mx2, My2, Mz2 = InternalLoads(nl2*L_cruise, T, W_cruise, abs(nl2)*D_cruise, nl2*M_cruise, n, y_points, Cl_DATA, AoA, sweep)


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



    '''
    plt.subplot(2, 2, 1)
    plt.plot(y_points, Vx, label='n=3.8')
    plt.plot(y_points, Vx1, label='cruise')
    plt.plot(y_points, Vx2, label='n=-1')
    plt.title('Shear Force Vx along Span')
    plt.xlabel('Spanwise Position [m]')
    plt.ylabel('Vx [N]')
    plt.xlim(left=0)
    plt.grid()
    plt.legend()
    '''

    plt.subplot(2, 2, 2)
    plt.plot(y_points, Vz, label=f'n={round(nmax, 2)}')
    plt.plot(y_points, Vz1, label='cruise')
    plt.plot(y_points, Vz2, label=f'n={round(nmin, 2)}')
    plt.title('Shear Force Vz along Span')
    plt.xlabel('Spanwise Position [m]')
    plt.ylabel('Vz [N]')
    plt.xlim(left=0)
    plt.legend()
    plt.grid()

    plt.subplot(2, 2, 3)
    plt.plot(y_points, Mx, label=f'n={round(nmax, 2)}')
    plt.plot(y_points, Mx1, label='cruise')
    plt.plot(y_points, Mx2, label=f'n={round(nmin, 2)}')
    plt.title('Bending Moment Mx along Span')
    plt.xlabel('Spanwise Position [m]')
    plt.ylabel('Mx [Nm]')
    plt.xlim(left=0)
    plt.legend()
    plt.grid()


    '''
    plt.subplot(2, 2, 4)
    plt.plot(y_points, Mz, label='n=3.8')
    plt.plot(y_points, Mz1, label='cruise')
    plt.plot(y_points, Mz2, label='n=-1')
    plt.title('Torque Mz along Span')
    plt.xlabel('Spanwise Position [m]')
    plt.ylabel('Mz [Nm]')
    plt.xlim(left=0)
    plt.legend()
    plt.grid()
    '''

    plt.subplot(2,2,4)
    plt.plot(y_points, My, label=f'n={round(nmax,2)}')
    plt.plot(y_points, My1, label='cruise')
    plt.plot(y_points, My2, label=f'n={round(nmin,2)}')
    plt.title('Torque My along Span')
    plt.xlabel('Spanwise Position [m]')
    plt.ylabel('My [Nm]')
    plt.legend()
    plt.grid()
    plt.tight_layout()

    plt.show()

if __name__ == '__main__':
    load_distribution_diagram(ac_data=aircraft_data)
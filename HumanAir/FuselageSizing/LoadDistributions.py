import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import integrate

n = 100 # resolution of data along half span 
#Mass in the wing
Batterymass_f = 0 #kg from 2nd weight estimation, say 0 if on the wing
Batterymass_w = 0 #kg if the battery is on the wing
Structuralmass = 0 #kg from 2nd weight estimation 
Enginemass = 0 #kg from propulsion team 
#Propellermass = 0
#Landinggearmass = 0

#Propeller Thrust
#Thrust_TO = 0 #N 
#Thrust_Cruise = 0 #N
T=0 #no engines on the wing

#Worst-case load factors from Vn diagram
Loadfactor = 0 

#Aerodynamic loads -- comes from the xflr5 files
def import_data(n):
    Excel_cruise = 0
    Excel_TO = 0
    data_cruise = pd.read_csv(Excel_cruise, index_col=None)
    data_TO = pd.read_csv(Excel_cruise, index_col=None)
    data_cruise = data_cruise.to_dict(orient='list')
    data_TO = data_TO.to_dict(orient="list")
    y_old = data_cruise['y-span']

    y_points = np.linspace(data_cruise['y-span'][0], data_cruise['y-span'][-1],n)
    for i in data_cruise.keys():
        data_cruise[i]=np.interp(y_points,y_old, data_cruise[i])
    for i in data_TO.keys():
        data_TO[i]=np.interp(y_points,y_old,data_TO[i])

    return data_cruise, data_TO, y_points 

#define the forces along half span
def forces(Excel_cruise, Excel_TO, y_points,n):
    #Lift and drag at cruise 
    Cl = Excel_cruise['Cl']
    Cdi = Excel_cruise['ICd'] #see if this is the case
    c = Excel_cruise['Chord']
    V = 60 #m/s cruise speed 
    rho = 0 #kg/m^3 density at altitude of cruise 
    dy = (Excel_cruise['y-span'][-1]-Excel_cruise['y-span'][0])/n
    Lcruise = Cl*0.5*rho*V**2*c #N/m
    Dcruise = Cdi*0.5*rho*V**2*c#N/m

    #Lift and drag at Take-off 
    Cl = Excel_TO['Cl']
    Cdi = Excel_TO['ICd']
    c = Excel_TO['Chord']
    V = 38 #m/s take off speed 
    rho = 1.225 #kg/m^3 density at sea level
    dy = (Excel_TO['y-span'][-1]-Excel_TO['y-span'][0])/n
    Ltakeoff = Cl*0.5*rho*V**2*c #N/m
    Dtakeoff = Cdi*0.5*rho*V**2*c#N/m

    #Weight distribution along span
    W_ave = (Structuralmass + Batterymass_w)/2*9.81/Excel_TO['y-span'][-1] #N/m
    W = W_ave*Excel_TO['Chord']/((Excel_TO['Chord'][0]+Excel_TO['Chord'][-1])/2) #ensuring that the weight distribution accurately reflects the variation in chord length along the wing span

    #Thrust force -- engines are not on the wing

    #Moment Distribution 
    V = 0 #m/s cruise speed
    rho = 0 #kg/m^3 density at cruise altitude
    M_cruise = Excel_cruise['CmAirf@chord/4']*0.5*rho*V**2*Excel_cruise['Chord']

    V = 0 #m/s take off speed
    rho = 1.225 #kg/m^3 density at sea level 
    M_TO = Excel_TO['CmAirf@chord/4']*0.5*rho*V**2*Excel_TO['Chord']

    return Lcruise, Ltakeoff, W, M_cruise, M_TO, Dcruise, Dtakeoff

def TestForces(Lcruise, Ltakeoff, Tcruise, W, Excel_TO):
    dy = (Excel_TO['y-span'][-1]- Excel_TO['y-span'][0]) / n
    c = Excel_TO['Chord']
    print('Ltot = ',np.trapz(Lcruise/9.81, Excel_TO['y-span']),' kg')
    print('Wtot = ',np.trapz(W/9.81,Excel_TO['y-span']),' kg')

sweep = 0 #add sweep angle 
def InternalLoads(L,T,W,D,M,n,y_points,data,sweep):
    b = data['y-span'][-1]*2
    Dtot = T-D # drag and thrust act on the x axis
    Vx = integrate.cumtrapz(np.flip(Dtot * b / (2 * n)))[::-1]
    Vz = integrate.cumtrapz(np.flip((-L + W) * b / (2 * n)))[::-1]
    Vx = np.append(Vx,[0])
    Vz = np.append(Vz,[0])

    #add the moment about x 
    Mx = integrate.cumtrapz(np.flip(Vz * b / (2 * n)))[::-1]
    Mx = np.append(Mx,[0])

    #add the torque function
    Mz = integrate.cumtrapz(np.flip(Vz*np.tan(sweep)*b/(2*n)))[::-1]
    #dy = (Excel_TO['y-span'][-1]- Excel_TO['y-span'][0]) / n
    #Mz = integrate.cumtrapz(Vz*np.tan(sweep)*(y_points[1:]-y_points[:-1]),dx=dy)
    Mz = np.append(Mz,[0])

    return Vx, Vz, Mx, Mz

Excel_cruise, Excel_TO, y_points = import_data(n)
Lcruise, Ltakeoff, W, M_cruise, M_TO, Dcruise, Dtakeoff = forces(Excel_cruise, Excel_TO, y_points, n)
Vx, Vz, Mx, Mz = InternalLoads(Lcruise, T, W, Dcruise, M_cruise, n, y_points, Excel_cruise) #for the cruise

plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(y_points, Vx)
plt.title('Shear Force Vx along span')
plt.xlabel('Spanwise Position (m)')
plt.ylabel('Vx (N)')

plt.subplot(2, 2, 2)
plt.plot(y_points, Vz)
plt.title('Shear Force Vz along span')
plt.xlabel('Spanwise Position (m)')
plt.ylabel('Vz (N)')

plt.subplot(2, 2, 3)
plt.plot(y_points, Mx)
plt.title('Bending Moment Mx along span')
plt.xlabel('Spanwise Position (m)')
plt.ylabel('Mx (Nm)')

plt.subplot(2, 2, 4)
plt.plot(y_points, Mz)
plt.title('Torque Mz along span')
plt.xlabel('Spanwise Position (m)')
plt.ylabel('Mz (Nm)')

plt.tight_layout()
plt.show()
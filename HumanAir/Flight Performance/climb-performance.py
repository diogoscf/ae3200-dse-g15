import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint

import helper
import aircraft


def plot_climb_rate(acf):
    h_max = 7000
    
    W = acf.W_MTO
    dT = 0    
    
    h_list = np.arange(0, h_max, 10)
    RC_max = np.zeros(len(h_list))
    P_a = np.zeros(len(h_list))
    P_r_min = np.zeros(len(h_list))

    service_ceiling = 0
    
    for i, h in enumerate(h_list):
        RC_max[i] = acf.RC_max(W, h, dT)
        P_a[i] = acf.P_a(h, dT)
        P_r_min[i] = acf.P_r_min(W, h, dT)

        # service ceiling defined as altitude with 100ft/min (0.508m/s) max climbrate
        # FAA Pilots Handbook of Aeronautical Knowledge, but it is not a fixed definition I think
        if RC_max[i] < 0.508 and service_ceiling < 0.01 and i > 0:
            service_ceiling = h_list[i-1]
    
    print(f"Sea level MTOW max climb rate: {RC_max[0]:.3f} m/s")
    print(f"Service ceiling: {service_ceiling:.0f} m")
    
    plt.figure(figsize=(10,7))
    plt.plot(RC_max, h_list, label="Maximum Climb Rate")
    plt.plot(P_a/100000, h_list, color="g", label="Power available /100000")
    plt.plot(P_r_min/100000, h_list, color="r", label="Minimum required power /100000")
    plt.axhline(y=service_ceiling, color='grey', label="Service Ceiling")        
    plt.xlabel(r"Climb Rate (m/s)")
    plt.ylabel("Altitude (m)")
    plt.ylim(0,h_max)
    plt.xlim(0,RC_max[0]*1.05)
    plt.legend()
    #plt.subplots_adjust(right=0.75)
    #plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.savefig("climb_rate.svg")
    plt.show()
    
def plot_climb_gradient(acf):
    h_max = 7000
    
    W = acf.W_MTO
    dT = 0    
    
    h_list = np.arange(0, h_max, 10)
    max_gradient = np.zeros(len(h_list))
    T = np.zeros(len(h_list))
    D = np.zeros(len(h_list))
    
    for i, h in enumerate(h_list):
        max_gradient[i] = acf.climb_angle_max(W, h, dT)
        T[i] = acf.P_a(h, dT)/acf.V_Dmin(W, h, dT)
        D[i] = acf.D_min(W)

    print(f"Sea level MTOW max climb gradient: {max_gradient[0]:.3f} %")
    
    plt.figure(figsize=(10,7))
    plt.plot(max_gradient, h_list, label="Maximum Climb Gradient")
    plt.plot(T/1000, h_list, color="g", label="Thrust /1000")
    plt.plot(D/1000, h_list, color="r", label="Drag /1000")
    plt.xlabel(r"Climb Gradient (%)")
    plt.ylabel("Altitude (m)")
    plt.ylim(0,h_max)
    plt.xlim(0,max_gradient[0]*1.05)
    plt.legend()
    #plt.subplots_adjust(right=0.75)
    #plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.savefig("climb_gradient.svg")
    plt.show()

if __name__ == "__main__":
    acf = aircraft.Aircraft()
    plot_climb_rate(acf)
    plot_climb_gradient(acf)
    
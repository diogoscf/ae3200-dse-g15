import numpy as np
import aircraft

import matplotlib.pyplot as plt

def plot_variation():
    
    plt.figure(figsize=(6,3))
    
    for i in [1,2,4]:
        if i == 2:
            # ***** design *****
            elevation   = 750
            temp_offset = 18
            slope       = 0
            surface     = "grass"
            label = "750m ISA+18"
        if i == 1:
            # **** sea level *****
            elevation   = 0
            temp_offset = 0
            slope       = 0
            surface     = "grass"
            label = "0m ISA+0"
        if i == 4:
            # **** 1500m *****
            elevation   = 1500
            temp_offset = 18
            slope       = 0
            surface     = "grass"
            label = "1500m ISA+18"
        if i == 3:
            # **** 1000m *****
            elevation   = 1000
            temp_offset = 18
            slope       = 0
            surface     = "grass"
            label = "1000m ISA+18"
        if i == 5:
            # **** 1000m *****
            elevation   = 2000
            temp_offset = 18
            slope       = 0
            surface     = "grass"
            label = "2000m ISA+18"
        
        acf = aircraft.Aircraft()
        W_list = np.linspace(acf.W_OE, acf.W_MTO, 25)
        
        
        acf = aircraft.Aircraft()
        
        TO_dist = []
        
        for W in W_list:
            TO_dist.append(acf.takeoff_ground_run(W, elevation, temp_offset, slope, surface)[0])
            
        plt.plot(W_list/9.80665, TO_dist, label=label)
            
    #plt.ylim(300, 900)
    plt.ylabel("Required runway length [m]")
    plt.xlabel("Gross takeoff weight [kg]")
    plt.legend()
    plt.grid()
    plt.savefig("plots/takeoff.svg")
    plt.show()
        
    print(f"OE: {acf.W_OE}")
    print(f"MTO: {acf.W_MTO}")

              
if __name__ == "__main__":
     plot_variation()

# we can reduce wing area to 28m^2 if we
# use 50% payload and 50% fuel for 750m ISA+15 takeoff
# decrease MTOW max ROC to 4 m/s
# increase MTOW stall speed to 31 m/s (CS23)


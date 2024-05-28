import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
sys.path.append(project_root)

from HumanAir.LoadingDiagram.Parameters import Parameters_ConvNoCanard as p
import HumanAir.LoadingDiagram.Equations as eq
import HumanAir.LoadingDiagram.Plotting as plot


class WP_WS:
    # ========== 0: Initialise WS ==========
    def __init__(self):
        self.step = 0.1
        self.final_value_WS = 1500
        self.WS=np.arange(0.0001,self.final_value_WS,self.step)
        self.ylst=np.linspace(0,1,int(self.final_value_WS/self.step))

        #ReferenceWS=[449.0307667, 602.1288757, 643.1831172, 725.2916003, 725.2916003, 578.0625962,
        #789.2055899, 780.3381114, 930.8329016, 925.2755192, 1274.01538, 1286.972761, 1236.320571,
        #987.5590745, 1321.946577, 1662.527136, 1437.958113, 1871.550807, 1837.537433]

        #ReferenceWP=[0.089508018, 0.082049016, 0.080130987, 0.068752535, 0.068752535, 0.058831561,
        #0.066557244, 0.069219534, 0.088135561, 0.078032631, 0.06465638, 0.064313168, 0.062743365,
        #0.066302235, 0.055426119, 0.066832653, 0.065302355, 0.093120238, 0.042435741]

#       ========== 1: Calculate Values ==========
        self.Stallspeed_x=eq.Stallspeedx(p.h_Land, p.V_stall, p.Clmax_Land) #Choose if we want to use Cl with clean config or landing config
        self.Takeoff_y=eq.Takeoff(p.TOP, self.WS, p.h_TO, p.Clmax_TO)
        self.Landing_x=eq.Landingx(p.Clmax_Land, p.h_Land, p.s_land, p.f)
        self.Cruise_y=eq.Cruise(p.eta_p, p.h_Cruise, p.Cdo, p.V_cruise, self.WS, p.A, p.e, p.CruisePower, p.CruiseWeight)
        self.Climbrate_y=eq.Climbrate(p.eta_p, p.A, p.e, p.h_TO, p.Cdo, p.climbrate, self.WS)
        self.Climbgradient_y=eq.Climbgradient(p.eta_p, self.WS, p.climbrate, p.V_climb, p.A, p.e, p.Clmax_clean, p.Cl_SafetyFactor, p.Cdo, p.h_TO)
        #Manouvering_y=eq.Manouvering(p.Cdo, p.h_Cruise, p.V_cruise, WS, p.nmax, p.A, p.e, p.eta_p)

        # Convert the x-values to y-values
        self.Landing_y = self.Landing_x * self.ylst
        self.Stallspeed_y = self.Stallspeed_x * self.ylst

    def calculate_optimal_point(self):
        # Create a matrix where each column is one of the y-values
        y_matrix = np.vstack((self.Stallspeed_y, self.Landing_y, self.Takeoff_y, self.Cruise_y, self.Climbrate_y, self.Climbgradient_y))

        # Find the minimum values across all curves for each WS
        min_envelope = np.min(y_matrix, axis=0)


        optimal_WP = np.max(min_envelope)
        optimal_WS = self.WS[0]

        # change this value to tune it but just god knows which is the optimal one
        WP_tolerance = 1e-3  # Define the tolerance for W/P change
        WS_tolerance = 0.1 * self.WS[np.where(min_envelope == optimal_WP)[0][0]]   # Define the tolerance for W/S change


        for y in range(len(min_envelope)-1):
            if np.abs(min_envelope[y] - optimal_WP) < WP_tolerance: # check if the numbers are increasing and update the optimal values
                if np.abs(self.WS[y] - optimal_WS) > WS_tolerance:
                    optimal_WP = min_envelope[y]
                    optimal_WS = self.WS[y]

            elif min_envelope[y] - optimal_WP > WP_tolerance: # check if the values for W/P are decreasing and update the optimal values
                optimal_WP = min_envelope[y]
                optimal_WS = self.WS[y]

        #print(np.max(min_envelope))
        # print(f"Optimal W/S = {optimal_WS}, Optimal W/P = {optimal_WP}")
        return optimal_WP, optimal_WS

# index=np.where(np.abs(Cruise_y-Takeoff_y)<0.0005)[0][0]
# print("INDEX:")
# print(index)
# print("W/S = ", WS[index])
# print("W/P = ", Cruise_y[index])
#

#  ========== 2: Plot Lines =========="""
    def plot(self, saving = None):
        plt.figure(figsize=(10,7))

        plot.Plotx(self.Stallspeed_x, self.ylst, "Stall Speed Requirement")
        plot.Plotx(self.Landing_x, self.ylst, "Landing Requirement")
        plot.Ploty(self.WS, self.Takeoff_y, "Take-off requirement")
        plot.Ploty(self.WS, self.Cruise_y, "Cruise Requirement")
        plot.Ploty(self.WS, self.Climbrate_y, "Climbrate Requirement")
        plot.Ploty(self.WS, self.Climbgradient_y, "Climbgradient Requirement")


        plt.fill_between(self.WS, self.Takeoff_y, 1, color='red', alpha=.1)
        plt.fill_between(self.WS, self.Cruise_y, 1, color='red', alpha=.1)
        plt.fill_between(self.WS, self.Climbrate_y, 1, color='red', alpha=.1)
        plt.fill_between(self.WS, self.Climbgradient_y, 1, color='red', alpha=.1)
        plt.fill_betweenx(self.ylst, self.Stallspeed_x, 2000, color='red', alpha=.1)
        plt.fill_betweenx(self.ylst, self.Landing_x, 2000, color='red', alpha=.1)
        #plt.scatter(ReferenceWS, ReferenceWP, color='orange', alpha=0.5, label="Reference Aircraft")
        #plt.scatter(WS[index],Cruise_y[index], label="Chosen Design Point", s=100, color='red')
        #plt.scatter(1040.95, 0.076651773, label="Cessna 206", s=100, color='blue')
        plt.xlabel(r"W/S (N/$m^2$)")
        plt.ylabel("W/P (N/W)")
        plt.ylim(0,0.4)
        plt.xlim(0,1500)
        plt.subplots_adjust(right=0.75)
        plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
        plt.title(p.name)

        if saving:
            plt.savefig("Performance_FlyingWing.svg")

        plt.show()
        
        
if __name__ == "__main__":
    wp = WP_WS()
    wp.plot()
    print(wp.calculate_optimal_point())
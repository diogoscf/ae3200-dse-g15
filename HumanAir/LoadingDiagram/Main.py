import numpy as np
import matplotlib.pyplot as plt
import Parameters.Parameters_ConvNoCanard as p
import Equations as eq
import Plotting as plot

"""========== 0: Initialise WS =========="""
WS=np.arange(0,1500,2)
ylst=np.linspace(0,1,500)

ReferenceWS=[449.0307667, 602.1288757, 643.1831172, 725.2916003, 725.2916003, 578.0625962,
789.2055899, 780.3381114, 930.8329016, 925.2755192, 1274.01538, 1286.972761, 1236.320571,
987.5590745, 1321.946577, 1662.527136, 1437.958113, 1871.550807, 1837.537433]

ReferenceWP=[0.089508018, 0.082049016, 0.080130987, 0.068752535, 0.068752535, 0.058831561,
0.066557244, 0.069219534, 0.088135561, 0.078032631, 0.06465638, 0.064313168, 0.062743365,
0.066302235, 0.055426119, 0.066832653, 0.065302355, 0.093120238, 0.042435741]

"""========== 1: Calculate Values =========="""
Stallspeed_x=eq.Stallspeedx(p.h_Land, p.V_stall, p.Clmax_Land) #Choose if we want to use Cl with clean config or landing config
Takeoff_y=eq.Takeoff(p.TOP, WS, p.h_TO, p.Clmax_TO)
Landing_x=eq.Landingx(p.Clmax_Land, p.h_Land, p.s_land, p.f)
Cruise_y=eq.Cruise(p.eta_p, p.h_Cruise, p.Cdo, p.V_cruise, WS, p.A, p.e, p.CruisePower, p.CruiseWeight)
Climbrate_y=eq.Climbrate(p.eta_p, p.A, p.e, p.h_TO, p.Cdo, p.climbrate, WS)
Climbgradient_y=eq.Climbgradient(p.eta_p, WS, p.climbrate, p.V_climb, p.A, p.e, p.Clmax_clean, p.Cl_SafetyFactor, p.Cdo, p.h_TO)
#Manouvering_y=eq.Manouvering(p.Cdo, p.h_Cruise, p.V_cruise, WS, p.nmax, p.A, p.e, p.eta_p)


index=np.where(np.abs(Cruise_y-Takeoff_y)<0.0005)[0][0]
print("W/S = ", WS[index])
print("W/P = ", Cruise_y[index])


"""========== 2: Plot Lines =========="""
plt.figure(figsize=(10,7))
plot.Plotx(Stallspeed_x, ylst, "Stall Speed Requirement")
plot.Plotx(Landing_x, ylst, "Landing Requirement")
plot.Ploty(WS, Takeoff_y, "Take-off requirement")
plot.Ploty(WS, Cruise_y, "Cruise Requirement")
plot.Ploty(WS, Climbrate_y, "Climbrate Requirement")
plot.Ploty(WS, Climbgradient_y, "Climbgradient Requirement")
plt.fill_between(WS, Takeoff_y, 1, color='red', alpha=.1)
plt.fill_between(WS, Cruise_y, 1, color='red', alpha=.1)
plt.fill_between(WS, Climbrate_y, 1, color='red', alpha=.1)
plt.fill_between(WS, Climbgradient_y, 1, color='red', alpha=.1)
plt.fill_betweenx(ylst, Stallspeed_x, 2000, color='red', alpha=.1)
plt.fill_betweenx(ylst, Landing_x, 2000, color='red', alpha=.1)
plt.scatter(ReferenceWS, ReferenceWP, color='orange', alpha=0.5, label="Reference Aircraft")
plt.scatter(WS[index],Cruise_y[index], label="Chosen Design Point", s=100, color='red')
plt.scatter(1040.95, 0.076651773, label="Cessna 206", s=100, color='blue')
plt.xlabel(r"W/S (N/$m^2$)")
plt.ylabel("W/P (N/W)")
plt.ylim(0,0.4)
plt.xlim(0,1500)
plt.subplots_adjust(right=0.75)
plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
plt.title(p.name)
plt.savefig("Performance_FlyingWing.svg")
plt.show()
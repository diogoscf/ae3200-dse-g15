import numpy as np
import matplotlib.pyplot as plt
import Parameters as p
import Equations as eq
import Plotting as plot

"""========== 0: Initialise WS =========="""
WS=np.arange(0,2000,2)
ylst=np.linspace(0,1,500)

"""========== 1: Calculate Values =========="""
Stallspeed_x=eq.Stallspeedx(p.h_Land, p.V_stall, p.Clmax_clean) #Choose if we want to use Cl with clean config or landing config
Takeoff_y=eq.Takeoff(p.TOP, WS, p.h_TO, p.Clmax_TO)
Landing_x=eq.Landingx(p.Clmax_Land, p.h_Land, p.s_land, p.f)
Cruise_y=eq.Cruise(p.eta_p, p.h_Cruise, p.Cdo, p.V_cruise, WS, p.A, p.e, p.CruisePower, p.CruiseWeight)
Climbrate_y=eq.Climbrate(p.eta_p, p.A, p.e, p.h_TO, p.Cdo, p.climbrate, WS)
Climbgradient_y=eq.Climbgradient(p.eta_p, WS, p.climbrate, p.V_climb, p.A, p.e, p.Clmax_clean, p.Cl_SafetyFactor, p.Cdo, p.h_TO)
#Manouvering_y=eq.Manouvering(p.Cdo, p.h_Cruise, p.V_cruise, WS, p.nmax, p.A, p.e, p.eta_p)


"""========== 2: Plot Lines =========="""
plt.figure()
plot.Plotx(Stallspeed_x, ylst, "Stall Speed Requirement")
plot.Plotx(Landing_x, ylst, "Landing Requirement")
plot.Ploty(WS, Takeoff_y, "Take-off requirement")
plot.Ploty(WS, Cruise_y, "Cruise Requirement")
plot.Ploty(WS, Climbrate_y, "Climbrate Requirement")
#plot.Ploty(WS, Climbgradient_y, "Climbgradient Requirement")
plt.legend()
plt.fill_between(WS, Takeoff_y, 1, color='red', alpha=.1)
plt.fill_between(WS, Cruise_y, 1, color='red', alpha=.1)
plt.fill_between(WS, Climbrate_y, 1, color='red', alpha=.1)
plt.fill_between(WS, Climbgradient_y, 1, color='red', alpha=.1)
plt.fill_betweenx(ylst, Stallspeed_x, 2000, color='red', alpha=.1)
plt.fill_betweenx(ylst, Landing_x, 2000, color='red', alpha=.1)
plt.xlabel("W/S (N/m^2)")
plt.ylabel("W/P (N/W)")
plt.ylim(0,0.4)
plt.xlim(0,1500)
plt.show()
import numpy as np
import aircraft


def vary_area():
        
    
    for S in np.arange(20, 40, step=1):
        
        acf = aircraft.Aircraft()
        
        reduced_weight = acf.W_OE + 0.5*acf.W_pl_no_pilot + 0.5*acf.W_MF
        
        acf.S = S
        
        climbrate = acf.RC_max(acf.W_MTO, 0, 0)
        gradient  = acf.climb_slope_max(acf.W_MTO, 0, 0)
        V_S0      = acf.stall_speed(acf.W_MTO, 0, 0)
        TO_dist   = acf.takeoff_ground_run(reduced_weight, 750, 18, 0, "grass")[0]
        land_dist = acf.landing_ground_distance(reduced_weight, 750, 18, 0, "grass", reversible_pitch=True)
        
        pass_cr = "pass" if climbrate >= 5 else "FAIL"
        pass_gr = "pass" if gradient >= 8.3 else "FAIL"
        pass_vs = "pass" if V_S0 <= 25 else "FAIL"
        pass_to = "pass" if TO_dist <= 500 else "FAIL"
        pass_ld = "pass" if land_dist <= 500 else "FAIL"
    
        
        print(f"""
              ****************** S = {S} ***********************
              Climb rate     = {climbrate:>8.2f} m/s    {pass_cr}
              Climb gradient = {gradient:>8.2f} %      {pass_gr}
              V_S0           = {V_S0:>8.2f} m/s    {pass_vs}
              Takeoff run    = {TO_dist:>8.0f} m      {pass_to}
              Landing roll   = {land_dist:>8.0f} m      {pass_ld}
              *************************************************
              """)
              
def vary_weight():
    
    acf = aircraft.Aircraft()
    
    acf.S = 28    

    for W in np.arange(acf.W_OE, acf.W_MTO+500, step=500):
        
        TO_dist   = acf.takeoff_ground_run(W, 750, 18, 0, "grass")[0]
        land_dist = acf.landing_ground_distance(W, 750, 18, 0, "grass", reversible_pitch=True)
        
        pass_to = "pass" if TO_dist <= 500 else "FAIL"
        pass_ld = "pass" if land_dist <= 500 else "FAIL"
    
        fuel_perc = (W - acf.W_OE - acf.W_pl_no_pilot*.5) / acf.W_MF * 100
        
        print(f"""
              ***** W = {W:.0f} *** Design payload 540 kg *******
              Fuel           = {fuel_perc:>8.2f} %
              Takeoff run    = {TO_dist:>8.0f} m      {pass_to}
              Landing roll   = {land_dist:>8.0f} m      {pass_ld}
              *************************************************
              """)
              
if __name__ == "__main__":
    vary_area()
    #vary_weight()

# we can reduce wing area to 28m^2 if we
# use 50% payload and 50% fuel for 750m ISA+15 takeoff
# decrease MTOW max ROC to 4 m/s
# increase MTOW stall speed to 31 m/s (CS23)

# TODO: met niels over takeoff en landing speeds hebben
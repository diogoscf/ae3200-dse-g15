import numpy as np
import aircraft

acf = aircraft.Aircraft()

reduced_weight = acf.W_MTO - (acf.W_pl_max-acf.W_pl_des)

for S in np.arange(20, 40, step=1):
    acf.S = S
    climbrate = acf.RC_max(acf.W_MTO, 0, 0)
    gradient  = acf.climb_slope_max(acf.W_MTO, 0, 0)
    V_S1      = acf.stall_speed(acf.W_MTO, 0, 0)
    TO_dist   = acf.takeoff_ground_run(reduced_weight, 750, 18, 0, "grass")[0]
    land_dist = acf.landing_ground_distance(reduced_weight, 750, 18, 0, "grass")

    
    print(f"""
          ****************** S = {S} **********************
          Climb rate -   = {climbrate:.2f} m/s
          Climb gradient = {gradient:.2f} %
          V_S1           = {V_S1:.2f} m/s
          Takeoff run    = {TO_dist:.0f} m
          Landing roll   = {land_dist:.0f} m
          *************************************************
          """)
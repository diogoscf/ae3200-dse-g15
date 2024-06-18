# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 10:41:53 2024

@author: Alex
"""

import mission_performance

num_legs = 4

tot_dist_nm = 600
leg_dist_m = tot_dist_nm * 1852 / num_legs

airfield_elevation = 750
cruise_altitude = 3000
loiter_altitude = 1200
temp_offset = 18
climb_rate = 2.5
cruise_speed = 60
loiter_duration = 115*60 # [s]

macf = mission_performance.MAircraft(airfield_elevation)
macf.dT = temp_offset


#
# flight
#

for i in range(num_legs):
    # takeoff
    macf.take_off(airfield_elevation, "grass", electric=i==0)
    
    # acceleration to ideal CL for climbing
    V_start_climb = macf.acf.V_Prmin(macf.W, macf.h, macf.dT)
    macf.fly_accelerate_const_alt(V_start_climb, electric=i==0)
    
    # climb
    macf.fly_const_CL_const_climbrate(climb_rate, cruise_altitude, electric=i==0)
    
    # accelerate to cruise speed
    macf.fly_accelerate_const_alt(cruise_speed, electric=i==0)
    
    # cruise
    macf.fly_const_V_const_alt(leg_dist_m, electric=i==0)
    
    if i < num_legs-1:
        # descent
        macf.fly_const_V_descent(airfield_elevation)
    
        # land
        macf.land()
    else:
        # descent to loiter alt
        macf.fly_const_V_descent(loiter_altitude)
        
        # decelerate or accelerate to loiter speed
        V_loiter = macf.acf.V_Prmin(macf.W, macf.h, macf.dT)
        macf.fly_accelerate_const_alt(V_loiter)
        
        # loiter
        macf.fly_const_CL_const_alt(loiter_duration)
        
        # descent
        macf.fly_const_V_descent(airfield_elevation)
        
        # land
        macf.land()


# plot and print

macf.plot_flight()
macf.print_energy_level()

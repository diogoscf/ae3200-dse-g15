# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 10:41:53 2024

@author: Alex
"""

# TODO: add fuel reserve requirement to report?

import mission_performance

better_batteries = False

better_battery_factor = 685 / 350

num_legs = 4

range_adjust = 0.725 # since distance is also covered in climb and descent we lower the cruise range
tot_dist_nm = 600 * range_adjust
leg_dist_m = tot_dist_nm * 1852 / num_legs

airfield_elevation = 750
cruise_altitude = 3000
loiter_altitude = 1200
temp_offset = 18
climb_rate = 2.5
cruise_speed = 60
loiter_duration = 75*60 # [s]

macf = mission_performance.MAircraft(airfield_elevation)
macf.dT = temp_offset

if better_batteries:
    macf.acf.max_bat_cap = macf.acf.max_bat_cap * better_battery_factor
    macf.bat_cap = macf.acf.max_bat_cap
    macf.b_lst[0] = macf.bat_cap
    macf.last_bat_cap = macf.bat_cap

#
# flight
#

for i in range(num_legs):
    # takeoff
    print(f"**************** Takeoff distance: {macf.acf.takeoff_ground_run(macf.W, 750, 18, 0, 'grass')[0]:.2f} m *******************")
    macf.take_off(airfield_elevation, "grass", electric=i==num_legs-2)
    
    # acceleration to ideal CL for climbing
    V_start_climb = macf.acf.V_Prmin(macf.W, macf.h, macf.dT)
    macf.fly_accelerate_const_alt(V_start_climb, electric=i==num_legs-2)
    
    # climb
    macf.fly_const_CL_const_climbrate(climb_rate, cruise_altitude, electric=i==num_legs-2)
    
    # accelerate to cruise speed
    macf.fly_accelerate_const_alt(cruise_speed, electric=i==num_legs-2)
    
    # cruise
    macf.fly_const_V_const_alt(leg_dist_m, electric=i>=num_legs-2) # if we have batcap left over this will use it also during fourth cruise
    
    if i < num_legs-1:
        # descent
        macf.fly_const_V_descent(airfield_elevation) # TODO: const climbrate descent
    
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
        macf.fly_const_V_descent(airfield_elevation) # TODO: const climbrate descent
        
        # land
        macf.land()


# plot and print

macf.plot_flight()
macf.print_energy_level()

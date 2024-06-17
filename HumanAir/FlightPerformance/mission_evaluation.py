# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 10:41:53 2024

@author: Alex
"""

import mission_performance

tot_dist_nm = 600 # 1400
leg_dist_m = tot_dist_nm * 1852 / 4

macf = mission_performance.MAircraft(750)
# macf.V = 30
# #macf._fly_const_CL_const_climbrate(2.5, 3000)
# macf._fly_const_V_const_alt(1000000)
# print(macf.ground_distance)
# print(macf.flight_time)
# print(macf.V)
# #macf.take_off(750, "grass", electric=False)
# #macf.climb(3000, 2.5)

macf.take_off(750, "grass", electric=True)
V_start_climb = macf.acf.V_Prmin(macf.W, macf.h, macf.dT)
macf.fly_accelerate_const_alt(V_start_climb, electric=True)
macf.fly_const_CL_const_climbrate(2.5, 3000, electric=True)
macf.fly_accelerate_const_alt(60, electric=True)
macf.fly_const_V_const_alt(leg_dist_m, electric=True)
macf.fly_const_V_descent(750)
macf.land()
for _ in range(2):
    macf.take_off(750, "grass", electric=False)
    V_start_climb = macf.acf.V_Prmin(macf.W, macf.h, macf.dT)
    macf.fly_accelerate_const_alt(V_start_climb, electric=False)
    macf.fly_const_CL_const_climbrate(2.5, 3000, electric=False)
    macf.fly_accelerate_const_alt(60, electric=False)
    macf.fly_const_V_const_alt(leg_dist_m, electric=False)
    macf.fly_const_V_descent(750)
    macf.land()
macf.take_off(750, "grass", electric=False)
V_start_climb = macf.acf.V_Prmin(macf.W, macf.h, macf.dT)
macf.fly_accelerate_const_alt(V_start_climb, electric=False)
macf.fly_const_CL_const_climbrate(2.5, 3000, electric=False)
macf.fly_accelerate_const_alt(60, electric=False)
macf.fly_const_V_const_alt(leg_dist_m, electric=False)
macf.fly_const_V_descent(1200)
V_loiter = macf.acf.V_Prmin(macf.W, macf.h, macf.dT)
macf.fly_accelerate_const_alt(V_loiter)
macf.fly_const_CL_const_alt(115*60)
macf.fly_const_V_descent(750)
macf.plot_flight()
macf.print_energy_level()

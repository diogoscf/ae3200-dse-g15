# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 10:41:53 2024

@author: Alex
"""

import mission_performance

macf = mission_performance.MAircraft()
macf.V = 30
#macf._fly_const_CL_const_climbrate(2.5, 3000)
macf._fly_const_V_const_alt(1000000)
print(macf.ground_distance)
print(macf.flight_time)
print(macf.V)
#macf.take_off(750, "grass", electric=False)
#macf.climb(3000, 2.5)
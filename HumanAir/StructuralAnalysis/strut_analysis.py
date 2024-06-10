"""Strut Analysis """ 

import numpy as np
from optimisation import get_deflection

#Strut Force - P 
def strut_force(MOI,y, M, A_strut0, E,l_strut, theta_strut, max_iter= 100, tol=1e-6):
    A = A_strut0
    P = 0
    for i in range(max_iter):
        #deflection at 40% of halfspan
        v = get_deflection(MOI,y,M,E)
        v_40 = v[round(0.4*len(v))]

        #deflection at an angle 
        delta = v_40/np.cos(theta_strut)

        #strut force
        P_new = delta*A*E/l_strut

        if abs(P_new - P) < tol:
            break

        #update
        P = P_new
        A = P*l_strut/(delta*E)

        V_strut = P
        Vz_strut = P*np.cos(theta_strut)
        Vy_strut = P*np.sin(theta_strut)

    return Vz_strut, Vy_strut, V_strut, A


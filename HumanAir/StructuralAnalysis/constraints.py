import numpy as np 
from LoadDistributions import Mx, Vy
E = 68e9 #youngs modulus aluminium 
I_stringer= 20.2e-12 #inertia of stringer (m^4)
A_stringer = 0.005 #area of stringer (m^2)
L_halfspan = 
L_eff = 0.092*L_halfspan
sigma_yield = 40e6 #[Pa]

def stringer_buckling_constraint(variables, Mx, A_stringer):

    P_cr = (np.pi)**2*E*I_stringer/L_eff**2 #critical force at which stringer will buckle
    
    #inertia
    h_max = h_15c/2
    t_spar = variables[:n_halfspan].flatten()
    t_skin = variables[n_halfspan : 2 * n_halfspan].flatten()
    no_stringers = variables[2 * n_halfspan :].flatten()
    I = get_I(t_spar, t_skin, no_stringers, A_stringer)

    #compressive stress
    sigma_bending = Mx*h_max/I
    sigma_axial = Vy/A_stringer
    sigma = sigma_bending+sigma_axial 

    P = sigma/(no_stringers*A_stringer)
    return P_cr - np.max(P)

def tensile_failure_constraint(Variables, Mc, A_stringer):
    #inertia
    h_max = h_15c/2
    t_spar = variables[:n_halfspan].flatten()
    t_skin = variables[n_halfspan : 2 * n_halfspan].flatten()
    no_stringers = variables[2 * n_halfspan :].flatten()
    I = get_I(t_spar, t_skin, no_stringers, A_stringer)

    #bending stress
    sigma = Mx*h_max/I
    return sigma_yield - np.max(sigma)

def twist_angle_constraint():
    return max_twist - twist

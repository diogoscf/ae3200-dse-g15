import numpy as np 
from LoadDistributions import Mx, Vy, Vy_strut
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

#shear buckling of panel 
def get_Q(t_spar):
    Q = h_15c**2/8*t_spar + h_50c**2/8*t_spar
    return Q

def shear_buckling_constraint(variables, V, nu, ks = 9.5):
    #inertia
    h_max = h_15c/2
    t_spar = variables[:n_halfspan].flatten()
    t_skin = variables[n_halfspan : 2 * n_halfspan].flatten()
    no_stringers = variables[2 * n_halfspan :].flatten()
    I = get_I(t_spar, t_skin, no_stringers, A_stringer)

    #critical shear stress
    b = h_50c #critical height between the two spars is shorter spar 
    shear_critical = (np.pi)**2*ks*E/(12*(1-nu)**2)*(t_spar/b)**2

    #first moment of area
    Q = get_Q(t_spar)

    #shear
    shear = V*Q/(I*t_spar)

    return shear_critical - np.max(shear)

#Buckling of a stiffened panel 
def get_A(t_skin, t_spar, no_stringers, A_stringer, spar_pos):
    A = h_15c*t_spar + h_50c*t_spar + no_stringers*A_stringer + 2*t_skin*(spar_pos[1]-spar_pos[0])
    return A

def stiffened_skin_buckling_constraint(variables, spar_pos, A_stringer E, nu, C_skin=4, C_stiffener = 0.425, t_stiffener = 2, b_stiffener = 25, C_we = 4):  
    '''Critical Stress '''

    #parameters
    t_skin = variables[n_halfspan : 2 * n_halfspan].flatten()
    no_stringers = variables[2 * n_halfspan :].flatten()
    b = 2*(spar_pos[1]-spar_pos[0])/no_stringers

    #crippling stress of the stiffener 
    alpha = 0.8 #empirical values
    n = 0.6 #empirical values
    stress_crip_ratio = alpha*(C_stiffener/sigma_yield*np.pi**2*E/(12*(1-nu**2))*(t_stiffener/b_stiffener)**2)**(1-n)
    
    if stress_crip_ratio < 1:
        stress_crippling = sigma_yield*stress_crip_ratio #check it cripples before it yields
    
    #calculate effective width of the stiffener (accounting for attachement between the stringer and the skin)
    we = t_skin/2*np.sqrt(C_we*np.pi**2/(12*(1-nu**2)))*np.sqrt(E/stress_crippling)

    #Calculate the initial skin buckling taking into account the 2we
    stress_critical_skin = C_skin*np.pi**2*E/(12*(1-nu**2))*(t_skin/(b-2*we))**2

    #final critical stress of the entire panel 
    stress_critical = (stress_critical_skin*t_skin*(b-2*we)+stress_crippling*(A_stringer+2*we*t_skin))/(t_skin*(b-2*we)+(A_stringer+2*we*t_skin))

    """Actual Stress"""
    #inertia 
    h_max = h_15c/2 # --> CHANGE THIS
    t_spar = variables[:n_halfspan].flatten()
    I = get_I(t_spar, t_skin, no_stringers, A_stringer)

    #bending stress
    stress_bending = Mx*h_max/I

    #area
    A = get_A(t_skin,t_spar,no_stringers,A_stringer,spar_pos)

    #axial stress 
    stress_axial = Vy_strut/A

    #final stress
    stress = np.abs(stress_axial)+stress_bending

    return stress_critical - np.max(stress)

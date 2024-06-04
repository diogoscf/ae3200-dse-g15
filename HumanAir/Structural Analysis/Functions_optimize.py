import numpy as np
from scipy.optimize import minimize
from scipy.integrate import cumtrapz
from LoadDistributions import Mx

"""
    Minimize the weight of the wing while optimizing for these three variables:
    - spar thickness (varying over the wing)
    - no. of stringers (equally spaced) (varies three times along the span)
    - skin thickness
    """
#Parameters
rho= 2710 #density of aluminium (kg/m^3)
n = 10 #number of segments, discretization
L = 17.61/2 #length of the halfspan (m)
h_15c = 0 #height of the spar at 15% of the chord, as a function of y -- NOTE: MAKE THEM ARRAYS like np.full of t_spar0
h_50c = 0 #height of the spar at 50% of the chord, as a function of y 
htot = h_15c + h_50c #total height, just for calculation ease, as a function of y
h_avemax = (htot/4) #averaged "max" height from the central line, as a function of y
w_top = 0 #width of top "straight" skin, as a function of y 
w_bottom = 0 #width of bottom "straight" skin, as a function of y
wtot = w_bottom+w_top #total width, just for calculation ease
A_stringer = 0 #area of stringer for 20mm L stringer
E = 68e9 #Youngs Modulus for aluminium (Pa)
sigma_yield = 40e6 #Yield strength for aluminium (Pa)


#Length of each segment 
delta_y = L/n

#Initial guess for thickness of spar
t_spar0 = np.full(n,0.01) #array with thickness per segment [m]
t_skin0 = 0.07 #constant thickness [m]
no_stringers0 = np.array(50,30,10) #just initial values 
t0 = 

#Objective function: minimize weight
def objective(t_spar, t_skin, no_stringers): 
    W_spar = np.sum(rho*t_spar*htot* delta_y)
    W_skin = np.sum(rho*t_skin**wtot*delta_y) 
    W_stringers = rho*A_stringer*(no_stringers[0]*0.4*L +no_stringers[1]*0.3*L+no_stringers[2]*0.3*L)
    weight = W_skin+W_spar+W_stringers
    return weight 

def get_deflection(I):
    y = np.linspace(0,L,n)
    M = Mx #get M(y) of most critical case from jennifer 
    integrand = M/I
    theta = -1/E*cumtrapz(integrand,y,initial =0)
    v= -1/E*cumtrapz(theta, y, initial= 0)
    return np.max(v)

def get_I(t_spar,t_skin,no_stringers):
    I_skin = t_skin*(w_top+w_bottom)*h_avemax**2
    I_spar = 1/12*t_spar*(h_15c**3+h_50c**3)
    I_stringers =A_stringer*no_stringers*h_avemax**2
    I = I_skin+I_spar+I_stringers
    return I

#Stiffness constraint function
def deflection_constraint(t_spar,t_skin,no_stringers):
    #maximum deflection 
    max_deflection = 1.7 #[m] - estimate
    #deflection 
    I = get_I(t_spar,t_skin,no_stringers)
    deflection = get_deflection(I)
    return max_deflection - deflection

def bending_stress_constraint(t_spar,t_skin, no_stringers):
    stresses = []
    for i in range(n):
        h_avemax = h_avemax[i]
        M = Mx[i]
        I = get_I(t_spar,t_skin,no_stringers)[i]
        stress = M*h_avemax/I
        stresses.append(sigma_yield - stress)
    return np.min(stresses)

#def get_A(t_spar,t_skin,no_stringers):
    A_skin = t_skin*(w_top+w_bottom)
    A_spar = t_spar*(htot)
    A_stringers = no_stringers*A_stringer
    A = A_skin+A_spar+A_stringers
    return A

#constraint on the axial stress (compressive) on the first segment carrying the strut
#def axial_stress(t_spar, t_skin, no_stringers):
    stresses = []
    for i in range(n):
        A = get_A(t_spar,t_skin, no_stringers)[i]
        force = Vy[i]
        stress = force/A
        stresses.append(sigma_yield - stress)
        return np.min(stresses)

#buckling constraint - due to the strut compressive stress
#def column_buckling_constraint(t_spar,t_skin,no_stringers):
    stresses = []
    for i in range(n):
        
    return np.min(stresses)

#Bounds for thickness (must be greater than zero)
bounds = [(0.001, 0.1)] * n

# Constraints dictionary
con1 = {'type': 'ineq', 'fun': bending_stress_constraint}
constraints = [con1]

# Perform the optimization
solution = minimize(objective, t0, method='SLSQP', bounds=bounds, constraints=constraints)

# Optimized thickness values
optimized_thickness = solution.x

# Print the results
print('Optimized thickness for each segment:', optimized_thickness)
print('Minimum weight:', solution.fun)
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import cumtrapz
from LoadDistributions import Mx
from TorsionalStiffness import TorsionalStiffness

"""
    Minimize the weight of the wing while optimizing for these three variables:
    - spar thickness (varying over the wing)
    - no. of stringers (equally spaced) (varies three times along the span)
    - skin thickness
    """

def calculate_sizes(percentage, total_size):
    sizes = np.array([int(total_size * p) for p in percentage])

    # Adjust for any rounding errors to ensure the total size is correct
    residual = total_size - np.sum(sizes)
    if residual != 0:
        # Distribute the residual to the middle elements
        middle = len(sizes) // 2
        sizes[middle] += residual // 2
        if len(sizes) % 2 == 0:
            sizes[middle - 1] += residual // 2
        else:
            sizes[middle + 1] += residual // 2
    return sizes

def create_segments(sizes, no_string):
    segments = [np.full(size, value) for size, value in zip(sizes, no_string)]
    seg = np.concatenate(segments)
    return seg.reshape((len(seg), 1))

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
'''
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
'''

########## input ###########
Sw = 34.56  # [m2]
taper_ratio = 0.4
AoA = -6 # [deg]
t1_spar = 0.010 # [m] thickness at the tip
t2_spar = 0.025 # [m] thickness at the root
t_skin = 0.007 # [m] thickness of skin
file_path = 'HumanAir/WingBox/airfoil.txt'
file_path_y = 'HumanAir\WingBox\Cl_DATA.txt'
A_str = 0.02
Cr = 2.5 # [m] root chord length
b = 19.93
x_pos = np.array([0.15, 0.5])
n = 100

# Initialise Class
torisonal_stiffness = TorsionalStiffness(file_path, file_path_y, Sw, taper_ratio, AoA, n, t1_spar, t2_spar, t_skin, x_pos, A_str, Cr,b)
h_mid, h_s1s2 = torisonal_stiffness.h_s1s2()
l_box_up, l_box_down = torisonal_stiffness.d_s1s2()

#Parameters
rho= 2710 #density of aluminium (kg/m^3)
n = 10 #number of segments, discretization
L = 17.61/2 #length of the halfspan (m)
h_15c = h_s1s2[0] #height of the spar at 15% of the chord, as a function of y -- NOTE: MAKE THEM ARRAYS like np.full of t_spar0
h_50c = h_s1s2[1] #height of the spar at 50% of the chord, as a function of y 
htot = h_15c + h_50c #total height, just for calculation ease, as a function of y
h_avemax = (htot/4) #averaged "max" height from the central line, as a function of y
w_top = l_box_up #width of top "straight" skin, as a function of y 
w_bottom = l_box_down #width of bottom "straight" skin, as a function of y
wtot = w_bottom+w_top #total width, just for calculation ease
A_stringer = 0 #area of stringer for 20mm L stringer
E = 68e9 #Youngs Modulus for aluminium (Pa)
sigma_yield = 40e6 #Yield strength for aluminium (Pa)


#Length of each segment 
delta_y = L/n


#Objective function: minimize weight
percentage = [0.1, 0.15, 0.25, 0.25, 0.15, 0.1]
no_string = [20, 30, 50, 50, 30, 20]

#Initial guess for thickness of spar
t_spar0 = torisonal_stiffness.t_spar() #array with thickness per segment [m]
t_skin0 = 0.07 * np.ones((len(t_spar0),1)) #constant thickness [m]
size = calculate_sizes(percentage, len(t_spar0) )
no_stringers0 = create_segments(size, no_string)
#print('t_spar0', t_spar0)
#print('t_skin0', t_skin0)
#print('no_stringers',  no_stringers0)

t0 = np.hstack((t_spar0, t_skin0, no_stringers0))

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
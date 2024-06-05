import numpy as np
from scipy.optimize import minimize
from scipy.integrate import cumtrapz
from LoadDistributions import Mx
from TorsionalStiffness import TorsionalStiffness
from Functions import import_data2

<<<<<<< HEAD
"""
    Minimize the weight of the wing while optimizing for these three variables:
    - spar thickness (varying over the wing)
    - no. of stringers (equally spaced) (varies three times along the span)
    - skin thickness
    """
#Parameters
rho= 2710 #density of aluminium (kg/m^3)
n = 38 #number of segments, discretization
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
=======
def calculate_sizes(percentage, total_size):
    sizes = np.array([int(total_size * p) for p in percentage])
    residual = total_size - np.sum(sizes)
    if residual != 0:
        middle = len(sizes) // 2
        sizes[middle] += residual // 2
        if len(sizes) % 2 == 0:
            sizes[middle - 1] += residual // 2
        else:
            sizes[middle + 1] += residual // 2
    return sizes
>>>>>>> 634daee2d4b2f8dfcf011aeda80f25201da64f89

def create_segments(sizes, no_string):
    segments = [np.full(size, value) for size, value in zip(sizes, no_string)]
    seg = np.concatenate(segments)
    return seg.reshape((len(seg), 1))

def objective(vars, htot, wtot, rho, b, n, A_stringer): 
    t_spar = vars[:n].reshape((n, 1))
    t_skin = vars[n:2*n].reshape((n, 1))
    no_stringers = vars[2*n:].reshape((n, 1))
    
    L = b / n
    W_spar = np.sum(rho * t_spar * htot * L)
    W_skin = np.sum(rho * t_skin * wtot * L) 
    W_stringers = np.sum(rho * A_stringer * no_stringers * L)
    weight = W_skin + W_spar + W_stringers
    return weight

def get_deflection(I, y, M):
    integrand = M / I
    theta = -1 / E * cumtrapz(integrand, y, initial=0)
    v = -1 / E * cumtrapz(theta, y, initial=0)
    return np.max(v)

def get_I(t_spar, t_skin, no_stringers, A_stringer):
    I_skin = t_skin * (w_top + w_bottom) * h_avemax**2
    I_spar = 1/12 * t_spar * (h_15c**3 + h_50c**3)
    I_stringers = A_stringer * no_stringers * h_avemax**2
    I = I_skin + I_spar + I_stringers
    return I

def deflection_constraint(vars, y, Mx, A_stringer, max_deflect):
    t_spar = vars[:n].reshape((n, 1))
    t_skin = vars[n:2*n].reshape((n, 1))
    no_stringers = vars[2*n:].reshape((n, 1))
    I = get_I(t_spar, t_skin, no_stringers, A_stringer)
    deflection = get_deflection(I, y, Mx)
    return max_deflect - deflection

<<<<<<< HEAD
#def bending_stress_constraint(t_spar,t_skin, no_stringers):
    stresses = []
    for i in range(n):
        h_avemax = h_avemax[i]
        M = Mx[i]
        I = get_I(t_spar,t_skin,no_stringers)[i]
        stress = M*h_avemax/I
        stresses.append(sigma_yield - stress)
    return np.min(stresses)
=======
########## Input ###########
Sw = 34.56  # [m^2]
taper_ratio = 0.4
AoA = -6  # [deg]
t1_spar = 0.010  # [m] thickness at the tip
t2_spar = 0.025  # [m] thickness at the root
t_skin = 0.007  # [m] thickness of skin
file_path = 'HumanAir/WingBox/airfoil.txt'
file_path_y = 'HumanAir/WingBox/Cl_DATA.txt'
A_str = 8e-5
Cr = 2.5  # [m] root chord length
b = 19.93
x_pos = np.array([0.15, 0.5])
Cl_DATA = import_data2('HumanAir/Structural Analysis/Cl_DATA.txt')
n = len(Cl_DATA[AoA]['coefficient'])
max_deflect = 2 # [m]
>>>>>>> 634daee2d4b2f8dfcf011aeda80f25201da64f89

# Initialize Torsional Stiffness Class
torsional_stiffness = TorsionalStiffness(file_path, file_path_y, Sw, taper_ratio, AoA, n, t1_spar, t2_spar, t_skin, x_pos, A_str, Cr, b)
h_mid, h_s1s2 = torsional_stiffness.h_s1s2()
l_box_up, l_box_down = torsional_stiffness.d_s1s2()
c, y = torsional_stiffness.chord()

# Parameters
rho = 2710  # density of aluminium (kg/m^3)
h_15c = h_s1s2[:, 0]  # height of the spar at 15% of the chord, as a function of y
h_50c = h_s1s2[:, 1]  # height of the spar at 50% of the chord, as a function of y
htot = (h_15c + h_50c).reshape((n, 1))  # total height, just for calculation ease, as a function of y
h_avemax = htot / 4  # averaged "max" height from the central line, as a function of y
w_top = l_box_up.reshape((n, 1))  # width of top "straight" skin, as a function of y 
w_bottom = l_box_down.reshape((n, 1))  # width of bottom "straight" skin, as a function of y
wtot = w_bottom + w_top  # total width, just for calculation ease
E = 68e9  # Young's Modulus for aluminium (Pa)
sigma_yield = 40e6  # Yield strength for aluminium (Pa)

percentage = [0.1, 0.15, 0.25, 0.25, 0.15, 0.1]
no_string = [20, 30, 50, 50, 30, 20]

# Initial guess for thickness of spar
t_spar0 = torsional_stiffness.t_spar()
# Constant skin thickness [m]
t_skin0 = t_skin * np.ones((len(t_spar0), 1)) 
size = calculate_sizes(percentage, len(t_spar0))
no_stringers0 = create_segments(size, no_string)

# Initial guess
t0 = np.hstack((t_spar0, t_skin0, no_stringers0))

# Define the bounds for each variable
bounds_t_spar = [(0.005, 0.1)] * n
bounds_t_skin = [(0.005, 0.1)] * n
bounds_no_stringers = [(1, 100)] * n
bounds = bounds_t_spar + bounds_t_skin + bounds_no_stringers

# Constraints dictionary
<<<<<<< HEAD
con1 = {'type': 'ineq', 'fun':deflection_constraint}
constraints = [con1]
=======
constraints = [{'type': 'ineq', 'fun': deflection_constraint, 'args': (y, Mx, A_str, max_deflect)}]
>>>>>>> 634daee2d4b2f8dfcf011aeda80f25201da64f89

# Perform the optimization
solution = minimize(objective, t0, args=(htot, wtot, rho, b, n, A_str), method='SLSQP', bounds=bounds, constraints=constraints)

# Optimized thickness values
optimized_thickness = solution.x
optimized_t_spar = optimized_thickness[:n]
optimized_t_skin = optimized_thickness[n:2*n]
optimized_no_stringers = optimized_thickness[2*n:]

# Print the results
print('Optimized thickness for spar:', optimized_t_spar)
print('Optimized thickness for skin:', optimized_t_skin)
print('Optimized number of stringers:', optimized_no_stringers)
print('Minimum weight:', solution.fun)

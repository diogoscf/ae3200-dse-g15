import numpy as np
from scipy.optimize import minimize
from scipy.integrate import cumulative_trapezoid
import sys
import os
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.StructuralAnalysis.WingStructure import WingStructure
from Functions import import_data2
from LoadDistributions import Mx


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


def create_segments(sizes, no_string):
    segments = [np.full(size, value) for size, value in zip(sizes, no_string)]
    seg = np.concatenate(segments)
    return seg.reshape((len(seg),))


def objective(vars, htot, wtot, rho, b, n, A_stringer):
    t_spar = vars[:n].flatten()
    t_skin = vars[n : 2 * n].flatten()
    no_stringers = vars[2 * n :].flatten()

    L = b / n
    W_spar = np.sum(rho * t_spar * htot * L)
    W_skin = np.sum(rho * t_skin * wtot * L)
    W_stringers = np.sum(rho * A_stringer * no_stringers * L)
    weight = W_skin + W_spar + W_stringers
    return weight


def get_deflection(I, y, M):
    integrand = M / I
    theta = -1 / E * cumulative_trapezoid(integrand, y, initial=0)
    v = -1 / E * cumulative_trapezoid(theta, y, initial=0)
    return np.max(v)


def get_I(t_spar, t_skin, no_stringers, A_stringer):
    I_skin = t_skin * (w_top + w_bottom) * h_avemax**2
    I_spar = 1 / 12 * t_spar * (h_15c**3 + h_50c**3)
    I_stringers = A_stringer * no_stringers * h_avemax**2
    I = I_skin + I_spar + I_stringers
    return I


def deflection_constraint(vars, y, Mx, A_stringer, max_deflect):
    t_spar = vars[:n_halfspan].flatten()
    t_skin = vars[n_halfspan : 2 * n_halfspan].flatten()
    no_stringers = vars[2 * n_halfspan :].flatten()
    I = get_I(t_spar, t_skin, no_stringers, A_stringer)
    deflection = get_deflection(I, y, Mx)
    return max_deflect - deflection


########## Input ###########
Sw = 34.56  # [m^2]
taper_ratio = 0.4
AoA = -6  # [deg]
t1_spar = 0.010  # [m] thickness at the tip
t2_spar = 0.025  # [m] thickness at the root
t_skin = 0.007  # [m] thickness of skin
file_path = "HumanAir/WingBox/airfoil.txt"
file_path_y = "HumanAir/WingBox/Cl_DATA.txt"
A_str = 8e-5
Cr = 2.5  # [m] root chord length
b = 19.93
x_pos = np.array([0.15, 0.5])
Cl_DATA = import_data2("HumanAir/StructuralAnalysis/Cl_DATA.txt")
n_fullspan = len(Cl_DATA[AoA]["coefficient"])
n_halfspan = n_fullspan // 2
max_deflect = 2  # [m]

# Initialize Torsional Stiffness Class
torsional_stiffness = WingStructure(
    file_path, file_path_y, Sw, taper_ratio, AoA, n_fullspan, t1_spar, t2_spar, t_skin, x_pos, A_str, Cr, b
)
h_mid, h_s1s2 = torsional_stiffness.h_s1s2()
l_box_up, l_box_down = torsional_stiffness.d_s1s2()
c, y = torsional_stiffness.chord_distribution()
h_mid, h_s1s2 = h_mid[h_mid.shape[0] // 2 :], h_s1s2[h_s1s2.shape[0] // 2 :]
l_box_up, l_box_down = l_box_up[l_box_up.shape[0] // 2 :], l_box_down[l_box_down.shape[0] // 2 :]
c, y = c[c.shape[0] // 2 :], y[y.shape[0] // 2 :]

# Parameters
rho = 2710  # density of aluminium (kg/m^3)
h_15c = h_s1s2[:, 0]  # height of the spar at 15% of the chord, as a function of y
h_50c = h_s1s2[:, 1]  # height of the spar at 50% of the chord, as a function of y
htot = (h_15c + h_50c).flatten()  # total height, just for calculation ease, as a function of y
h_avemax = htot / 4  # averaged "max" height from the central line, as a function of y
w_top = l_box_up.flatten()  # width of top "straight" skin, as a function of y
w_bottom = l_box_down.flatten()  # width of bottom "straight" skin, as a function of y
wtot = w_bottom + w_top  # total width, just for calculation ease
E = 68e9  # Young's Modulus for aluminium (Pa)
sigma_yield = 40e6  # Yield strength for aluminium (Pa)

percentage = [0.5, 0.3, 0.2]
no_string = [50, 30, 20]

# Initial guess for thickness of spar
t_spar0 = torsional_stiffness.t_spar().flatten()
t_spar0 = t_spar0[t_spar0.shape[0] // 2 :]
# Constant skin thickness [m]
t_skin0 = t_skin * np.ones((len(t_spar0)))
size = calculate_sizes(percentage, len(t_spar0))
no_stringers0 = create_segments(size, no_string)
# Initial guess
t0 = np.hstack((t_spar0, t_skin0, no_stringers0)).flatten()

# Define the bounds for each variable
bounds_t_spar = [(0.005, 0.1)] * n_halfspan
bounds_t_skin = [(0.005, 0.1)] * n_halfspan
bounds_no_stringers = [(1, 100)] * n_halfspan
bounds = bounds_t_spar + bounds_t_skin + bounds_no_stringers

# Constraints dictionary
constraints = [{"type": "ineq", "fun": deflection_constraint, "args": (y, Mx, A_str, max_deflect)}]

# Perform the optimization
solution = minimize(
    objective, t0, args=(htot, wtot, rho, b, n_halfspan, A_str), method="SLSQP", bounds=bounds, constraints=constraints
)

# Optimized thickness values
optimized_thickness = solution.x
optimized_t_spar = optimized_thickness[:n_halfspan]
optimized_t_skin = optimized_thickness[n_halfspan : 2 * n_halfspan]
optimized_no_stringers = optimized_thickness[2 * n_halfspan :]

# Print the results
print("Optimized thickness for spar:", optimized_t_spar)
print("Optimized thickness for skin:", optimized_t_skin)
print("Optimized number of stringers:", optimized_no_stringers)
print("Minimum weight:", solution.fun)

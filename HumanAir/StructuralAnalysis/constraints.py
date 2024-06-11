import numpy as np
from optimisation import constraint_data_from_variables, MOI_from_variables, get_axial_stresses, get_axial_forces

def stringer_buckling_constraint(variables, Mx, Vy, A_stringer, hmax, E, I_stringer, L_eff, MOI_args):
    P_cr = (np.pi) ** 2 * E * I_stringer / L_eff**2  # critical force at which stringer will buckle
    MOI, stringers = constraint_data_from_variables(variables, *MOI_args)[:2]

    # compressive stress
    P = get_axial_forces(Mx, Vy, MOI, hmax, stringers, A_stringer)
    # print(P_cr - np.max(P), flush=True)
    return (P_cr - np.max(P)) / P_cr
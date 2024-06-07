import numpy as np
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from HumanAir.aircraft_data import aircraft_data
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from scipy.linalg import eig





def Divergence(K_theta, CL_alpha, a, b, rho):
    """
    Calculate the divergence speed of an aircraft based on the given parameters. 
    The typical section is the airfoil cross-section of the wing at 70% of the wing span.

    Parameters
    ----------
    K_theta : float
        The torsional stiffness of the full 3D wing from root to typical section.
    CL_alpha : float
        The lift curve slope of the full 3D wing wing.
    a : float
        Distance from half-chord to the elastic axis of the typical section airfoil. non-dimensionalized by the half-chord length.
    b : float
        Half-chord length of the typical section (not wingspan!!!)

    Returns
    -------
    V_div : float
        The divergence speed of the aircraft.
    """
    q = K_theta / ((2 * b) * CL_alpha * (1/2 + a) * b) # where S = 2 * b since we analyse per unit span.
    V_div = np.sqrt(2 * q / rho)
    return V_div


def Reversal(c, K_theta, CL_alpha, CM_ac_beta, b, rho):
    """
    Calculate the flutter speed of an aircraft based on the given parameters.

    Parameters
    ----------
    c : float
        The distance from the half-chord  to the hinge of the aileron        
    K_theta : float
        The torsional stiffness of the full 3D wing from root to typical section.
    CL_alpha : float
        The lift curve slope of the full 3D wing wing.
    CM_ac_beta : float
        The change in pitching moment about the typical section aerodynamic section with deflection of the control surface (aileron).
    b : float
        Half-chord length of the typical section (not wingspan!!!)

    Returns
    -------
    V_rev : float
        The aileron reversal speed of the aircraft.
    """
    CL_beta = 2 * (np.sqrt(1 - c**2) + np.arccos(c)) # The lift variation with deflection of the control surface (aileron).
    q = - CL_beta * K_theta / (CL_alpha * CM_ac_beta * 2 * b *(2 * b))
    V_rev = np.sqrt(2 * q / rho)
    return V_rev


def rho_altitude(rho):
    """
    calculate the altitude form a given air density.

    Parameters
    ----------
    rho : float
        The air density.

    Returns
    -------
    altitude : float
        The altitude corresponding to the given air density.
    """
    altitude = (- np.log(rho / 1.225) * (287.05 * 288.15) / 9.80665) / (1 + np.log(rho / 1.225) * (287.05 * (-0.0065)/ 9.80665))
    return altitude




def Flutter(K_h, K_theta, rho_arr, V_arr, CL_alpha, b, a, x_theta, m_arr, I_theta):
    """
    Calculate the flutter speed of an aircraft based on the given parameters.

    Parameters
    ----------
    K_h : float
        The torsional stiffness of the full 3D wing from root to typical section.
    K_theta : float
        The torsional stiffness of the full 3D wing from root to typical section.
    rho_arr : numpy array
        Array containing all air densities you would like to analyse aircraft flutter for.
    V_arr : numpy array
        Array containing all aircraft velocities from 0 to dive speed * 1.15 (+15% for certification, according to Dr. Roeland De Breuker).
    CL_alpha : float
        The lift curve slope of the full 3D wing wing.
    b : float
        Half-chord length of the typical section (not wingspan!!!)
    a : float
        Distance from half-chord to the elastic axis of the typical section airfoil. non-dimensionalized by the half-chord length.
    x_theta: float
        Distance from the elastic axis to the COG of the main airfoil component (not the aileron).
    m_arr : numpy array
        Array containing all aircraft mass configurations you would like to analyse aircraft flutter for (fully loaded, empty, etc.)
    I_theta : float
        The torsional moment of inertia of the wing about the elastic axis.
        
    Returns
    -------
    V_flut : float
        The flutter speed of the aircraft.
    """
    results = []
    for rho in rho_arr:
        for m in m_arr:
            eigenvalues = []
            for V in V_arr:
                # Define matrices
                CL_alpha_dot =  CL_alpha * (1-a) # The rate of change of the lift for a rate of change in angle of attack.
                q = 1/2 * rho * V**2
                S_theta = m * x_theta * b
                K_s = np.array([[K_h, 0], [0, K_theta]])
                K_a = np.array([[0, -q * (2 * b) * CL_alpha], [0, -q * (2 * b) * CL_alpha * (1/2 + a)]])
                C_a = np.array([[q * (2 * b) * CL_alpha * 1/V, q * (2 * b) * CL_alpha_dot * b/V], [q * (2 * b) * CL_alpha * 1/V (1/2 + a) * b, q * (2 * b) * b/V * b * CL_alpha * a * (1/2 - a)]])
                M_s = np.array([[m, S_theta], [S_theta, I_theta]])

                # Define the function representing the equation M_s * p^2 - C_a * p + (K_s - K_a) = 0
                def equation(p):
                    term1 = M_s * p**2  # M_s * p^2
                    term2 = C_a * p  # C_a * p
                    term3 = K_s - K_a  # K_s - K_a
                    result = term1 - term2 + term3
                    return result.flatten()  # Return a flattened array for fsolve

                # Initial guess for the scalar p
                p0 = 1

                # Solve the equation
                solution = fsolve(equation, p0)
                eigenvalues.append(solution)
            
            # Cycle through all eigenvalue vectors, and find the index of the vector when the real value turns positive, this happens at the flutter speed for the current configuration.
            for i in range(len(eigenvalues)):
                if eigenvalues[i][0] > 0:
                    V_flut = V_arr[i]
                    flutter = True
                    break
                else:
                    flutter = False
            
            if flutter: # Save the configuration and flutter speed
                result = {
                    "V_max_flight [m/s]": V_arr[-1]*1.15, 
                    "altitude [m]": rho_altitude(rho),
                    "mass_configuration [kg]": m,
                    "V_flutter [m/s]": V_flut, 
                    "a": a,
                    "b": b,
                    "x_theta": x_theta,
                }
                results.append(result)

    if len(results) == 0:
        return "No flutter detected for the given configurations."
    else:
        return results
                

def flutter_diagram(K_h, K_theta, rho, V_arr, CL_alpha, b, a, x_theta, m, I_theta):
    """
    Calculate the flutter speed of an aircraft based on the given parameters.

    Parameters
    ----------
    K_h : float
        The torsional stiffness of the full 3D wing from root to typical section.
    K_theta : float
        The torsional stiffness of the full 3D wing from root to typical section.
    rho : float
        The aircraft altitude you would like to analyse aircraft flutter for.
    V_arr : numpy array
        Array containing all aircraft velocities from 0 to dive speed * 1.15 (+15% for certification, according to Dr. Roeland De Breuker).
    CL_alpha : float
        The lift curve slope of the full 3D wing wing.
    b : float
        Half-chord length of the typical section (not wingspan!!!)
    a : float
        Distance from half-chord to the elastic axis of the typical section airfoil. non-dimensionalized by the half-chord length.
    x_theta: float
        Distance from the elastic axis to the COG of the main airfoil component (not the aileron).
    m : float
        The aircraft mass configuration you would like to analyse aircraft flutter for (fully loaded, empty, etc.)
    I_theta : float
        The torsional moment of inertia of the wing about the elastic axis.
        
    Returns
    -------
    V_flut : float
        The flutter speed of the aircraft.
    """
    eigenvalues_lst = []
    for V in V_arr:
        # Define matrices
        CL_alpha_dot =  CL_alpha * (1-a) # The rate of change of the lift for a rate of change in angle of attack.
        q = 1/2 * rho * V**2
        S_theta = m * x_theta * b
        K_s = np.array([[K_h, 0], [0, K_theta]])
        K_a = np.array([[0, -q * (2 * b) * CL_alpha], [0, -q * (2 * b) * CL_alpha * (1/2 + a)]])
        C_a = np.array([[q * (2 * b) * CL_alpha * 1/V, q * (2 * b) * CL_alpha_dot * b/V], [q * (2 * b) * CL_alpha * 1/V * (1/2 + a) * b, q * (2 * b) * b/V * b * CL_alpha * a * (1/2 - a)]])
        M_s = np.array([[m, S_theta], [S_theta, I_theta]])

        # print("K_s: ", K_s)
        # print("K_a: ", K_a)
        # print("C_a: ", C_a)
        # print("M_s: ", M_s)
        
        # Solve M_s * p^2 - C_a * p + (K_s - K_a) = 0
        # p1= (C_a + np.sqrt(C_a**2 - 4 * M_s @ (K_s - K_a)))/(2 * M_s)
        # p2= (C_a - np.sqrt(C_a**2 - 4 * M_s @ (K_s - K_a)))/(2 * M_s)
        # print("p1: ", p1)
        # print("p2: ", p2)
        # solution = np.array([p1, p2])
        
        # Define the matrices A and B
        A = -np.linalg.inv(M_s) @ C_a

        B = np.linalg.inv(M_s) @ (K_s - K_a)

        # Construct the system matrix F
        n = A.shape[0]
        I = np.eye(n)
        zero_matrix = np.zeros_like(A)

        F = np.block([
            [zero_matrix, I],
            [-B, -A]
        ])

        # Compute the eigenvalues
        eigenvalues, _ = eig(F)

        # Print the eigenvalues
        # print("Eigenvalues: ", eigenvalues)s



        # def equation(p):
        #     term1 = M_s * p**2  # M_s * p^2
        #     term2 = C_a * p  # C_a * p
        #     term3 = K_s - K_a  # K_s - K_a
        #     result = term1 - term2 + term3
        #     return result.flatten()  # Return a flattened array for fsolve

        # # Initial guess for the scalar p
        # p0 = 1

        # # Solve the equation
        # solution = fsolve(equation, p0)
        eigenvalues_lst.append(eigenvalues)
        print("Eigenvalues: ", eigenvalues)

    # Plot real parts of eigenvalue pairs in p against V
    real_eigenvalues_1 = [eigenvalue[0].real for eigenvalue in eigenvalues_lst]
    real_eigenvalues_2 = [eigenvalue[1].real for eigenvalue in eigenvalues_lst]
    # real_eigenvalues_3 = [eigenvalue[2].real for eigenvalue in eigenvalues_lst]
    # real_eigenvalues_4 = [eigenvalue[3].real for eigenvalue in eigenvalues_lst]
    plt.plot(V_arr, real_eigenvalues_1, label="Real part 1")
    plt.plot(V_arr, real_eigenvalues_2, label="Real part 2")
    # plt.plot(V_arr, real_eigenvalues_3, label="Real part 3")
    # plt.plot(V_arr, real_eigenvalues_4, label="Real part 4")
    plt.legend()
    plt.xlabel("V [m/s]")
    plt.ylabel("Re(p)")
    plt.title("Flutter diagram - Real parts of eigenvalues")
    plt.show()
    
    # plot imaginary parts of eigenvalue pairs in p against V
    # imaginary_eigenvalues_1 = [eigenvalue[0].imag for eigenvalue in eigenvalues_lst]
    # imaginary_eigenvalues_2 = [eigenvalue[1].imag for eigenvalue in eigenvalues_lst]
    imaginary_eigenvalues_3 = [eigenvalue[2].imag for eigenvalue in eigenvalues_lst]
    imaginary_eigenvalues_4 = [eigenvalue[3].imag for eigenvalue in eigenvalues_lst]
    # plt.plot(V_arr, imaginary_eigenvalues_1, label="Imaginary part 1")
    # plt.plot(V_arr, imaginary_eigenvalues_2, label="Imaginary part 2")
    plt.plot(V_arr, imaginary_eigenvalues_3, label="Imaginary part 3")
    plt.plot(V_arr, imaginary_eigenvalues_4, label="Imaginary part 4")
    plt.legend()
    plt.xlabel("V [m/s]")
    plt.ylabel("Im(p)")
    plt.title("Flutter diagram - Imaginary parts of eigenvalues")
    plt.show()

    return None



if __name__ == "__main__":
    # Mass parameters
    ma        = 1.567
    mf        = 0
    Icg_theta = 1
    Icg_beta = .01

    # Geometric parameters
    b       = 0.127
    a       = -0.5
    c       = .5
    x_theta = -.5
    x_beta = .1
    S       = c*1

    # Stiffness parameters
    K_h     = 2818.8
    K_theta = 37.3
    I_theta = ma*(x_theta*b)**2+mf*(c-a+x_beta)**2*b**2+Icg_theta+Icg_beta
    rho       = 1.225
    V         = 23
    alpha0    = 5*np.pi/180
    C_M_AC    = 0
    C_L_alpha = 2*np.pi


    flutter_diagram(K_h, K_theta, rho, np.linspace(np.finfo(float).eps, 150, 100), C_L_alpha, b, a, x_theta, ma+mf, I_theta)
    # DO NOT START THE LINSPACE AT 0, THIS WILL CAUSE A DIVISION BY ZERO ERROR IN THE C_a MATRIX



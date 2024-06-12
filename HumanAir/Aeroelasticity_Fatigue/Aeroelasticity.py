import numpy as np
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
# from HumanAir.aircraft_data import aircraft_data
# from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# from scipy.linalg import eig
from sympy import symbols, Matrix, det, Poly


"""
NOTE:
- We only analyze the part of the wing after the strut, so that we only have to deal with a cantilever beam.
We thus argue the strut blocks all vibration from occurring between it and the fuselage.
So the typical section is at 70% of the wing span after the strut.

- Nevermind, it is better not to do this, consider the full wing, and get conservative results :).
If we fail due to flutter, we can revise based on this, but untill then just stay conservative.

- The typical section is the airfoil cross-section of the wing at 70% of the wing span.

"""


def Divergence(K_theta, CL_alpha, a, B, rho):
    """
    Calculate the divergence speed of an aircraft based on the given parameters.
    The typical section is the airfoil cross-section of the wing at 70% of the wing span.

    Parameters
    ----------
    K_theta : float
        The torsional stiffness of the full 3D wing from root to typical section.
    CL_alpha : float
        The lift curve slope of the typical section airfoil.
    a : float
        Distance from half-chord to the elastic axis of the typical section airfoil.
        Non-dimensionalized by the half-chord length.
    B : float
        Half-chord length of the typical section (not wingspan!!!)
    rho : float
        The air density for the given operating regime.

    Returns
    -------
    V_div : float
        The divergence speed of the aircraft.
    """
    q = K_theta / ((2 * B) * CL_alpha * (1 / 2 + a) * B)  # where S = 2 * b since we analyse per unit span.
    V_div = np.sqrt(2 * q / rho)
    return V_div


def Reversal(h, K_theta, CL_alpha, CM_ac_beta, B, rho):
    """
    Calculate the flutter speed of an aircraft based on the given parameters.

    Parameters
    ----------
    h : float
        The distance from the half-chord to the hinge of the aileron
    K_theta : float
        The torsional stiffness of the full 3D wing from root to typical section.
    CL_alpha : float
        The lift curve slope of the typical section airfoil.
    CM_ac_beta : float
        The change in pitching moment about the typical section aerodynamic section
        with deflection of the control surface (aileron).
    B : float
        Half-chord length of the typical section (not wingspan!!!)
    rho : float
        The air density for the given operating regime.

    Returns
    -------
    V_rev : float
        The aileron reversal speed of the aircraft.
    """
    CL_beta = 2 * (
        np.sqrt(1 - h**2) + np.arccos(h)
    )  # The lift variation with deflection of the control surface (aileron).
    q = -CL_beta * K_theta / (CL_alpha * CM_ac_beta * 2 * B * (2 * B))
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
    altitude = (-np.log(rho / 1.225) * (287.05 * 288.15) / 9.80665) / (
        1 + np.log(rho / 1.225) * (287.05 * (-0.0065) / 9.80665)
    )
    return np.round(altitude, 1)


def flutter_analysis(m, I_theta, S_theta, rho_arr, K_h, K_theta, C_L_alpha, S, a, B, V_arr):
    """
    Calculate the flutter speed of an aircraft for various configurations (m_arr) and operating regimes (rho_arr).

    Parameters
    ----------
    m : float
        Mass of the airfoil = mass of the wing per unit span
        (divide the total mass of the wing by the wingspan to get this value).
    I_theta : float
        The torsional moment of inertia of the typical section about the elastic axis.
    S_theta : float
         Static moment related to the elastic axis (m * x_theta * b)
         (the moment due to the wing weight that constantly acts about the elastic axis of the wing).
    rho_arr : numpy array
        Array containing all operating regimes (air densities) you would like to analyse aircraft flutter for.
    K_h : float
        The bending stiffness of the full 3D wing from root to typical section.
    K_theta : float
        The torsional stiffness of the full 3D wing from root to typical section.
    C_L_alpha : float
        The lift curve slope of the typical section airfoil.
    S : float
        The wing surface area.
    a : float
        Distance from half-chord to the elastic axis of the typical section airfoil,
        non-dimensionalized by the half-chord length.
    B : float
        Half-chord length of the typical section (not wingspan!!!)
    V_arr : numpy array
        Array containing all airspeeds you would like to analyse aircraft flutter for. Should be 0 to V_dive * 1.15
        (+15% for certification, according to Dr. Roeland De Breuker)

    Returns
    -------
    Identifier : int
        The index of the configuration with the lowest flutter speed.
    V_flut_dict : dict
        Dictionary containing the flutter speed of all configurations analyzed,
        together with the relevant configuration information.
    """
    V_flut_dict = {}
    index = 0
    for rho in rho_arr:
        # Verbose
        if index % 2 == 0:
            print(f"Analysed {index}/{len(rho_arr)} configurations.")

        # Define matrices
        M_s = np.array([[m, S_theta], [S_theta, I_theta]])
        K_s = np.array([[K_h, 0], [0, K_theta]])
        K_a = np.array([[0, -S * C_L_alpha], [0, S * C_L_alpha * (0.5 + a) * B]])
        C_a = np.array([[-S * C_L_alpha, 0], [S * C_L_alpha * (0.5 * a) * B, 0]])
        # C_a = np.array([[-S*C_L_alpha, -S*C_L_alpha*(1-a)*B],
        # [S*C_L_alpha*(0.5*a)*B, S*(1-a)*C_L_alpha*B**2*(0.5+a) + S*B**2*C_L_alpha*a*(0.5-a)]])

        # Define the symbolic variable
        ps = symbols("ps")

        # Calculate eigenfrequencies with changing flow velocity
        q = 1 / 2 * rho * V_arr**2
        V = V_arr
        p_save = []

        for i in range(len(q)):
            Kae = Matrix(K_s - q[i] * K_a)
            A = ps**2 * Matrix(M_s) - ps * q[i] / V[i] * Matrix(C_a) + Kae
            DA = det(A)
            characteristic_equation = Poly(DA, ps).all_coeffs()
            p = np.roots([float(coeff) for coeff in characteristic_equation])
            # The first two eigenvalues are complex conjugates,
            # and so are the last two eigenvalues in the 4 eigenvalue array p.
            # For flutter analysis we are only interested in when the real part of the eigenvalue becomes positive,
            # and when the imaginary parts get close together.
            # The metric of becoming real is not affected by complex conjugates, and it does not matter
            # whether the imaginary parts come close together when negative or
            # when positive. All that is important is that we select the same conjugates
            # (that is [0,0] or [1,1], not [0,1] or [1,0] where 1 is the conjugate of one e-value.
            # We can therefore discard the second and fourth eigenvalues in the array p
            # (which are the conjugates of p[0] and p[2] respectively).
            p_save.append([p[0], p[2]])
        # Convert p_save to a numpy array for easier manipulation
        p_save = np.array(p_save)

        # Separate real parts of p_save, as these drive the flutter speed
        real_p_save = np.real(p_save)

        # Compute the point where the real part of the eigenvalue becomes positive for V != 0
        V_flut = 0
        for i in range(len(real_p_save)):
            V_flut1 = np.inf
            V_flut2 = np.inf
            if (
                real_p_save[i][0] > 0 and V[i] > 0.1
            ):  # 0.1 is arbitrary to avoid the first point where the real part is positive
                # If one of the eigenvalues really instantly gets a positive real part,
                # the np.round(x, 0) function i use at the end
                # accounts for this and rounds the V_flut back down to V_flut = 0, which it should be.
                # However, in case the eigenvalue
                # does not get positive real part instantly and just starts at 0
                # (which the program interperets as positive),
                # V_flut will not be set = 0 as well.
                # I now realize i also could have said real_p_save[i][0]>0.01, but i think my method is more precise
                # (what if 0.005, my method would catch this instantly).
                V_flut1 = V[i]
                break
            if (
                real_p_save[i][1] > 0 and V[i] > 0.1
            ):  # 0.1 is arbitrary to avoid the first point where the real part is positive
                V_flut2 = V[i]
                break
        if V_flut1 != 0 or V_flut2 != 0:
            V_flut = min(V_flut1, V_flut2)

        # Store Results
        V_flut_dict[index] = {"V_flut [m/s]": V_flut, "rho [kg/m^3]": rho}

        # Increment counter
        index += 1

    # Find the configuration with the lowest flutter speed
    V_flut_min = np.inf
    for key in V_flut_dict:
        if V_flut_dict[key]["V_flut [m/s]"] < V_flut_min:
            V_flut_min = V_flut_dict[key]["V_flut [m/s]"]
            Identifier = key
    if V_flut_min == np.inf:  # So if none of the eigenvalues got a positive real part after V = 0
        Identifier = None

    return Identifier, V_flut_dict


def flutter_diagram(m, I_theta, S_theta, rho, K_h, K_theta, C_L_alpha, S, a, B, V_arr):
    """
    Generate a flutter diagram for the the given aircraft configuration

    Parameters
    ----------
    m : float
        Mass of the airfoil = mass of the wing per unit span
        (divide the total mass of the wing by the wingspan to get this value).
    I_theta : float
        The torsional moment of inertia of the typical section about the elastic axis.
    S_theta : float
        Static moment related to the elastic axis (m * x_theta * b)
        (the moment due to the wing weight that constantly acts about the elastic axis of the wing).
    rho : float
        The air density for the given operating regime.
    K_h : float
        The bending stiffness of the full 3D wing from root to typical section.
    K_theta : float
        The torsional stiffness of the full 3D wing from root to typical section.
    C_L_alpha : float
        The lift curve slope of the typical section airfoil.
    S : float
        The wing surface area.
    a : float
        Distance from half-chord to the elastic axis of the typical section airfoil,
        non-dimensionalized by the half-chord length.
    B : float
        Half-chord length of the typical section (not wingspan!!!)
    V_arr : numpy array
        Array containing all airspeeds you would like to analyse aircraft flutter for. Should be 0 to V_dive * 1.15
        (+15% for certification, according to Dr. Roeland De Breuker)

    Returns
    -------
    None
    """
    # Define matrices
    M_s = np.array([[m, S_theta], [S_theta, I_theta]])
    K_s = np.array([[K_h, 0], [0, K_theta]])
    K_a = np.array([[0, -S * C_L_alpha], [0, S * C_L_alpha * (0.5 + a) * B]])
    C_a = np.array([[-S * C_L_alpha, 0], [S * C_L_alpha * (0.5 * a) * B, 0]])
    # C_a = np.array([[-S*C_L_alpha, -S*C_L_alpha*(1-a)*B],
    # [S*C_L_alpha*(0.5*a)*B, S*(1-a)*C_L_alpha*B**2*(0.5+a) + S*B**2*C_L_alpha*a*(0.5-a)]])

    # Define the symbolic variable
    ps = symbols("ps")

    # Calculate eigenfrequencies with changing flow velocity
    q = 1 / 2 * rho * V_arr**2
    V = V_arr
    p_save = []

    for i in range(len(q)):
        Kae = Matrix(K_s - q[i] * K_a)
        A = ps**2 * Matrix(M_s) - ps * q[i] / V[i] * Matrix(C_a) + Kae
        DA = det(A)
        characteristic_equation = Poly(DA, ps).all_coeffs()
        p = np.roots([float(coeff) for coeff in characteristic_equation])
        # The first two eigenvalues are complex conjugates,
        # and so are the last two eigenvalues in the 4 eigenvalue array p.
        # For flutter analysis we are only interested in when the real part of the eigenvalue becomes positive,
        # and when the imaginary parts get close together.
        # The metric of becoming real is not affected by complex conjugates,
        # and it does not matter whether the imaginary parts come close together when negative or
        # when positive. All that is important is that we select the same conjugates
        # (that is [0,0] or [1,1], not [0,1] or [1,0] where 1 is the conjugate of one e-value.
        # We can therefore discard the second and fourth eigenvalues in the array p
        # (which are the conjugates of p[0] and p[2] respectively).
        p_save.append([p[0], p[2]])
    # Convert p_save to a numpy array for easier manipulation
    p_save = np.array(p_save)

    # Separate real and imaginary parts of p_save
    real_p_save = np.real(p_save)
    imag_p_save = np.imag(p_save)

    # Compute the point where the real part of the eigenvalue becomes positive for V != 0
    V_flut = 0
    for i in range(len(real_p_save)):
        V_flut1 = np.inf
        V_flut2 = np.inf
        if (
            real_p_save[i][0] > 0 and V[i] > 0.1
        ):  # 0.1 is arbitrary to avoid the first point where the real part is positive
            V_flut1 = V[i]
            break
        if (
            real_p_save[i][1] > 0 and V[i] > 0.1
        ):  # 0.1 is arbitrary to avoid the first point where the real part is positive
            V_flut2 = V[i]
            break
    if V_flut1 != 0 or V_flut2 != 0:
        V_flut = min(V_flut1, V_flut2)

    # Plot flutter diagram: frequency vs dynamic pressure
    plt.figure()
    # plt.plot(q, imag_p_save, '.', linewidth=2)
    # plt.xlabel('Dynamic pressure [N/m^2]')
    plt.plot(V, imag_p_save, ".", linewidth=2)
    plt.xlabel("Airspeed [m/s]")
    plt.ylabel("Frequency ω [1/s]")
    plt.ylim(0, np.max(imag_p_save))
    plt.gca().tick_params(labelsize=14)
    plt.legend()
    plt.show()

    # Plot flutter diagram: damping vs dynamic pressure
    plt.figure()
    # plt.plot(q, real_p_save, '.', linewidth=2)
    # plt.xlabel('Dynamic pressure [N/m^2]')
    plt.plot(V, real_p_save, ".", linewidth=2)
    plt.xlabel("Airspeed [m/s]")
    plt.ylabel("Damping σ [1/s]")
    plt.axhline(0, color="k", linewidth=2)
    plt.gca().tick_params(labelsize=14)
    # plot the flutter speed
    if V_flut != 0:
        plt.plot(V_flut, 0, "ro")
        plt.text(V_flut * 1.05, -3, f"V_flut = {np.round(V_flut, 0)}", fontsize=12)
    plt.legend()
    plt.show()

    # Plot flutter diagram: imaginary part vs real part of p_save
    plt.figure()
    plt.plot(real_p_save, imag_p_save, "r.", linewidth=2)
    plt.xlabel("Damping σ [1/s]")
    plt.ylabel("Frequency ω [1/s]")
    plt.axvline(0, color="k", linewidth=2)
    plt.gca().tick_params(labelsize=14)
    plt.legend()
    plt.show()

    return None


def static_aeroelasticity(K_h, K_theta, S, C_L_alpha, a, B, q, alpha0, C_M_AC):
    """
    Analyze the static aeroelasticity of an aircraft wing.

    Parameters
    ----------
    K_h : float
        The bending stiffness of the full 3D wing from root to typical section.
    K_theta : float
        The torsional stiffness of the full 3D wing from root to typical section.
    S : float
        The wing surface area.
    C_L_alpha : float
        The lift curve slope of the typical section airfoil.
    a : float
        Distance from half-chord to the elastic axis of the typical section airfoil,
        non-dimensionalized by the half-chord length.
    B : float
        Half-chord length of the typical section (not wingspan!!!)
    q : float
        The dynamic pressure at which you would like to analyze the aeroelastic deformation of the wing.
    alpha0 : float
        The zero lift angle of attack of the typical section.
    C_M_AC : float
        The pitching moment coefficient about the aerodynamic center of the typical section.

    Returns
    -------
    None
    """
    # Define matrices
    K_s = np.array([[K_h, 0], [0, K_theta]])
    K_a = np.array([[0, -S * C_L_alpha], [0, S * C_L_alpha * (0.5 + a) * B]])
    F_st_alpha0 = np.array([[-S * C_L_alpha], [S * C_L_alpha * (0.5 + a) * B]])
    F_st_M_AC = np.array([[0], [S * C_M_AC * 2 * B]])

    # Iterative aeroelastic solution
    tol = 1
    iter = 0
    x_old = np.array([np.finfo(float).eps, np.finfo(float).eps])
    x_save = []

    while tol > 1e-6:
        iter += 1
        x_new = np.linalg.solve(K_s, q * F_st_alpha0 * (alpha0 + x_old[1]) + q * F_st_M_AC)
        tol = abs(1 - np.linalg.norm(x_new) / np.linalg.norm(x_old))
        x_old = x_new
        x_save.append(x_old)

    x_save = np.array(x_save).T[0]

    # Monolithic aeroelastic solution
    x_mono = np.linalg.solve(K_s - q * K_a, q * F_st_alpha0 * alpha0 + q * F_st_M_AC)

    print(x_save)
    # Plot convergence history
    plt.figure()
    plt.plot(range(1, iter + 1), x_save[1, :] * 180 / np.pi, linewidth=2)
    plt.plot([0, iter], [x_mono[1] * 180 / np.pi, x_mono[1] * 180 / np.pi], linewidth=2)
    plt.legend(["Partitioned", "Monolithic"])
    plt.xlabel("Number of iterations", fontsize=14)
    plt.ylabel("Twist rotation [deg]", fontsize=14)
    plt.axis([0, iter, 0, 7])
    plt.show()

    plt.figure()
    plt.plot(range(1, iter + 1), x_save[0, :], linewidth=2)
    plt.plot([0, iter], [x_mono[0], x_mono[0]], linewidth=2)
    plt.xlabel("Number of iterations", fontsize=14)
    plt.ylabel("Heave displacement [m]", fontsize=14)
    plt.legend(["Partitioned", "Monolithic"])
    plt.axis([0, iter, -5, 0])
    plt.show()


def static_trimmed_aeroelasticity(K_h, K_theta, S, C_L_alpha, a, V, q, W, C_M_AC):
    """
    Analyze the static trimmed aeroelasticity of an aircraft wing.

    Parameters
    ----------
    K_h : float
        The bending stiffness of the full 3D wing from root to typical section.
    K_theta : float
        The torsional stiffness of the full 3D wing from root to typical section.
    S : float
        The wing surface area.
    C_L_alpha : float
        The lift curve slope of the typical section airfoil.
    a : float
        Distance from half-chord to the elastic axis of the typical section airfoil,
        non-dimensionalized by the half-chord length.
    B : float
        Half-chord length of the typical section (not wingspan!!!)
    q : float
        The dynamic pressure at which you would like to analyze the aeroelastic deformation of the wing.
    W : float
        The weight of the total aircraft in [N]
    C_M_AC : float
        The pitching moment coefficient about the aerodynamic center of the typical section.

    Returns
    -------
    None
    """
    # Define matrices
    K_s = np.array([[K_h, 0], [0, K_theta]])
    K_a = np.array([[0, -S * C_L_alpha], [0, S * C_L_alpha * (0.5 + a) * V]])
    F_st_alpha0 = np.array([[-S * C_L_alpha], [S * C_L_alpha * (0.5 + a) * V]])
    F_st_M_AC = np.array([[0], [S * C_M_AC * 2 * V]])

    # Iterative aeroelastic trim solution
    tol1 = 1
    iter1 = 0
    alpha0_new = 1 * np.pi / 180
    alpha0_old = np.array([alpha0_new])
    x_old = np.array([[np.finfo(float).eps], [np.finfo(float).eps]])
    iter2_save = []
    alpha0_save = [alpha0_old]
    x_save = [x_old]

    while tol1 > 1e-6:
        tol2 = 1
        iter2 = 0
        iter1 += 1
        while tol2 > 1e-6:
            iter2 += 1
            x_new = np.linalg.solve(K_s, q * F_st_alpha0 * (alpha0_new + x_old[1]) + q * F_st_M_AC)
            tol2 = abs(1 - np.linalg.norm(x_new) / np.linalg.norm(x_old))
            x_old = x_new
            x_save.append(x_old)

        iter2_save.append(iter2)
        Lift = q * S * C_L_alpha * (x_new[1] + alpha0_old)
        alpha0_new = W / Lift * alpha0_old
        tol1 = abs(1 - alpha0_new / alpha0_old)
        alpha0_save.append(alpha0_new)
        alpha0_old = alpha0_new

    x_save = np.array(x_save).T
    alpha0_save = np.array(alpha0_save)

    # Monolithic aeroelastic solution
    U = np.block(
        [
            [K_s - q * K_a, np.array([[q * S * C_L_alpha], [-q * S * C_L_alpha * (0.5 + a) * V]])],
            [np.array([0, q * S * C_L_alpha]), q * S * C_L_alpha],
        ]
    )
    V = np.array([0, q * F_st_M_AC[1][0], W])
    x_mono = np.linalg.solve(U, V)

    # Plot convergence history
    plt.figure()
    plt.plot(range(len(alpha0_save)), alpha0_save * 180 / np.pi, linewidth=2)
    plt.plot([0, len(alpha0_save)], [x_mono[2] * 180 / np.pi, x_mono[2] * 180 / np.pi], linewidth=2)
    plt.xlabel("Number of iterations", fontsize=14)
    plt.ylabel("Rigid angle of attack [deg]", fontsize=14)
    plt.legend(["Partitioned", "Monolithic"])
    plt.gca().tick_params(labelsize=14)
    plt.axis([0, len(alpha0_save), 0, (np.max(alpha0_save)) * 180 / np.pi + 1])
    plt.show()

    plt.figure()
    plt.plot(range(len(x_save[0][1])), x_save[0][1, :] * 180 / np.pi, linewidth=2)
    plt.plot([0, len(x_save[0][1])], [x_mono[1] * 180 / np.pi, x_mono[1] * 180 / np.pi], linewidth=2)
    plt.xlabel("Number of iterations", fontsize=14)
    plt.ylabel("Twist rotation [deg]", fontsize=14)
    plt.legend(["Partitioned", "Monolithic"])
    plt.gca().tick_params(labelsize=14)
    plt.axis([0, len(x_save[0][1]), 0, 7])
    plt.show()

    plt.figure()
    plt.plot(range(len(x_save[0][0])), x_save[0][0, :], linewidth=2)
    plt.plot([0, len(x_save[0][0])], [x_mono[0], x_mono[0]], linewidth=2)
    plt.xlabel("Number of iterations", fontsize=14)
    plt.ylabel("Heave displacement [m]", fontsize=14)
    plt.legend(["Partitioned", "Monolithic"])
    plt.gca().tick_params(labelsize=14)
    plt.axis([0, len(x_save[0][0]), -5, 0])
    plt.show()


if __name__ == "__main__":
    analyze = (
        "flutter"  # "divergence", "reversal", "flutter", "static aeroelasticity" or "static trimmed aeroelasticity"
    )

    if analyze == "divergence":
        # Geometric Parameters
        a = -0.5
        B = 0.127

        # Stiffness Parameters
        K_theta = 37.3

        # Aerodynamic Parameters
        C_L_alpha = 2 * np.pi
        rho = 1.225

        # Calculate divergence speed
        V_div = Divergence(K_theta, C_L_alpha, a, B, rho)
        print(f"Divergence speed: {V_div} m/s")

    if analyze == "reversal":
        # Geometric Parameters
        h = 0.5  # TODO: change to correct value
        B = 0.127

        # Stiffness Parameters
        K_theta = 37.3

        # Aerodynamic Parameters
        C_L_alpha = 2 * np.pi
        CM_ac_beta = 0
        rho = 1.225

        # Calculate reversal speed
        V_rev = Reversal(h, K_theta, C_L_alpha, CM_ac_beta, B, rho)
        print(f"Reversal speed: {V_rev} m/s")

    if analyze == "flutter":
        # Mass Parameters
        m = 1.567 + 0
        I_theta = 0.01347
        S_theta = 0.08587

        # Geometric Parameters
        a = -0.5
        B = 0.127
        c = 0.5
        S = 1 * c

        # Stiffness Parameters
        K_h = 2818.8
        K_theta = 37.3

        # Aerodynamic Parameters
        rho = 1.225
        rho_arr = np.linspace(0.95, 1.225, 20)[::-1]
        C_L_alpha = 2 * np.pi

        # Velocity Profile:
        eps = np.finfo(float).eps  # Smallest positive number such that 1.0 + eps != 1.0
        # q = np.linspace(eps, 150, 100)
        # V_arr = np.sqrt(2*q/rho)
        V = 23
        V_arr = np.linspace(eps, 60, 400)  # np.linspace(eps, 15.6492159 , 100)
        q = 1 / 2 * rho * V**2

        # flutter_v2(m, I_theta, S_theta, rho, K_h, K_theta, C_L_alpha, S, a, b, V_arr)
        Identifier, V_flut_dict = flutter_analysis(
            m, I_theta, S_theta, rho_arr, K_h, K_theta, C_L_alpha, S, a, B, V_arr
        )

        # If no flutter:
        if Identifier is None:
            print("########################################################################################")
            print("GOOD NEWS: No flutter occurs over all mass configurations and operating regimes")
            print("########################################################################################")

        # If flutter:
        else:
            # Plot flutter diagram for this configuration
            print("########################################################################################")
            print("rho: ", V_flut_dict[Identifier]["rho [kg/m^3]"])
            print("Altitude: ", rho_altitude(V_flut_dict[Identifier]["rho [kg/m^3]"]))
            print("########################################################################################")
            flutter_diagram(
                m, I_theta, S_theta, V_flut_dict[Identifier]["rho [kg/m^3]"], K_h, K_theta, C_L_alpha, S, a, B, V_arr
            )

    if analyze == "static aeroelasticity":
        # Mass Parameters
        m = 1.567 + 0
        I_theta = 0.01347
        S_theta = 0.08587

        # Geometric Parameters
        a = -0.5
        B = 0.127
        c = 0.5
        S = 1 * c

        # Stiffness Parameters
        K_h = 2818.8
        K_theta = 37.3

        # Aerodynamic Parameters
        rho = 1.225
        C_L_alpha = 2 * np.pi
        alpha0 = 5 * np.pi / 180
        C_M_AC = 0

        # Velocity Profile:
        V = 23  # function only analyzes one velocity, one dynamic pressure
        q = 1 / 2 * rho * V**2

        # Static Aeroelasticity
        static_aeroelasticity(K_h, K_theta, S, C_L_alpha, a, B, q, alpha0, C_M_AC)

    if analyze == "static trimmed aeroelasticity":
        # Mass Parameters
        m = 1.567 + 0
        I_theta = 0.01347
        S_theta = 0.08587
        w = 30

        # Geometric Parameters
        a = -0.5
        B = 0.127
        c = 0.5
        S = 1 * c

        # Stiffness Parameters
        K_h = 2818.8
        K_theta = 37.3

        # Aerodynamic Parameters
        rho = 1.225
        C_L_alpha = 2 * np.pi
        alpha0 = 5 * np.pi / 180
        C_M_AC = 0

        # Velocity Profile:
        V = 23
        q = 1 / 2 * rho * V**2

        # Static Trimmed Aeroelasticity
        static_trimmed_aeroelasticity(K_h, K_theta, S, C_L_alpha, a, B, q, w, C_M_AC)

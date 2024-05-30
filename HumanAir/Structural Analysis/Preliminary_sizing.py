
def PrelimSizing(Vx, Vy, Mx, My, Mz, chords):
    # Material properties
    sigma = 324000000 # Pa yield strength 2024 Al
    tau = 283000000 # Pa yield shear strength 2024 Al
    tc = 0.15 # thickness over chord ratio of airfoil
    t = chords * tc

    Ixx = abs(Mx) * t / (2 * sigma)
    Iyy = abs(My) * chords / (2 * sigma)
    # assuming two area points at the thickness extremities:
    A = Ixx * 4 / t ** 2
    # determining the necessary spar thickness for shear stresses due to internal shear force:
    Q = 2 * A * 0.5 * t
    ts = abs(Vy) * Q / (Iyy * tau) # spar thickness needed to withstand shear stress due to bending

    # thickness of wingbox for torque
    twb = abs(Mz)/(2 * chords**2 * (0.4*0.15) * tau)

    #plt.plot(y_points, A)
    #plt.plot(y_points, Ixx)
    #plt.show()
    return Ixx, Iyy, ts, twb
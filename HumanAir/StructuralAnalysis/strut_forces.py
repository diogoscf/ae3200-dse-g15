import numpy as np

def read_points(L_cruise, W_cruise, W_fuel,idxs,load_factor):
    a = np.max(L_cruise)
    b = np.min(L_cruise)
    c = np.min(W_cruise)
    d = W_cruise[idxs[1]]
    e = d+W_fuel[idxs[0]]
    f = np.max(W_cruise)
    
    return a,b,c,d,e,f

def strut_force(E,L,h_fuselage,v_max,a,b,c,d,e,f,t_skin,h_spar1_MAC,h_spar2_MAC,x_spar1_MAC, x_spar2_MAC):
    #give different a,b,d,e,f for different load cases
    wa0 = a - b
    wb = b
    wc = c
    wd0 = d -c
    we = e - c
    wf0 = f - e

    """compatibility equation: 
            vmax = vA + vB - vC - vD - vE - v_strut """

    #inertia
    havg_spar_MAC = (h_spar1_MAC+h_spar2_MAC)/2
    I = 1/12*(x_spar2_MAC-x_spar1_MAC)*havg_spar_MAC**3 - 1/12*(x_spar2_MAC-x_spar1_MAC-2*t_skin)*(havg_spar_MAC-2*t_skin)**3
    #I = I*10e-12 #[m^4] if arguments are given in mm

    #deflections 
    vA = wa0*L**4/(30*E*I)
    vB = wb*L**4/(8*E*I)
    vC = wc*L**4/(8*E*I)
    vD = wd0*(0.6*L)**4/(30*E*I) +(0.4*L)*(0.6*L)**3/(24*E*I)
    vE = we*(0.4*L)**4/(8*E*I) +(0.6*L)*we*(0.4*L)**3/(6*E*I)
    vF = wf0*(0.4*L)**4/(30*E*I)+(0.6*L)*wf0*(0.4*L)**3/(24*E*I)
    v_strut_factor = (0.4*L)**3/(3*E*I)+(0.6*L)*(0.4*L)**2/(2*E*I)


    #Strut Force
    Vz_strut = (v_max-vA-vB+vC+vD+vE+vF)/v_strut_factor
    angle = np.arctan(h_fuselage/(0.4*L))
    Vy_strut = Vz_strut/np.tan(angle)

    return Vz_strut, Vy_strut
import numpy as np
import matplotlib.pyplot as plt


def fatigue_life(Sult, Smax, alpha, K_t = None):
    """
    Generate high-level estimate for the fatigue life of a material using the S-N curve.

    Parameters
    ----------
    Sult : float
        The ultimate tensile strength of the material in [MPa].
    Smax : float
        The maximum stress the aircraft should carry without plastic deformation, according to CS23 in [MPa].
    alpha : float
        Correction factor on the ultimate tensile strength, necessary for the final asymptote of the S-N curve.
    K_t : float
        The fatigue stress concentration factor.
    Returns
    -------
    Nf : float
        The fatigue life of the material.

    """
    # General parameters
    Sult = Sult * 1e6
    Smax = Smax * 1e6
    N_min = 100
    N_max = 1e6
    
    if Smax >= Sult or Smax <= alpha * Sult:
        raise ValueError('The maximum stress should be between 0 and the ultimate tensile strength.')


    if K_t == None:
        # From N = 0 to N = N_min
        x1 = np.linspace(0, np.log(100), 2)
        S1 = np.ones(len(x1))*np.log(Sult)
        
        # From N = N_min to N = N_max
        x2 = np.linspace(np.log(N_min), np.log(N_max), 2)
        S2 = np.log(Sult) + (np.log(alpha * Sult) - np.log(Sult))/(np.log(N_max) - np.log(N_min)) * (x2- np.log(N_min))
        
        # From N = N_max to N = infinity
        x3 = np.linspace(np.log(N_max), 20, 2)
        S3 = np.ones(len(x2))* np.log(alpha * Sult)

        # Compute intersection point
        x = (np.log(Smax) - np.log(Sult)) / ((np.log(alpha * Sult) - np.log(Sult))/(np.log(N_max) - np.log(N_min))) + (np.log(N_min))

        # Lifetime in number of cycles
        Nf = int(np.exp(x)/4) # divide by 4 to account for the coarseness of the approach

        # Plot the Basquin curve
        plt.figure()
        plt.plot(np.exp(x1), np.exp(S1), 'b')
        plt.plot(np.exp(x2), np.exp(S2), 'b')
        plt.plot(np.exp(x3), np.exp(S3), 'b')
        # Draw horizontal and vertical lines to intersection point, and note the values on the axes
        # Plot intersection point
        plt.plot(np.exp(x), Smax, 'ro')
        plt.plot([np.exp(x), np.exp(x)], [0, Smax], 'r--')
        plt.plot([0, np.exp(x)], [Smax, Smax], 'r--')
        # Plot the coodinate of the intersection point next to the point, with a slight offset using numerical values
        plt.text(np.exp(x), Smax, '({:.2e}, {:.2e})'.format(np.exp(x), Smax), fontsize=12, verticalalignment='bottom')
        # put the lifetime in a textbox on the top right corner
        plt.text(0.95, 0.95, f'lifetime: {Nf} flights', transform=plt.gca().transAxes, fontsize=14,
        verticalalignment='top', horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.5))
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Number of cycles')
        plt.ylabel('Stress')
        plt.title('Basquin curve')
        plt.show()
    
    if K_t is not None: 
        # From N = 0 to N = N_min
        x1 = np.linspace(0, np.log(100), 2)
        S1 = np.ones(len(x1))*np.log(Sult)
        
        # From N = N_min to N = N_max
        x2 = np.linspace(np.log(N_min), np.log(N_max), 2)
        S2 = np.log(Sult) + (np.log(alpha * Sult / K_t) - np.log(Sult))/(np.log(N_max) - np.log(N_min)) * (x2- np.log(N_min))
        
        # From N = N_max to N = infinity
        x3 = np.linspace(np.log(N_max), 20, 2)
        S3 = np.ones(len(x2))* np.log(alpha * Sult / K_t)

        # Compute intersection point
        x = (np.log(Smax) - np.log(Sult)) / ((np.log(alpha * Sult / K_t) - np.log(Sult))/(np.log(N_max) - np.log(N_min))) + (np.log(N_min))

        # Lifetime in number of cycles
        Nf = int(np.exp(x)/4) # divide by 4 to account for the coarseness of the approach

        # Plot the Basquin curve
        plt.figure()
        plt.plot(np.exp(x1), np.exp(S1), 'b')
        plt.plot(np.exp(x2), np.exp(S2), 'b')
        plt.plot(np.exp(x3), np.exp(S3), 'b')
        # Draw horizontal and vertical lines to intersection point, and note the values on the axes
        # Plot intersection point
        plt.plot(np.exp(x), Smax, 'ro')
        plt.plot([np.exp(x), np.exp(x)], [0, Smax], 'r--')
        plt.plot([0, np.exp(x)], [Smax, Smax], 'r--')
        # Plot the coodinate of the intersection point next to the point, with a slight offset using numerical values
        plt.text(np.exp(x), Smax, '({:.2e}, {:.2e})'.format(np.exp(x), Smax), fontsize=12, verticalalignment='bottom')
        # put the lifetime in a textbox on the top right corner
        plt.text(0.95, 0.95, f'lifetime: {Nf} flights', transform=plt.gca().transAxes, fontsize=14,
        verticalalignment='top', horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.5))
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Number of cycles')
        plt.ylabel('Stress')
        plt.title('Basquin curve')
        plt.show()



    test = False
    # if test: 
    #     # From N = 0 to N = N_min
    #     x1 = np.linspace(0, np.log(100), 2)
    #     S1 = np.ones(len(x1))*np.log(Sult)
        
    #     # From N = N_min to N = N_max
    #     x2 = np.linspace(np.log(N_min), np.log(N_max), 2)
    #     S2 = np.log(Sult) + (np.log(alpha * Sult) - np.log(Sult))/(np.log(N_max) - np.log(N_min)) * (x2- np.log(N_min))
        
    #     # From N = N_max to N = infinity
    #     x3 = np.linspace(np.log(N_max), 20, 2)
    #     S3 = np.ones(len(x2))* np.log(alpha * Sult)

    #     # Compute intersection point
    #     x = (np.log(Smax) - np.log(Sult)) / ((np.log(alpha * Sult) - np.log(Sult))/(np.log(N_max) - np.log(N_min))) + (np.log(N_min))

    #     # Plot the Basquin curve
    #     plt.figure()
    #     plt.plot(np.exp(x1), np.exp(S1), 'b')
    #     plt.plot(np.exp(x2), np.exp(S2), 'b')
    #     plt.plot(np.exp(x3), np.exp(S3), 'b')
    #     # Draw horizontal and vertical lines to intersection point, and note the values on the axes
    #     # Plot intersection point
    #     plt.plot(np.exp(x), Smax, 'ro')
    #     plt.plot([np.exp(x), np.exp(x)], [0, Smax], 'r--')
    #     plt.plot([0, np.exp(x)], [Smax, Smax], 'r--')
    #     # Plot the coodinate of the intersection point next to the point, with a slight offset using numerical values
    #     plt.text(np.exp(x), Smax, '({:.2e}, {:.2e})'.format(np.exp(x), Smax), fontsize=12, verticalalignment='bottom')

    #     plt.xscale('log')
    #     plt.yscale('log')
    #     plt.xlabel('Number of cycles')
    #     plt.ylabel('Stress')
    #     plt.title('Basquin curve')


    #     # From N = 0 to N = N_min
    #     x1 = np.linspace(0, np.log(100), 2)
    #     S1 = np.ones(len(x1))*np.log(Sult)
        
    #     # From N = N_min to N = N_max
    #     x2 = np.linspace(np.log(N_min), np.log(N_max), 2)
    #     S2 = np.log(Sult) + (np.log(alpha * Sult / K_t) - np.log(Sult))/(np.log(N_max) - np.log(N_min)) * (x2- np.log(N_min))
        
    #     # From N = N_max to N = infinity
    #     x3 = np.linspace(np.log(N_max), 20, 2)
    #     S3 = np.ones(len(x2))* np.log(alpha * Sult / K_t)

    #     # Compute intersection point
    #     x = (np.log(Smax) - np.log(Sult)) / ((np.log(alpha * Sult / K_t) - np.log(Sult))/(np.log(N_max) - np.log(N_min))) + (np.log(N_min))

    #     # Plot the Basquin curve
        
    #     plt.plot(np.exp(x1), np.exp(S1), 'b')
    #     plt.plot(np.exp(x2), np.exp(S2), 'b')
    #     plt.plot(np.exp(x3), np.exp(S3), 'b')
    #     # Draw horizontal and vertical lines to intersection point, and note the values on the axes
    #     # Plot intersection point
    #     plt.plot(np.exp(x), Smax, 'ro')
    #     plt.plot([np.exp(x), np.exp(x)], [0, Smax], 'r--')
    #     plt.plot([0, np.exp(x)], [Smax, Smax], 'r--')
    #     # Plot the coodinate of the intersection point next to the point, with a slight offset using numerical values
    #     plt.text(np.exp(x), Smax, '({:.2e}, {:.2e})'.format(np.exp(x), Smax), fontsize=12, verticalalignment='bottom')

    #     plt.xscale('log')
    #     plt.yscale('log')
    #     plt.xlabel('Number of cycles')
    #     plt.ylabel('Stress')
    #     plt.title('Basquin curve')
    #     plt.show()



if __name__ == '__main__':
    fatigue_life(500, 300, 0.3)
    fatigue_life(500, 300, 0.3, 2.5)
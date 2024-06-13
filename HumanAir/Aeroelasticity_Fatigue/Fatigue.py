import numpy as np
import matplotlib.pyplot as plt


def fatigue_life(Sult=None, alpha=None, Smax=None, K_t=None, verification=False, Experimental_SN=False):
    """
    Generate high-level estimate for the fatigue life of a material using the S-N curve.

    Parameters
    ----------
    Sult : float
        The ultimate tensile strength of the material in [MPa].
    alpha : float
        Correction factor on the ultimate tensile strength, necessary for the final asymptote of the S-N curve.
    Smax : float
        The maximum stress the aircraft should carry without plastic deformation, according to CS23 in [MPa].
    K_t : float
        The fatigue stress concentration factor.
    verification: bool
        If True, the verification on the function will be ran only
    experimental_SN: bool
        If True, the experimental data will be used to extract the S-N curve
    Returns
    -------
    Nf : float
        The fatigue life of the material.

    """
    if not verification:
        if not Experimental_SN:
            # General parameters
            Sult = Sult * 1e6
            if Smax is not None:
                Smax = Smax * 1e6
            N_min = 100
            N_max = 1e6

            if Smax is not None:
                if Smax >= Sult or Smax <= alpha * Sult:
                    raise ValueError("The maximum stress should be between 0 and the ultimate tensile strength.")

            if K_t is None:
                # From N = 0 to N = N_min
                x1 = np.linspace(0, np.log(100), 2)
                S1 = np.ones(len(x1)) * np.log(Sult)

                # From N = N_min to N = N_max
                x2 = np.linspace(np.log(N_min), np.log(N_max), 2)
                S2 = np.log(Sult) + (np.log(alpha * Sult) - np.log(Sult)) / (np.log(N_max) - np.log(N_min)) * (
                    x2 - np.log(N_min)
                )

                # From N = N_max to N = infinity
                x3 = np.linspace(np.log(N_max), 20, 2)
                S3 = np.ones(len(x2)) * np.log(alpha * Sult)

                # Merge all coordinates
                x = np.concatenate((x1, x2, x3))
                S = np.concatenate((S1, S2, S3))

                # Compute intersection point
                if Smax is not None:
                    i = (np.log(Smax) - np.log(Sult)) / (
                        (np.log(alpha * Sult) - np.log(Sult)) / (np.log(N_max) - np.log(N_min))
                    ) + (np.log(N_min))

                    # Lifetime in number of cycles
                    Nf = int(np.exp(i) / 4)  # divide by 4 to account for the coarseness of the approach

                # Plot the S-N curve
                plt.figure()
                plt.plot(np.exp(x), np.exp(S), "b")
                if Smax is not None:
                    # Plot intersection point
                    plt.plot(np.exp(i), Smax, "ro")
                    plt.plot([np.exp(i), np.exp(i)], [0, Smax], "r--")
                    plt.plot([0, np.exp(i)], [Smax, Smax], "r--")
                    # Plot the coodinate of the intersection point next to the point, with a slight offset using numerical values
                    plt.text(
                        np.exp(i),
                        Smax,
                        "({:.2e}, {:.2e})".format(np.exp(i), Smax),
                        fontsize=12,
                        verticalalignment="bottom",
                    )
                    # put the lifetime in a textbox on the top right corner
                    plt.text(
                        0.95,
                        0.95,
                        f"lifetime: {Nf} flights",
                        transform=plt.gca().transAxes,
                        fontsize=14,
                        verticalalignment="top",
                        horizontalalignment="right",
                        bbox=dict(facecolor="white", alpha=0.5),
                    )
                plt.xscale("log")
                plt.yscale("log")
                plt.xlabel("Number of cycles")
                plt.ylabel("Stress")
                plt.title("S-N curve")
                plt.show()

            if K_t is not None:
                # From N = 0 to N = N_min
                x1 = np.linspace(0, np.log(100), 2)
                S1 = np.ones(len(x1)) * np.log(Sult)

                # From N = N_min to N = N_max
                x2 = np.linspace(np.log(N_min), np.log(N_max), 2)
                S2 = np.log(Sult) + (np.log(alpha * Sult / K_t) - np.log(Sult)) / (np.log(N_max) - np.log(N_min)) * (
                    x2 - np.log(N_min)
                )

                # From N = N_max to N = infinity
                x3 = np.linspace(np.log(N_max), 20, 2)
                S3 = np.ones(len(x2)) * np.log(alpha * Sult / K_t)

                # Merge all coordinates
                x = np.concatenate((x1, x2, x3))
                S = np.concatenate((S1, S2, S3))

                # Compute intersection point
                if Smax is not None:
                    i = (np.log(Smax) - np.log(Sult)) / (
                        (np.log(alpha * Sult) - np.log(Sult)) / (np.log(N_max) - np.log(N_min))
                    ) + (np.log(N_min))

                    # Lifetime in number of cycles
                    Nf = int(np.exp(i) / 4)  # divide by 4 to account for the coarseness of the approach

                # Plot the S-N curve
                plt.figure()
                plt.plot(np.exp(x), np.exp(S), "b")
                if Smax is not None:
                    # Plot intersection point
                    plt.plot(np.exp(i), Smax, "ro")
                    plt.plot([np.exp(i), np.exp(i)], [0, Smax], "r--")
                    plt.plot([0, np.exp(i)], [Smax, Smax], "r--")
                    # Plot the coodinate of the intersection point next to the point, with a slight offset using numerical values
                    plt.text(
                        np.exp(i),
                        Smax,
                        "({:.2e}, {:.2e})".format(np.exp(i), Smax),
                        fontsize=12,
                        verticalalignment="bottom",
                    )
                    # put the lifetime in a textbox on the top right corner
                    plt.text(
                        0.95,
                        0.95,
                        f"lifetime: {Nf} flights",
                        transform=plt.gca().transAxes,
                        fontsize=14,
                        verticalalignment="top",
                        horizontalalignment="right",
                        bbox=dict(facecolor="white", alpha=0.5),
                    )
                plt.xscale("log")
                plt.yscale("log")
                plt.xlabel("Number of cycles")
                plt.ylabel("Stress")
                plt.title("S-N curve")
                plt.show()

        if Experimental_SN:
            Smax = Smax * 1e6
            data = np.array(
                [
                    (10.0, 482.63),
                    (20.0, 482.63),
                    (50.0, 482.63),
                    (70.0, 482.63),
                    (100.0, 420.28),
                    (200.0, 325.43),
                    (500.0, 241.32),
                    (1000.0, 198.91),
                    (2000.0, 168.92),
                    (5000.0, 142.3),
                    (7000.0, 135.83),
                    (10000.0, 120.66),
                    (20000.0, 99.46),
                    (50000.0, 80.64),
                    (100000.0, 71.15),
                    (200000.0, 64.45),
                    (500000.0, 58.5),
                    (1000000.0, 55.5),
                    (2000000.0, 53.38),
                    (5000000.0, 51.5),
                    (10000000.0, 50.55),
                    (20000000.0, 49.88),
                    (50000000.0, 49.28),
                    (100000000.0, 48.99),
                    (200000000.0, 48.77),
                    (500000000.0, 48.59),
                    (1000000000.0, 48.49),
                ]
            )
            data[:, 1] = data[:, 1] * 1e6
            if Smax is not None:
                if Smax >= data[0, 1] or Smax <= data[-1, 1]:
                    raise ValueError("The maximum stress should be between 0 and the ultimate tensile strength.")

            # Linearly interpolate between curve points to get the Nf for the given Smax
            if Smax is not None:
                for i in range(len(data)):
                    if Smax <= data[i, 1] and Smax >= data[i + 1, 1]:
                        Nf = int(
                            np.exp(
                                np.log(data[i, 0])
                                + (np.log(data[i + 1, 0]) - np.log(data[i, 0]))
                                / (np.log(data[i + 1, 1]) - np.log(data[i, 1]))
                                * np.log((Smax / data[i, 1]))
                            )
                        )
                        break

            # plot the experimental data curve
            plt.figure()
            plt.plot(data[:, 0], data[:, 1])
            plt.xscale("log")
            plt.yscale("log")
            plt.xlabel("Number of cycles")
            plt.ylabel("Stress")
            plt.title("S-N curve")
            if Smax is not None:
                plt.plot(Nf, Smax, "ro")
                plt.plot([Nf, Nf], [0, Smax], "r--")
                plt.plot([0, Nf], [Smax, Smax], "r--")
                # Plot the coodinate of the intersection point next to the point, with a slight offset using numerical values
                plt.text(Nf, Smax, "({:.2e}, {:.2e})".format(Nf, Smax), fontsize=12, verticalalignment="bottom")
                # put the lifetime in a textbox on the top right corner
                plt.text(
                    0.95,
                    0.95,
                    f"lifetime: {Nf} flights",
                    transform=plt.gca().transAxes,
                    fontsize=14,
                    verticalalignment="top",
                    horizontalalignment="right",
                    bbox=dict(facecolor="white", alpha=0.5),
                )
            plt.show()

    # test = False
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

    if verification:
        data = np.array(
            [
                (10.0, 482.63),
                (20.0, 482.63),
                (50.0, 482.63),
                (70.0, 482.63),
                (100.0, 420.28),
                (200.0, 325.43),
                (500.0, 241.32),
                (1000.0, 198.91),
                (2000.0, 168.92),
                (5000.0, 142.3),
                (7000.0, 135.83),
                (10000.0, 120.66),
                (20000.0, 99.46),
                (50000.0, 80.64),
                (100000.0, 71.15),
                (200000.0, 64.45),
                (500000.0, 58.5),
                (1000000.0, 55.5),
                (2000000.0, 53.38),
                (5000000.0, 51.5),
                (10000000.0, 50.55),
                (20000000.0, 49.88),
                (50000000.0, 49.28),
                (100000000.0, 48.99),
                (200000000.0, 48.77),
                (500000000.0, 48.59),
                (1000000000.0, 48.49),
            ]
        )
        data[:, 1] = data[:, 1] * 1e6
        Sult = 500 * 1e6
        alpha = 0.1
        N_min = 100
        N_max = 1e6
        # From N = 0 to N = N_min
        x1 = np.linspace(0, np.log(100), 2)
        S1 = np.ones(len(x1)) * np.log(Sult)

        # From N = N_min to N = N_max
        x2 = np.linspace(np.log(N_min), np.log(N_max), 2)
        S2 = np.log(Sult) + (np.log(alpha * Sult) - np.log(Sult)) / (np.log(N_max) - np.log(N_min)) * (
            x2 - np.log(N_min)
        )

        # From N = N_max to N = infinity
        x3 = np.linspace(np.log(N_max), 20, 2)
        S3 = np.ones(len(x2)) * np.log(alpha * Sult)

        # Merge all coordinates
        x = np.concatenate((x1, x2, x3))
        S = np.concatenate((S1, S2, S3))

        # Plot the S-N curve
        plt.figure()
        plt.plot(np.exp(x), np.exp(S), "b")
        # plot the experimental data curve
        plt.plot(data[:, 0], data[:, 1])
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Number of cycles")
        plt.ylabel("Stress")
        plt.title("S-N curve")
        plt.legend(["Basquin law approximation", "Experimental data"])
        plt.show()





if __name__ == "__main__":
    fatigue_life(Sult=500, alpha=0.1, Smax=120, verification=False, Experimental_SN=False)
    fatigue_life(Sult=500, alpha=0.1, Smax=120, verification=False, Experimental_SN=True)
    fatigue_life(verification=True)

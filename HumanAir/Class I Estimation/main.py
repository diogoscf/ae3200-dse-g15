import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures



"MTOW vs Payload"
def plot_MTOW_vs_Payload(combined = True, show = None, save = None):

    # single propeller
    MTOW_single_propeller = np.array([680.388, 997.903, 1065.942, 1202.02, 1202.02, 952.544, 1315.418, 1315.418, 1339.912]) # kg
    Payload_single_propeller = np.array([180.529, 318.421, 326.1320, 332.9360, 324.3180, 145.6030, 494.4150, 383.2850, 335.6580]) # kg

    # double propeller
    MTOW_double_propeller = np.array([3342.068582, 4399.845989, 3243.185446, 1814.36948, 2190.851147, 3810.175908, 3077.62423, 4501.904272, 4628.910136]) # kg
    Payload_double_propeller = np.array([583.7733802, 668.5951534, 594.6595971, 332.4832072, 537.9605508, 861.825503, 544.310844, 725.747792, 889.0410452]) # kg

    # combine data points for plotting
    if combined:
        MTOW_ = np.concatenate((MTOW_single_propeller, MTOW_double_propeller))
        Payload_ = np.concatenate((Payload_single_propeller, Payload_double_propeller))
        MTOW_ = [MTOW_]
        Payload_ = [Payload_]

    else:
        MTOW_ = [MTOW_single_propeller, MTOW_double_propeller]
        Payload_ = [Payload_single_propeller, Payload_double_propeller]

    # keep track of the combined parameter in order to know how to save the images
    step = 0

    for MTOW in MTOW_:
        for Payload in Payload_:

            # plot data points
            plt.scatter(Payload, MTOW)
            x_data = np.linspace(min(Payload), max(Payload), 10000000)

            # show linear regression line
            regressor = LinearRegression()
            regressor.fit(Payload.reshape(-1, 1), MTOW)
            plt.plot(x_data, regressor.predict(x_data.reshape(-1,1)), color='red', alpha=0.5, markersize=1)


            # show polynomial regression line
            poly_reg = PolynomialFeatures(degree=2)

            # fit_transform() takes the input and returns the transformed input
            payload_poly = poly_reg.fit_transform(Payload.reshape(-1, 1))
            x_data_poly = poly_reg.fit_transform(x_data.reshape(-1, 1))

            # fit the transformed input to the linear regression model
            regressor_poly = LinearRegression()
            regressor_poly.fit(payload_poly, MTOW)
            plt.plot(x_data, regressor_poly.predict(x_data_poly), color='green', alpha=0.5)


            # show exponential regression line
            regressor_exp = LinearRegression()
            regressor_exp.fit(Payload.reshape(-1, 1), np.log(MTOW))
            plt.plot(x_data, np.exp(regressor_exp.predict(x_data.reshape(-1, 1))), color='blue', alpha=0.5)


            # show logarithmic regression line
            regressor_log = LinearRegression()
            regressor_log.fit(np.log(Payload).reshape(-1, 1), MTOW)
            plt.plot(x_data, regressor_log.predict(np.log(x_data).reshape(-1, 1)), color='orange', alpha=0.5)

            # adding the grid
            plt.grid()


            # adding labels
            plt.ylabel('MTOW')
            plt.xlabel('Payload')

            if show:
                plt.show()

            # save the images based on the case
            if save:
                step += 1
                if step == 1 and combined:
                    plt.savefig("Combined_MTOW_vs_Payload.svg")

                elif step == 1 and not combined:
                    plt.savefig("Single_MTOW_vs_Payload.svg")

                elif step == 4 and not combined:
                    plt.savefig("Double_MTOW_vs_Payload.svg")


plot_MTOW_vs_Payload(combined = True, show = True)
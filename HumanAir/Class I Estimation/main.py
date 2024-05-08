import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
"MTOW vs Payload"

# single propeller
MTOW_single_propeller = np.array([680.388, 997.903, 1065.942, 1202.02, 1202.02, 952.544, 1315.418, 1315.418, 1339.912]) # kg
Payload_single_propeller = np.array([180.529, 318.421, 326.1320, 332.9360, 324.3180, 145.6030, 494.4150, 383.2850, 335.6580]) # kg

# double propeller
MTOW_double_propeller = np.array([3342.068582, 4399.845989, 3243.185446, 1814.36948, 2190.851147, 3810.175908, 3077.62423, 4501.904272, 4628.910136]) # kg
Payload_double_propeller = np.array([583.7733802, 668.5951534, 594.6595971, 332.4832072, 537.9605508, 861.825503, 544.310844, 725.747792, 889.0410452]) # kg

# combine data points for plotting

MTOW = np.concatenate((MTOW_single_propeller, MTOW_double_propeller))
Payload = np.concatenate((Payload_single_propeller, Payload_double_propeller))

# plot data points
plt.scatter(MTOW, Payload)
x_data = np.linspace(min(MTOW), max(MTOW), 1000000)

# show linear regression line
regressor = LinearRegression()
regressor.fit(MTOW.reshape(18, 1), Payload)
plt.plot(x_data, regressor.predict(x_data.reshape(-1,1)), color='red', alpha=0.5)


# show polynomial regression line
poly_reg = PolynomialFeatures(degree=2)
MTOW_poly = poly_reg.fit_transform(MTOW.reshape(18, 1))
x_data_poly = poly_reg.fit_transform(x_data.reshape(-1, 1))
regressor_poly = LinearRegression()
regressor_poly.fit(MTOW_poly, Payload)
plt.plot(x_data, regressor_poly.predict(x_data_poly), color='green', alpha=0.5)


# show exponential regression line
regressor_exp = LinearRegression()
regressor_exp.fit(MTOW.reshape(18, 1), np.log(Payload))
plt.plot(x_data, np.exp(regressor_exp.predict(x_data.reshape(-1, 1))), color='blue', alpha=0.5)


# show logarithmic regression line
regressor_log = LinearRegression()
regressor_log.fit(np.log(MTOW).reshape(-1, 1), Payload)
plt.plot(x_data, regressor_log.predict(np.log(x_data).reshape(-1, 1)), color='purple', alpha=0.5)
plt.show()
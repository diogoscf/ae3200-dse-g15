import sys
import os
import json
import numpy as np

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
sys.path.append(project_root)

from HumanAir.LoadingDiagram.Parameters import Parameters_ConvNoCanard as p
from Class_I_Weight_Estimation import WeightEstm as WeightEstimation
from HumanAir.LoadingDiagram.Main import WP_WS
from HumanAir.CO2Calculator.conceptual_co2 import calculate_co2_reduction as co2



# loading the json file
with open('../Configurations/design.json', 'r') as f:
    dict = json.load(f)

# initialise the loading diagram to get W/P and W/S
dict['W/P'], dict['W/S'] = WP_WS().calculate_optimal_point()
print(dict['W/P'], dict['W/S'])

# initialise the class I weight estimation
WeightEstm = WeightEstimation(dict)

# initialise the bat percentage
bat = np.arange(0, 0.18, 0.001)
coeff_exp, coeff_pol = WeightEstm.PolynomialRegression(bat)

# run the regression to find power required cruise
for step in range(len(bat)):

    # calculate the power required cruise
    dict['P_req_cruise_W'] = dict['P_cruise/P_TO'] * np.exp(coeff_exp[1]) * np.exp(coeff_exp[0] * bat[step])

    dict['E_bat_Wh'] = dict['P_req_cruise_W'] * dict['endurance'] / dict['P_cruise/P_TO'] * bat[step]
    print("STEP:", step)
    print("Battery Percentage:", bat[step]*100)
    print("Energy:",dict["E_bat_Wh"])
    print("Power required cruise",dict['P_req_cruise_W'])
    # print(bat[step]*100)
    co2_ratio = co2(ac_data=dict)
    #print(np.round(co2_ratio*100,2))
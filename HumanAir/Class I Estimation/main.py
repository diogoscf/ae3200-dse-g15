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



# loading the json file
with open('../Configurations/conventional - Nicholas.json', 'r') as f:
    dict = json.load(f)

# initialise the loading diagram to get W/P and W/S
dict['W/P'], dict['W/S'] = WP_WS().calculate_optimal_point()

# initialise the class I weight estimation
WeightEstm = WeightEstimation(dict)

# return the weight estimation values
MTOW, MTOW_cont, OEW_cont, PowertrainWeight_cont, BatteryWeight_cont, FuelWeight_cont, WingWeight_cont, Payload_cont = WeightEstm.Iterations(bat[1])



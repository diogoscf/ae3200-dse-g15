import sys
import os
import json
import numpy as np
#from tqdm import tqdm

# Integration in progress v2
# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(script_dir,'..', '..'))
sys.path.append(project_root)

from HumanAir.LoadingDiagram.Parameters import Parameters_ConvNoCanard as p
from Class_I_Weight_Estimation import WeightEstm as WeightEstimation
from HumanAir.LoadingDiagram.Main import WP_WS
from HumanAir.CO2_Calculator.conceptual_co2 import calculate_co2_reduction as co2


# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Construct the absolute path to the design.json file
design_json_path = os.path.join(script_dir, 'Configurations', 'design.json')

# Print the absolute path for debugging
print(f"Looking for design.json at: {os.path.abspath(design_json_path)}")
#C:\Users\nicho\Documents\GitHub\ae3200-dse-g15\HumanAir\Configurations # use this for vs code
#'../Configurations/design.json', 'r' # use this for pycharm

# Attempt to open the file
with open('../Configurations/design.json', 'r') as f:
    dict = json.load(f)



def Generate(p, dict, run=False):
    A_lst = np.arange(7.5, 9.51, 0.5)
    eta_p_lst = np.arange(0.8, 0.851, 0.05)
    Clmax_clean_lst = np.arange(1.6, 2.21, 0.2)
    Clmax_TO_lst = np.arange(2, 2.61, 0.2)
    Clmax_Land_lst = np.arange(2, 2.61, 0.2)
    Cd0_lst = np.arange(0.026, 0.0281, 0.002)
    V_cruise_lst = np.arange(60, 63.1, 1)
    climbrate_lst = np.arange(3.5, 4.51, 0.5)

    total_iterations = (len(A_lst) * len(eta_p_lst) * len(Clmax_clean_lst) *
                        len(Clmax_TO_lst) * len(Clmax_Land_lst) * len(Cd0_lst) *
                        len(V_cruise_lst) * len(climbrate_lst))

    current_iteration = 0
    if run:
        idx=-1
        dict_iterations= {}
        for A in A_lst:
            for eta_p in eta_p_lst:
                for Clmax_clean in Clmax_clean_lst:
                    for Clmax_TO in Clmax_TO_lst:
                        for Clmax_Land in Clmax_Land_lst:
                            for Cd0 in Cd0_lst:
                                for V_cruise in V_cruise_lst:
                                    
                                    dict["endurance"]=1111200/V_cruise/3600
                                    for climbrate in climbrate_lst:
                                        current_iteration+=1

                                        if current_iteration%500==0: print("Iteration:"+str(current_iteration)+"/"+str(total_iterations))
                                        
                                        p.A = A
                                        p.eta_p = eta_p
                                        p.Clmax_clean = Clmax_clean
                                        p.Clmax_TO = Clmax_TO
                                        p.Clmax_Land = Clmax_Land
                                        p.Cdo = Cd0
                                        p.V_cruise = V_cruise
                                        p.climbrate = climbrate
                                        

                                        # initialise the loading diagram to get W/P and W/S
                                        dict['W/P'], dict['W/S'] = WP_WS().calculate_optimal_point()

                                        # initialise the class I weight estimation
                                        WeightEstm = WeightEstimation(dict)

                                        # initialise the bat percentage
                                        bat = np.arange(0, 0.18, 0.001)
                                        coeff_exp, coeff_pol = WeightEstm.PolynomialRegression(bat)

                                        # run the regression to find power required cruise
                                        co2_ratio_max=0
                                        
                                        ok=0
                                        for step in range(len(bat)):

                                            # calculate the power required cruise
                                            dict['P_req_cruise_W'] = dict['P_cruise/P_TO'] * np.exp(coeff_exp[1]) * np.exp(coeff_exp[0] * bat[step])

                                            dict['E_bat_Wh'] = dict['P_req_cruise_W'] * dict['endurance'] / dict['P_cruise/P_TO'] * bat[step]

                                            
                                            co2_ratio = co2(ac_data=dict)
                                            if co2_ratio * 100 > 50 and ok==0: 
                                                
                                                print(co2_ratio*100)
                                                idx+=1
                                                ok=1



                                            if co2_ratio * 100 > co2_ratio_max and co2_ratio * 100 >50 and dict['E_bat_Wh']<250000:
                                                CO2 = co2_ratio

                                                dict_iterations[str(idx)] = {}
                                                dict_iterations[str(idx)]['A'] = A
                                                dict_iterations[str(idx)]['eta_p'] = eta_p
                                                dict_iterations[str(idx)]['Clmax_clean'] = Clmax_clean
                                                dict_iterations[str(idx)]['Clmax_TO'] = Clmax_TO
                                                dict_iterations[str(idx)]['Clmax_Land'] = Clmax_Land
                                                dict_iterations[str(idx)]['Cd0'] = Cd0
                                                dict_iterations[str(idx)]['V_cruise'] = V_cruise
                                                dict_iterations[str(idx)]['climbrate'] = climbrate
                                                dict_iterations[str(idx)]['CO2'] = np.round(CO2*100, 2)
                                                dict_iterations[str(idx)]['W/P'] = dict['W/P']
                                                dict_iterations[str(idx)]['W/S'] = dict['W/S']
                                                dict_iterations[str(idx)]['bat'] = bat[step]

                                                co2_ratio_max=co2_ratio

        file_name = 'data_iterations.json'
        with open(file_name, 'w') as file:
            json.dump(dict_iterations, file, indent=4)
                                    
                                    

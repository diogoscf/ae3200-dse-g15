import sys
import os
import json
import numpy as np
import logging
import colorlog


"Dear Programmer Please do not remove this line, it is very important for the correct function of the main program"

handler = colorlog.StreamHandler()
handler.setFormatter(colorlog.ColoredFormatter(
    "%(log_color)s%(levelname)s:%(message)s",
    log_colors={
        'DEBUG': 'cyan',
        'INFO': 'green',
        'WARNING': 'yellow',
        'ERROR': 'red',
        'CRITICAL': 'bold_red',
    }
))
logger = colorlog.getLogger()
logger.addHandler(handler)
logger.setLevel(logging.INFO)

# initialise the logging


# Integration in progress v2
logging.info(' Starting the program')

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(script_dir,'..', '..'))
sys.path.append(project_root)

from HumanAir.LoadingDiagram.Parameters import Parameters_ConvNoCanard as p
from Class_I_Weight_Estimation import WeightEstm as WeightEstimation
from HumanAir.LoadingDiagram.Main import WP_WS
from HumanAir.CO2_Calculator.conceptual_co2 import calculate_co2_reduction as co2
from HumanAir.Weights_and_CG.weight_fractions import find_lg, iterate_cg_lg
from HumanAir.AerodynamicDesign.Aerodynamics_Main import aerodynamic_design

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Construct the absolute path to the design.json file
design_json_path = os.path.join(script_dir,'..', 'Configurations', 'design.json')

# Print the absolute path for debugging
logging.info(f" Looking for design.json at: {os.path.abspath(design_json_path)}")


#"c:\\Users\\nicho\\Documents\\GitHub\\ae3200-dse-g15\\HumanAir\\Configurations\\design.json" # use this for vs code
#'../Configurations/design.json', 'r' # use this for pycharm

# Attempt to open the file
with open(design_json_path, 'r') as f:
    dict = json.load(f)

logging.info(" Opening design.json successful")


def Generate(p, dict, run=False):

    # tune the parameters with a reasonable range
    A_lst = np.arange(7.5, 9.51, 0.5)
    eta_p_lst = np.arange(0.8, 0.851, 0.05)
    Clmax_clean_lst = np.arange(1.6, 2.21, 0.2)
    Clmax_TO_lst = np.arange(2, 2.61, 0.2)
    Clmax_Land_lst = np.arange(2, 2.61, 0.2)
    Cd0_lst = np.arange(0.026, 0.0281, 0.002)
    V_cruise_lst = np.arange(60, 63.1, 1)
    climbrate_lst = np.arange(3.5, 4.51, 0.5)

    # calculate the total numbers of iterations
    total_iterations = (len(A_lst) * len(eta_p_lst) * len(Clmax_clean_lst) *
                        len(Clmax_TO_lst) * len(Clmax_Land_lst) * len(Cd0_lst) *
                        len(V_cruise_lst) * len(climbrate_lst))

    # initialise the iteration counter
    current_iteration = 0

    # run condition to not run the loop by mistake
    if run:
        idx = -1
        dict_iterations= {}
        for A in A_lst:
            for eta_p in eta_p_lst:
                for Clmax_clean in Clmax_clean_lst:
                    for Clmax_TO in Clmax_TO_lst:
                        for Clmax_Land in Clmax_Land_lst:
                            for Cd0 in Cd0_lst:
                                for V_cruise in V_cruise_lst:

                                    # update the endurance based on v cruise
                                    dict["endurance"]=1111200/V_cruise/3600

                                    for climbrate in climbrate_lst:
                                        current_iteration+=1

                                        # print the iteration number every 500 steps
                                        if current_iteration%500==0: print("Iteration:"+str(current_iteration)+"/"+str(total_iterations))

                                        # update the parameters
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
                                        co2_ratio_max = 0

                                        # set the condition to find the first point with co2 ratio > 50
                                        ok = 0
                                        for step in range(len(bat)):

                                            # calculate the power required cruise
                                            dict['P_req_cruise_W'] = dict['P_cruise/P_TO'] * np.exp(coeff_exp[1]) * np.exp(coeff_exp[0] * bat[step])
                                            dict['E_bat_Wh'] = dict['P_req_cruise_W'] * dict['endurance'] / dict['P_cruise/P_TO'] * bat[step]

                                            # calculate the co2 ratio for the specific combination of parameters
                                            co2_ratio = co2(ac_data=dict)

                                            if co2_ratio * 100 > 50 and ok == 0:

                                                idx += 1
                                                ok = 1



                                            if co2_ratio * 100 > co2_ratio_max and co2_ratio * 100 >50 and dict['E_bat_Wh']<250000: # the <250000 condition is for the battery to be able to be charged
                                                CO2 = co2_ratio

                                                # save the combination of parameters
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

        # save the json file with all possible design options
        file_name = 'data_iterations.json'
        with open(file_name, 'w') as file:
            json.dump(dict_iterations, file, indent=4)
                                    
def calculate_weighted_score(point_data, weights):

    # setting the worst optimal data
    worst_optimal_data = {'A': 9.5,
                    'eta_p': 0.85,
                    'Clmax_clean': 2.2,
                    'Clmax_TO': 2.6,
                    'Clmax_Land': 2.6,
                    'Cd0': 0.026,
                    'V_cruise': 60,
                    'climbrate': 3.5,
                    'bat': 0.12,
                    'CO2': 50}

    # initialing the score
    score = 0

    # calculate the score
    for key, weight in weights.items():
        score += (worst_optimal_data[key] - point_data.get(key, 0))/worst_optimal_data[key] * weight # get a normalised value for the contribution of each key
    return score

def find_optimal_design(maximum_weight_battery = 1000, weights = None, CO2_threshold = 50, design_points = None, printing = False):

    optimum_design_points = {key: value for key, value in design_points.items() if value['CO2'] > CO2_threshold}
    scores = [(key, calculate_weighted_score(value, weights)) for key, value in optimum_design_points.items()]

    sorted_design_points = sorted(scores, key=lambda x: x[1], reverse=True)

    step = 0
    value = np.inf
    minimum = np.inf

    while value > maximum_weight_battery:

        optimum_design_option = design_points[sorted_design_points[step][0]]

        # updating the design.json file with the optimum design option
        dict['Aero']['AR'] = optimum_design_option['A']
        dict['Power_prop']['eta_p'] = optimum_design_option['eta_p']
        dict['Aero']['CLmax_clean'] = optimum_design_option['Clmax_clean']
        dict['Aero']['CLmax_TO'] = optimum_design_option['Clmax_TO']
        dict['Aero']['CLmax_Land'] = optimum_design_option['Clmax_Land']
        dict['Aero']['CD0'] = optimum_design_option['Cd0']
        dict['Performance']['Vc_m/s'] = optimum_design_option['V_cruise']
        dict['Performance']['climb_rate'] = optimum_design_option['climbrate']
        dict['Power_prop']['bat'] = optimum_design_option['bat']
        dict['Performance']['W/S_N/m2'] = optimum_design_option['W/S']
        dict['Performance']['W/P_N/W'] = optimum_design_option['W/P']
        dict['Performance']['CO2'] = optimum_design_option['CO2']

        value = WeightEstimation(dict).Iterations(dict['Power_prop']['bat'])[4]

        if value < minimum:
            minimum = value

        step += 1

    # printing option
    if printing:
        print(sorted_design_points[step])
        print(design_points[str(sorted_design_points[step][0])])
        print(WeightEstimation(dict).Iterations(dict['Power_prop']['bat']))

if __name__ == '__main__':

    # change this to run the iterations generator
    run = False
    if run:
        logging.info(" Starting generating the new possible design points. This may take a while.")
        Generate(p, dict, run)

    logging.info(" Getting the data from the design point options")

    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the absolute path to the data_iterations.json file
    design_json_path = os.path.join(script_dir, 'Configurations', 'data_iterations.json')

    # Print the absolute path for debugging
    logging.info(f" Looking for data_iterations.json at: {os.path.abspath(design_json_path)}")

    # Construct the absolute path to the design.json file
    data_iterations_json_path = os.path.join(script_dir,'..', 'Configurations', 'data_iterations.json')

    with open(data_iterations_json_path, 'r') as f:
        design_points = json.load(f)

    logging.info(" Opening data_iterations.json successful")

    # Weights for each key
    # increasing score: positive weight
    # decreasing score: negative weight
    # no influence: weight = 0

    weights = {
        'A': +0.05,
        'eta_p': +0.1,
        'Clmax_clean': +0.15,
        'Clmax_TO': +0,
        'Clmax_Land': +0,
        'Cd0': -0.35,
        'V_cruise': -0.1,
        'climbrate': -0.15,
        'bat': +0,
        'CO2': -0.1
    }

    maximum_weight_battery = 1000
    CO2_threshold = 50
    printing = False

    find_optimal_design(maximum_weight_battery=maximum_weight_battery, weights=weights, CO2_threshold=CO2_threshold, design_points=design_points, printing=printing)
    logging.info(f" Finding the optimal design point with a maximum battery weight of {maximum_weight_battery}[kg] with a CO2 threshold of {CO2_threshold}[%] successful")
    logging.info(" Calculating the weight components")

    # getting the weight components
    WeightEstimation = WeightEstimation(dict)
    component_weights = WeightEstimation.Iterations(dict['Power_prop']['bat'])

    # print the component weights
    print(f"Component weights: MTOW {round(component_weights[1],2)}[kg],"
          f" OEW {round(component_weights[2],2)}[kg],"
          f" Powertrain {round(component_weights[3],2)}[kg],"
          f" Battery {round(component_weights[4],2)}[kg],"
          f" Fuel {round(component_weights[5],2)}[kg],"
          f" Wing {round(component_weights[6],2)}[kg],"
          f" Wpl_des {round(component_weights[6],2)}[kg]")

    print('Total weight:', round(component_weights[1],2), '[kg] including contingency')
    print('Contingency:', (round((dict['Contingency']-1)*100,0)),"%")

    logging.info(" Calculating the weight components successful")

    logging.info(" Calculating the Xcg excursion")

    #return wcg, CGlist, xlemac
    wcg, CGlist, xlemac = iterate_cg_lg(ac_datafile = dict)

    print(f"Xcg Range is between': {round(min(CGlist), 2)} and {round(max(CGlist), 2)} [m]")

    logging.info(" Calculating the Xcg excursion successful")
    logging.info(" Calculating the MAC")

    # this line is to comply with nicholas's mood, please dont remove it as it calculates the mac
    mac = aerodynamic_design(dict, checkwingplanform=False, checkflowparameters=False, checkstability=False, checkhsplanform=False)
    print(f'MAC: {round(mac,2)} [m]')
    logging.info(" Calculating the MAC successful")
    logging.info(" Calculating the aerodynamic design")

    # initialise the checking paramaters
    check_flow_parameter = False
    check_stability = False
    check_wing_planform = False
    check_horizontal_stabilizer_planform = False

    # updating the xcg aft and xcg front
    dict['Stability']['Cg_Aft'] = (round(max(CGlist), 2) - dict['Stability']['XLEMAC_m']) / mac
    dict['Stability']['Cg_Front'] = (round(min(CGlist), 2) - dict['Stability']['XLEMAC_m']) / mac

    # run the aerodynamic design
    aerodynamic_design(dict, checkwingplanform=check_wing_planform, checkflowparameters=check_flow_parameter, checkstability=check_stability, checkhsplanform=check_horizontal_stabilizer_planform)

    logging.info(" Calculating the aerodynamic design successful")

    logging.info("Saving the modified design.json file")
    print("Finally, I am free")


    # Construct the absolute path to the design.json file
    design_json_path = os.path.join(script_dir,'..', 'Configurations', 'design.json')

    # save the modified design.json file and hope it doesnt break
    with open(design_json_path, 'w') as f:
        json.dump(dict, f, indent=4)
    logging.info("Program finished successfully")





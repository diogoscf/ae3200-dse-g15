import sys
import os
import json
import numpy as np
import logging
import colorlog

# integration v2 dont touch this please
"Dear Programmer Please do not remove this line, it is very important for the correct function of the main program"

def setup_logging():
    handler = colorlog.StreamHandler()
    handler.setFormatter(colorlog.ColoredFormatter(
        "%(log_color)s%(levelname)s:%(message)s",
        log_colors={
            'DEBUG': 'cyan',
            'INFO': 'green',
            'WARNING': 'yellow',
            'ERROR': 'red',
            'CRITICAL': 'bold_red'
        }
    ))
    logger = colorlog.getLogger()
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(script_dir,'..'))
sys.path.append(project_root)

from HumanAir.LoadingDiagram.Parameters import Parameters_ConvNoCanard as p
from HumanAir.Class_I_Weight.Class_I_Weight_Estimation import WeightEstm as WeightEstimation
from HumanAir.LoadingDiagram.Main import WP_WS
from HumanAir.CO2_Calculator.conceptual_co2 import calculate_co2_reduction_average_flight as co2
from HumanAir.Weights_and_CG.weight_fractions import find_lg, iterate_cg_lg
from HumanAir.AerodynamicDesign.Aerodynamics_Main import aerodynamic_design
from HumanAir.FinancialAnalysis.conceptual_financial_analysis import hourly_operating_cost
from HumanAir.Class_II_Weight.Class_II_Weight import RunClassII

def load_json_file(file_name):

    "Getting the design.json file"
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the absolute path to the design.json file
    design_json_path = os.path.join(script_dir,'..', 'HumanAir','Configurations', file_name)

    # Print the absolute path for debugging
    logging.info(f" Looking for {file_name} at: {os.path.abspath(design_json_path)}")

    # Attempt to open the file
    with open(design_json_path, 'r') as f:
        dict = json.load(f)

    logging.info(f" Opening {file_name} successful")

    return dict

"Generating the design points"
def Generate(p, dict, run=False):

    # tune the parameters with a reasonable range
    A_lst = np.arange(7.0, 18.51, 0.5)
    eta_p_lst = np.arange(0.8, 0.851, 0.05)
    Clmax_clean_lst = np.arange(1.6, 2.21, 0.2)
    Clmax_TO_lst = np.arange(2, 2.61, 0.2)
    Clmax_Land_lst = np.arange(2, 2.61, 0.2)
    Cd0_lst = np.arange(0.028, 0.0321, 0.002)
    V_cruise_lst = np.arange(60, 65.1, 1)
    climbrate_lst = np.arange(2.5, 5.01, 0.5)

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
                                    dict["Performance"]["endurance"]=1111200/V_cruise/3600

                                    for climbrate in climbrate_lst:
                                        current_iteration+=1

                                        # print the iteration number every 500 steps
                                        if current_iteration%500==0: 
                                            logging.info(" Iteration: "+str(current_iteration)+"/"+str(total_iterations))
                                            #print()

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
                                        dict['Performance']['W/P_N/W'], dict['Performance']['W/S_N/m2'] = WP_WS().calculate_optimal_point()

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
                                            dict['Power_prop']['P_req_cruise_W'] = dict['Performance']['P_cruise/P_TO'] * np.exp(coeff_exp[1]) * np.exp(coeff_exp[0] * bat[step])
                                            dict['Power_prop']['E_bat_Wh'] = dict['Power_prop']['P_req_cruise_W'] * dict['Performance']['endurance'] / dict['Performance']['P_cruise/P_TO'] * bat[step]

                                            # calculate the co2 ratio for the specific combination of parameters
                                            co2_ratio = co2(ac_data=dict)

                                            if co2_ratio * 100 > 25 and ok == 0:

                                                idx += 1
                                                ok = 1



                                            if co2_ratio * 100 > co2_ratio_max and co2_ratio * 100 >25 and dict['Power_prop']['E_bat_Wh']<300000 and  dict['Power_prop']['P_req_cruise_W']/0.8<250000: # the <250000 condition is for the battery to be able to be charged
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
                                                dict_iterations[str(idx)]['W/P'] = dict['Performance']['W/P_N/W']
                                                dict_iterations[str(idx)]['W/S'] = dict['Performance']['W/S_N/m2']
                                                dict_iterations[str(idx)]['bat'] = bat[step]

                                                co2_ratio_max=co2_ratio

        # save the json file with all possible design options
        data_iterations_json_path = os.path.join(script_dir, '..', "HumanAir", 'Configurations', 'data_iterations.json')
        with open(data_iterations_json_path, 'w') as file:
            json.dump(dict_iterations, file, indent=4)

"Calculating the weighted score to find the best design point"
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
                    'bat': 0.01,
                    'CO2': 30}

    # initialing the score
    score = 0

    # calculate the score
    for key, weight in weights.items():
        score += (worst_optimal_data[key] - point_data.get(key, 0))/worst_optimal_data[key] * weight # get a normalised value for the contribution of each key
    return score

"Finding the optimal design point with the highest score and the lowest battery weight"
def find_optimal_design(maximum_weight_battery=1000, weights=None, CO2_threshold=50, design_points=None, printing=False, step=0):
    optimum_design_points = {key: value for key, value in design_points.items() if value['CO2'] > CO2_threshold}
    scores = [(key, calculate_weighted_score(value, weights)) for key, value in optimum_design_points.items()]

    sorted_design_points = sorted(scores, key=lambda x: x[1], reverse=True)

    value = np.inf
    minimum = np.inf

    while value > maximum_weight_battery or design_points[sorted_design_points[step][0]]['Cd0'] < 0.027:
        step += 1

        optimum_design_option = design_points[sorted_design_points[step][0]]

        # Updating the design.json file with the optimum design option
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

        weight_estimation = WeightEstimation(dict)  # Correct instantiation of the WeightEstimation class
        iteration_results = weight_estimation.Iterations(dict['Power_prop']['bat'])
        value = iteration_results[4]
        MTOW = iteration_results[1]
        dict['Power_prop']['P_req_TO_W'] = 9.81 * MTOW / dict['Performance']['W/P_N/W']
        dict['Power_prop']['P_req_cruise_W'] = 9.81 * 0.8 * MTOW / dict['Performance']['W/P_N/W']
        dict['Power_prop']['E_bat_Wh'] = dict['Power_prop']['P_req_cruise_W'] * dict['Performance']['endurance'] / dict['Performance']['P_cruise/P_TO'] * dict['Power_prop']['bat']

        if value < minimum:
            minimum = value

    # Printing option
    if printing:
        print(sorted_design_points[step])
        print(design_points[str(sorted_design_points[step][0])])
        print(weight_estimation.Iterations(dict['Power_prop']['bat']))


"Main Function to run the program"
if __name__ == '__main__':
    # initialise the logging
    setup_logging()

    # Integration in progress v2
    logging.info(' Starting the program')

    # initialise the design.json file
    dict = load_json_file('design.json')

    run_generate = False
    run_classI = True
    run_classII = False

    if run_generate:
        logging.info(" Starting generating the new possible design points. This may take a while.")
        Generate(p, dict, run_generate)

    if run_classI:
        logging.info(" Getting the data from the design point options")

        design_points = load_json_file('data_iterations.json')

        weights = {
            'A': +0.,
            'eta_p': +0.1,
            'Clmax_clean': +0.05,
            'Clmax_TO': +0,
            'Clmax_Land': +0,
            'Cd0': -0.45,
            'V_cruise': -0.05,
            'climbrate': -0.1,
            'bat': -0.15,
            'CO2': -0.1
        }

        maximum_weight_battery = 1000
        CO2_threshold = 25
        printing = False

        # set up that the optimal stability range is not yet set
        find_optimal_stability = False
        # initialise from which step to start to search from the design_iterations.json
        step = 0

        logging.info(" Starting the search for optimal stability range")

        while not find_optimal_stability:
            find_optimal_design(maximum_weight_battery=maximum_weight_battery, weights=weights, CO2_threshold=CO2_threshold, design_points=design_points, printing=printing, step=step)
            logging.info(f" Finding the optimal design point with a maximum battery weight of {maximum_weight_battery}[kg] with a CO2 threshold of {CO2_threshold}[%] successful")
            logging.info(" Calculating the weight components")

            weight_estimation = WeightEstimation(dict)
            component_weights = weight_estimation.Iterations(dict['Power_prop']['bat'])

            print(f"Component weights:"
                f" OEW {round(component_weights[2], 2)}[kg],"
                f" Powertrain {round(component_weights[3], 2)}[kg],"
                f" Battery {round(component_weights[4], 2)}[kg],"
                f" Fuel {round(component_weights[5], 2)}[kg],"
                f" Wing {round(component_weights[6], 2)}[kg],"
                f" Wpl_des {round(component_weights[7], 2)}[kg]")

        print('Total weight:', round(component_weights[1], 2), '[kg] including contingency')
        print('Contingency:', (round((dict['Contingency'] - 1) * 100, 0)), "%")
        dict["Weights"]["MTOW_N"] = 9.81 * round(component_weights[1], 2)
        dict["Weights"]["OEW_N"] = 9.81 * round(component_weights[2], 2)
        dict["Weights"]["Wptr_N"] = 9.81 * round(component_weights[3], 2)
        dict["Weights"]["Wbat_N"] = 9.81 * round(component_weights[4], 2)
        dict["Weights"]["Wfuel_N"] = 9.81 * round(component_weights[5], 2)
        dict["Weights"]["Ww_N"] = 9.81 * round(component_weights[6], 2)
        dict["Weights"]["W_L_N"] = 9.81 * (round(component_weights[1], 2)-round(component_weights[5], 2))
        
        logging.info(" Calculating the weight components successful")

            # set up the condition to set up the range where the cg of the wing is with report of the mac
            PERCENTAGE = np.arange(-0.1, 0.51, 0.1)
            logging.info(" Starting the search for the optimal stability range in terms of where to position the cg of the wing")

            # iterate over the percentage to find the optimal stability range
            for pct in PERCENTAGE:

                logging.info(" Calculating the Xcg excursion")

                wcg, CGlist, xlemac = iterate_cg_lg(ac_datafile=dict, PERCENTAGE=pct)

                # dont remove this line as it complies with nicholas's mood
                dict['Stability']['XLEMAC_m']=xlemac
                mac = aerodynamic_design(dict, checkwingplanform=False, checkflowparameters=False, checkstability=False, checkhsplanform=False)

                dict['Stability']['Cg_Aft'] = (round(max(CGlist), 2) - dict['Stability']['XLEMAC_m']) / mac
                dict['Stability']['Cg_Front'] = (round(min(CGlist), 2) - dict['Stability']['XLEMAC_m']) / mac

                logging.info(" Prepare to check the stability")

                # dont remove this line as it complies with nicholas's mood
                aerodynamic_design(dict, checkwingplanform=False, checkflowparameters=False, checkstability=True, checkhsplanform=False)

                print("Is stability satisfied at a X_LEMAC "+ str(round(dict['Stability']['XLEMAC_m'], 2))+ " [m]" + "|"  + "[Y/N]: ")
                answer = input()

                if answer.lower() == "y":
                    break

            if answer.lower() == "y":

                # print the range of the cg
                print(f"Xcg Range is between: {round(min(CGlist), 2)} and {round(max(CGlist), 2)} [m]")

                logging.info(" Calculating the Xcg excursion successful")
                logging.info(" Calculating the MAC")

                # dont remove this line as it complies with nicholas's mood
                mac = aerodynamic_design(dict, checkwingplanform=False, checkflowparameters=False, checkstability=False,
                                        checkhsplanform=False)
                print(f'MAC: {round(mac, 2)} [m]')

                logging.info(" Calculating the MAC successful")
                logging.info(" Calculating the aerodynamic design")

                find_optimal_stability = True

                aerodynamic_design(dict, checkwingplanform=True, checkflowparameters=False, checkstability=False, checkhsplanform=True)

                logging.info(" Calculating the aerodynamic design successful")
                logging.info(" Calculating the hourly price")

                # calculating the hourly cost
                cost = hourly_operating_cost("maf_mission_graph.csv")

                print(f"Cost: {round(cost, 2)} [US$]")

                logging.info(" Calculating the hourly price successful")
                logging.info(" Saving the modified design.json file")

                # calculating if we get the copium batteries how much co2 reduction would increase
                dict["Power_prop"]["E_bat_Wh"] = 685 / 350 * dict["Power_prop"]["E_bat_Wh"]
                print("Reduction with future expected battery technology: " + str(round(co2(ac_data=dict) * 100, 2)) + "[%]")

            design_json_path = os.path.join(script_dir, '..', "HumanAir", 'Configurations', 'design.json')
            logging.info(" Design.json saved at: " + design_json_path)

                with open(design_json_path, 'w') as f:
                    json.dump(dict, f, indent=4)
    
    if run_classII:
        logging.info(" Calculate Class II Weight Groups")
        dict=RunClassII(dict,check=True)

        design_json_path = os.path.join(script_dir, '..', "HumanAir", 'Configurations', 'design.json')
        logging.info(" Design.json saved at: " + design_json_path)

        with open(design_json_path, 'w') as f:
            json.dump(dict, f, indent=4)

        logging.info(" Program finished successfully")






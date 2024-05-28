import numpy as np
import re
import pandas as pd
# Define the path to the text file
# file_path = 'data.txt'

''' DATA_EXAMPLE'''

def import_data(file_path):    
    # Initialise lists
    column_names = []
    data_rows = []
    
    # Read the text file
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        header_found = False
        for line in lines:
            if re.match(r'^\s*alpha\s+Beta\s+CL\s+CDi\s+CDv\s+CD\s+CY\s+Cl\s+Cm\s+Cn\s+Cni\s+QInf\s+XCP\s*$', line):
                column_names = re.findall(r'\S+', line.strip())
                header_found = True
                continue

            if header_found and re.match(r'^\s*-?\d+', line):
                data_rows.append(re.findall(r'-?\d+\.\d+', line.strip()))
                
    df = pd.DataFrame(data_rows, columns=column_names)
    df = df.astype(float)
    return df


'''WING_SPAN'''

def import_data2(file_path):
    data = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        angle_file = lines[0]
        angles = re.findall(r'VLM1 -\s*-?\d+\.?\d*', angle_file)
        angles = [float(re.search(r'VLM1 -\s*(-?\d+\.?\d*)', angle).group(1)) for angle in angles]


        for line in lines[1:]:
            values = line.split()  
            for i in range(0, len(values), 2):  
                angle_index = i // 2 
                angle = angles[angle_index]
                y_positions = float(values[i]) 
                lift_distribution = float(values[i + 1]) 
                
                if angle not in data:
                    data[angle] = {"y_positions": [], "lift_distribution": []}

                data[angle]["y_positions"].append(y_positions)
                data[angle]["lift_distribution"].append(lift_distribution)
    return data
import numpy as np
import re

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
                    data[angle] = {"y_span": [], "coefficient": []}

                data[angle]["y_span"].append(y_positions)
                data[angle]["coefficient"].append(lift_distribution)
    return data

def chord(Sw, taper_ratio, Cl_DATA, AoA, n):
    b = Cl_DATA[AoA]['y_span'][-1] * 2  
    # Generate spanwise coordinate points
    y = np.linspace(Cl_DATA[AoA]['y_span'][0], Cl_DATA[AoA]['y_span'][-1], n)  # n is the number of nodes
    # Calculate the chord distribution
    chord_length = 2 * Sw / (1 + taper_ratio) / b * (1 - (1 - taper_ratio) * np.abs(2 * y / b))
    return chord_length, y
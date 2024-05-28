import numpy as np
import re
import pandas as pd
# Define the path to the text file

def import_data(file_path):
    # file_path = 'data.txt'
    
    # Initialise lists
    column_names = []
    data_rows = []
    
    # Read the text file
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        header_found = False
        for line in lines:
            if re.match(r'^\s*\S+\s\S+', line):
                column_names = re.findall(r'\S+', line.strip())
                header_found = True
                continue
            
            if header_found and re.match(r'^\s*-?\d+', line):
                data_rows.append(re.findall(r'-?\d+\.\d+', line.strip()))
    
    df = pd.DataFrame(data_rows, columns=column_names)
    df = df.astype(float)
    return df


        
    
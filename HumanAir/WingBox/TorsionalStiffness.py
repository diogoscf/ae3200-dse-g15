import numpy as np
from Functions import chord, import_data2
import pandas as pd
import matplotlib.pyplot as plt

x_pos = np.array([0.15, 0.5])

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def import_data(file_path):
    df = pd.read_csv(file_path, sep='\s+', header=None, names=['x', 'y'], skiprows=1)
    return df

def interp(x, x_given, y_given):
    i = find_nearest(x_given, x)
    if i == 0:
        return y_given[0]
    elif i == len(x_given):
        return y_given[-1]
    else:
        x1, x2 = x_given[i - 1], x_given[i]
        y1, y2 = y_given[i - 1], y_given[i]
        return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
    
def spars(chord, x_pos):
    chord = np.array(chord).reshape((len(chord),1))
    return chord * x_pos 

def airfoil_division(df):
    for i in range(len(df)):
        if df['x'][i] == 0:
            idx = i
    df_1 = np.array(df.iloc[:idx,:].values.tolist())
    df_2 = np.array(df.iloc[idx:,:].values.tolist())
    return df_1, df_2


def y_spar(df, x_spar, chord_length):
    y_spar = []
    for i in range(len(x_spar)):
        for j in range(len(x_spar[i])):
            y_spar.append(interp(x_spar[i][j]/chord_length[i], df[:,0], df[:,1])*chord_length[i])
    y_spar = np.array(y_spar)
    return y_spar.reshape((len(x_spar), 2))        

def MOI(y, chord, spars):
    return

########## input ###########
Sw = 39  # [m2]
taper_ratio = 0.4
Cl_DATA = import_data2('HumanAir\WingBox\Cl_DATA.txt')
AoA = -6
n = 3

######## Execution ########
# Assuming chord function takes the arguments as described and returns a chord length
chord_length, y = chord(Sw, taper_ratio, Cl_DATA, AoA, n)
print(chord_length)
# Calculate the spars positions based on the chord length
x_spar = spars(chord_length, x_pos)

# airfoil data
df = import_data('HumanAir/WingBox/airfoil.txt')
df_up, df_down = airfoil_division(df)
y_spar_up = y_spar(df_up, x_spar, chord_length)
y_spar_down = y_spar(df_down, x_spar, chord_length)

print(y_spar_down)
print(y_spar_up)

interp(0.15, df_up[:,0], df_up[:,1])
#print(y_spar_up)
x_pos_flat = np.array(x_spar.flatten())

'''
plt.plot(df_down[:,0], df_down[:,1])
plt.plot(df_up[:,0], df_up[:,1])
plt.plot(x_spar[0]/chord_length[0], y_spar_down[0]/chord_length[0])
plt.plot(x_spar[0]/chord_length[0], y_spar_up[0]/chord_length[0])
plt.show()
'''

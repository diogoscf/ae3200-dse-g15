import matplotlib.pyplot as plt
import numpy as np

def Plotx(xvalue, ylst, name: str):
    xlst=np.ones(len(ylst))*xvalue
    plt.plot(xlst, ylst, label=name)

def Ploty(xlst, ylst, name: str):
    plt.plot(xlst, ylst, label=name)
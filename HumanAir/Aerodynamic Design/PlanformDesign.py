import numpy as np
import matplotlib.pyplot as plt

class Planform:
    def __init__(self, MTOW, WS, AR, Taper, QuarterChordSweep):
        self.MTOW = MTOW
        self.WS = WS
        self.AR = AR
        self.Taper = Taper
        self.QuarterChordSweep = QuarterChordSweep
        self.g = 9.81

    def WingSurfaceArea(self):
        return self.MTOW*self.g/self.WS

    def WingSpan(self):
        return np.sqrt(self.WingSurfaceArea()*self.AR)

    def RootChord(self):
        return 2*self.WingSurfaceArea()/((1+self.Taper)*self.WingSpan())

    def TipChord(self):
        return self.Taper*self.RootChord()

    def MAC(self):
        return (2*self.RootChord()/3)*(1+self.Taper+self.Taper**2)/(1+self.Taper)

    def MAC_y(self):
        return (self.WingSpan()/6)*(1+2*self.Taper)/(1+self.Taper)


    def PlotWingPlanform(self):
        plt.figure()
        plt.plot([0,0], [self.RootChord(),0],color='black')                                                                             #Plotting Root Chord
        tip_x_qc=self.RootChord()-self.RootChord()/4-np.tan(np.deg2rad(self.QuarterChordSweep))*self.WingSpan()/2
        plt.plot([self.WingSpan()/2, self.WingSpan()/2], [tip_x_qc+self.TipChord()/4, tip_x_qc-3*self.TipChord()/4],color='black')      #Plotting Tip Chord
        plt.plot([0, self.WingSpan()/2], [self.RootChord(),tip_x_qc+self.TipChord()/4],color='black')                                   #Plot Leading Edge
        plt.plot([0,self.WingSpan()/2], [0, tip_x_qc-3*self.TipChord()/4], color='black')                                               #Plot Trailing Edge
        plt.plot([0,self.WingSpan()/2], [self.RootChord()-self.RootChord()/4, tip_x_qc])                                                #Plot Quarter Chord
        plt.axis('equal')
        plt.show()

if __name__ == "__main__":
    WingConventional=Planform(1946.34, 618, 9.35, 0.4, 15)


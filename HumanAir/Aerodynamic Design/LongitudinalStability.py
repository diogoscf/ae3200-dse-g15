import numpy as np
import matplotlib.pyplot as plt
class LongitudinalStability:
    def __init__(self, CLh, CLah, Xcgh, XLEMAC, VhV, Wing, CMac, Xac, SM, deda, FuselageLength, CgAft, CgFwd, airfoil, Clalphah = 5.93792624, Clalphaah = 6.237426906):
        self.CLh = CLh
        self.CLah = CLah
        self.FuselageLength = FuselageLength
        self.lh=Xcgh*FuselageLength-(XLEMAC+0.4*Wing.MAC())
        self.VhV=VhV
        self.c=Wing.MAC()
        self.CMac=CMac
        self.Xac=Xac
        self.dxcg=0.01
        self.Xcg=np.arange(0,1+self.dxcg,self.dxcg)
        self.SM=SM
        self.Clalphah=Clalphah
        self.Clalphaah=Clalphaah
        self.deda=deda
        self.CgFwd=CgFwd
        self.CgAft=CgAft
        self.Airfoil=airfoil
        print("Longitudinal Stability Initialized")

    def CMac_Wing(self):
        return self.Airfoil['Cm_0']*()

    def Stability(self):
        return (self.Xcg-self.Xac+self.SM)/(self.Clalphah/self.Clalphaah*(1-self.deda)*self.lh/self.c*self.VhV**2)

    def Stability_NoMargin(self):
        return (self.Xcg - self.Xac) / (self.Clalphah / self.Clalphaah * (1 - self.deda) * self.lh / self.c * self.VhV ** 2)

    def Controllability(self):
        return (self.Xcg+self.CMac/self.CLah-self.Xac)/(self.CLh/self.CLah*self.lh/self.c*self.VhV**2)

    def ShS(self):
        Cgmin=self.CgFwd
        Cgmax=self.CgAft

        ShSStability=np.interp(Cgmax, self.Xcg, self.Stability())
        ShSControllability=np.interp(Cgmin, self.Xcg, self.Controllability())

        return np.max([ShSStability, ShSControllability])

    def Plotting(self):
        plt.figure()
        plt.plot(self.Xcg, self.Stability(), label="Stability", color='green', linestyle='solid')
        plt.plot(self.Xcg, self.Stability_NoMargin(), label="Stability (No Margin)", color='green', linestyle='dashed')
        plt.plot(self.Xcg, self.Controllability(), label="Controllability", color='blue', linestyle='solid')
        plt.plot([self.CgFwd, self.CgAft], [self.ShS(), self.ShS()], label="Cg Excursion", color='red', linestyle='solid')
        plt.legend()
        plt.show()





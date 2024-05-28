import numpy as np
import matplotlib.pyplot as plt
import json

class WeightEstm:
    def __init__(self,dict):

        self.dict = dict

    def OEW_prime (self):
        return self.dict["A"]*self.dict["MTOW"]+self.dict["B"]
    
    def PowertrainWeight(self,bat):
        return 9.81*self.dict["MTOW"]/1000/self.dict["W/P"]/self.dict["eta_gb"]/self.dict["eta_ptr"]*(1-bat)/self.dict["ptr_specific_power"]
    
    def BatteryWeight(self,bat):
        return 9.81*self.dict["MTOW"]/self.dict["W/P"]*self.dict["endurance"]*bat/self.dict["bat_specific_energy"]/self.dict["eta_bat"]/self.dict["DoD"]/self.dict["eta_EM"]
    
    def FuelWeight(self,bat):
        return 9.81*self.dict["MTOW"]/self.dict["W/P"]*(1-bat)*self.dict["endurance"]/self.dict["fuel_specific_energy"]/self.dict["eta_thermal_engine"]
    
    def WingWeight(self):
        return self.dict["Aw"]*self.dict["MTOW"]+self.dict["Bw"]
    
    #returs MTOW w/o cont, MTOW w cont, OEW w cont, Wptr w cont, Wbat w cont, Wf w fuel, Payload w cont
    def Iterations(self,bat):
        MTOW_new=0
        MTOW_old=self.dict["MTOW"]
        ok=False

        while np.abs((MTOW_new-self.dict["MTOW"])/self.dict["MTOW"])>0.02:

            if ok: 
                self.dict["MTOW"]=MTOW_new
            
            OEW_prime=self.OEW_prime()
            PowertrainWeight=self.PowertrainWeight(bat)
            BatteryWeight=self.BatteryWeight(bat)
            FuelWeight=self.FuelWeight(bat)
            WingWeight=self.WingWeight()

            MTOW_new=OEW_prime+PowertrainWeight+BatteryWeight+FuelWeight+WingWeight+self.dict["Payload"]

            ok=True
        
        self.dict["MTOW"]=MTOW_old
        return(MTOW_new,self.dict["cont"]*MTOW_new,self.dict["cont"]*OEW_prime,self.dict["cont"]*PowertrainWeight,self.dict["cont"]*BatteryWeight,self.dict["cont"]*FuelWeight,self.dict["cont"]*WingWeight,self.dict["cont"]*self.dict["Payload"])

    def PolynomialRegression(self):

        bat=np.arange(0,0.18,0.001)
        lst_P=[]

        for pbat in bat:
            row=self.Iterations(pbat)
            lst_P.append(9.81*row[1]/self.dict["W/P"])
    

        coeff_exp=np.polyfit(bat,np.log(lst_P),1)
        coeff_pol=np.polyfit(bat,lst_P,2)
        y_pol=coeff_pol[0]*bat**2+coeff_pol[1]*bat+coeff_pol[2]
        y_exp=np.exp(coeff_exp[1])*np.exp(coeff_exp[0]*bat)
        #plt.scatter(bat,lst_P)
        #plt.plot(bat,y_exp)
        #plt.plot(bat,y_pol)
        #plt.show()


if __name__ == "__main__":
    file_path='HumanAir/Configurations/conventional - Nicholas.json'

    data=WeightEstm(file_path)

    bat=0.11
    row=data.Iterations(bat)
    print(row)
    data.PolynomialRegression()


            
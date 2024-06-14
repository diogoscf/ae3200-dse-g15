from scipy.integrate import quad
import numpy as np
import aircraft
from helper import density
import matplotlib.pyplot as plt

def integrate(f, x0, x1):
    y, abserr = quad(f, x0, x1)
    if abserr > 0.000001:
        raise Exception("Integration error too large")
    return y

g = 9.80665

class MAircraft:
    """
    Mission-Aircraft class.
    Wraps Aircraft object for use for mission analysis.
    """
    
    def __init__(self, acf=aircraft.Aircraft()):
        self.acf = acf
        
        #
        # constants
        #
        
        self.dT = self.acf.dT_default # ISA temperature offset used throuhout flight [deg C]
        self.S  = self.acf.S() # clean config
        
        
        #
        # variables for mission analysis
        #
        

        self.V = 0 # current true air speed [m/s]
        self.h = 0 # current geopotential altitude [m]
        
        self.payload = self.acf.W_pl_no_pilot # current payload [N]
        self.fuel = self.acf.W_MF # current fuel weight [N]
        self.bat_cap = self.acf.max_bat_cap # current battery capacity [Wh]
        
        self.ground_distance = 0 # distance flown [m]
        self.flight_time = 0 # time since start takeoff [s]
    
    
    def _get_W(self):
        """ Returns current gross weight [N] """
        return self.acf.W_OE + self.payload + self.fuel
    W = property(fget=_get_W)

    def _get_rho(self):
        """ Returns air density in current conditions """
        return density(self.h, self.dT)
    rho = property(fget=_get_rho)
    
    def _get_CL(self):
        """ Returns current CL assuming small flight path angle (L = W) """
        return self.W / (0.5 * self.rho * self.V**2 * self.S)
    CL = property(fget=_get_CL)
    
    def _get_D(self):
        """ Returns current drag force for clean config assuming small flight path angle """
        return 0.5 * self.rho * self.V**2 * self.S * self.acf.CD(self.CL)
    D = property(fget=_get_D)
    
    def P_r(self):
        """ Power required for steady level flight """
        return self.D * self.V
    
    def _fly_accelerate_const_alt(self, V_target, electric=False, time_step=.5):
        """ Accelerate with max const power """
        if self.V > V_target:
            raise Exception("This function cannot decelerate")
        elif np.isclose(self.V, V_target):
            return
        
        t_lst = []
        v_lst = []
        s_lst = []
        
        while self.V < V_target:
            T = self.acf.T(self.V, self.h, self.dT, electric=electric)
            a = g * (T - self.D) / self.W # const acceleration
            
            if electric:
                self.bat_cap -= self.acf.bat_cap_rate(P_shaft=self.acf.P_shaft(self.h, self.dT, electric=True)) * time_step
            else:
                self.fuel -= self.acf.fuel_rate(P_shaft=self.acf.P_shaft(self.h, self.dT, electric=False)) * time_step

            self.ground_distance += (self.V + a*time_step/2) * time_step # avg speed over dt
            self.V += a * time_step
            self.flight_time += time_step
            
            t_lst.append(self.flight_time)
            v_lst.append(self.V)
            s_lst.append(self.ground_distance)
            
        # plt.figure()
        # plt.plot(t_lst, v_lst, 'g')
        # plt.show()
        
        # plt.figure()
        # plt.plot(t_lst, s_lst, 'b')
        # plt.show()
            
        
        
    def _fly_const_CL_const_climbrate(self, RC, h_target, electric=False, time_step=30):
        """ Climb at constant CL (current CL) and constant climb rate """
        if self.h > h_target:
            raise Exception("This function cannot descent")
        elif np.isclose(self.h, h_target):
            return
        
        t_lst = []
        v_lst = []
        s_lst = []
        h_lst = []
        
        CL = self.CL # will be kept constant
        
        while self.h < h_target:
            # with small angle approximation P = RC*W + P_r
            # and P_r = V * (D + W/g*a)
            
            # the target speed we need to accelerate or decelerate to the next instance
            rho2 = density(self.h+RC*time_step, self.dT)
            V2 = np.sqrt(self.W / (0.5 * rho2 * self.S * CL))
            a = (V2 - self.V) / time_step
            P = RC * self.W + self.V * (self.D + self.W/g*a)
            P_shaft = P / self.acf.prop_eff(self.V, self.h, self.dT)
            
            if electric:
                self.bat_cap -= self.acf.bat_cap_rate(P_shaft=P_shaft) * time_step
            else:
                self.fuel -= self.acf.fuel_rate(P_shaft=P_shaft) * time_step

            V_avg = self.V + (V2 - self.V)/2 # avg speed over dt assuming const a
            V_g   = np.sqrt(V_avg**2 - RC**2) # ground speed (so we are taking into account flight path angle here)
            self.ground_distance += V_g * time_step
            self.V += a * time_step
            self.flight_time += time_step
            self.h += RC * time_step
            
            t_lst.append(self.flight_time)
            v_lst.append(self.V)
            s_lst.append(self.ground_distance)
            h_lst.append(self.h)
            
        plt.figure()
        plt.plot(t_lst, v_lst, 'g')
        plt.show()
        
        plt.figure()
        plt.plot(t_lst, s_lst, 'b')
        plt.show()
        
        plt.figure()
        plt.plot(t_lst, h_lst, 'r')
        plt.show()
            
    def _fly_const_V_const_alt(self, target_distance, time_step=60, electric=False):
        """ Fly keeping current altitude and speed constant """
        
        w_lst = []
        s_lst = []
        
        distance = 0
        
        prop_eff = self.acf.prop_eff(self.V, self.h, self.dT)

        while distance < target_distance:
            P_shaft = self.V * self.D / prop_eff
            
            if electric:
                self.bat_cap -= self.acf.bat_cap_rate(P_shaft=P_shaft) * time_step
            else:
                self.fuel -= self.acf.fuel_rate(P_shaft=P_shaft) * time_step

            distance += self.V * time_step
            self.flight_time += time_step
            
            s_lst.append(distance)
            w_lst.append(self.W)
        
        self.ground_distance += distance
            
        plt.figure()
        plt.plot(s_lst, w_lst, 'g')
        plt.show()
        
    def _fly_const_CL_const_alt(self, CL, time):
        foo = 1 # TODO: hier verder

    def take_off(self, elevation, surface, electric=False):
        # weight not variable here
            
        dt, self.V, _ = self.acf.takeoff_ground_run(self.W(), elevation, self.dT, 0, surface, electric=electric, calc_time=True)
        
        if electric:
            self.bat_cap -= self.acf.bat_cap_rate(TO=True) * dt
        else:
            self.fuel -= self.acf.fuel_rate(TO=True) * dt
            
        print(dt)
            
        self.flight_time += dt
        self.h = elevation # neglect what little alt has been gained
        
        
    def climb(self, h_target, rate, electric=False):
        V_start_climb = self.acf.V_Prmin(self.W, self.h, self.dT)
        
        # V_start_climb = self.acf.V_Prmin(self.W(), self.h, self.dT)
        # if self.V < V_start_climb:
        #     # accelerate to optimal climb speed at full max cont power
        #     # delta t = integral of dV/a
        #     rho, S = self.rho(), self.acf.S()
        #     def distance(V):
        #         CL = self.W() / (0.5 * rho * V**2 * S)
        #         D = 0.5 * rho * V**2 * S * self.acf.CD(CL)
        #         T = self.acf.T(V, self.h, self.dT, electric=electric)
        #         a = 9.80665 * (T-D) / self.W()
        #         return V/a
        #     ds = integrate(distance, self.V, V_start_climb)
            
        #     def time(V):
        #         CL = self.W() / (0.5 * rho * V**2 * S)
        #         D = 0.5 * rho * V**2 * S * self.acf.CD(CL)
        #         T = self.acf.T(V, self.h, self.dT, electric=electric)
        #         a = 9.80665 * (T-D) / self.W()
        #         return 1/a
        #     dt = integrate(time, self.V, V_start_climb)

        #     print(ds)
        #     print(dt)
        #     print(self.acf.fuel_rate(P_shaft=self.acf.P_shaft(self.h, self.dT, electric=False)) * dt)
        #     if electric:
        #         self.bat_cap -= self.acf.bat_cap_rate(P_shaft=self.acf.P_shaft(self.h, self.dT, electric=True)) * dt
        #     else:
        #         self.fuel -= self.acf.fuel_rate(P_shaft=self.acf.P_shaft(self.h, self.dT, electric=False)) * dt
            
        #     self.V = V_start_climb
        #     self.time_flown += dt
            
        # elif not np.isclose(self.V, V_start_climb):
        #     raise Exception(f"Climb function cannot decelerate aircraft - V1 = {self.V:.2f} m/s V2 = {V_start_climb:.2f}")
        
        
        
    def cruise(acf, ground_distance, h=None, dT=None):
        # TODO: check if flying at max CL/CD is possible
        foo = 1
        
    def loiter(acf, t, h=None, dT=None):
        foo = 2
        
    def descent(acf, h=None, dT=None, rate=2.5):
        ground_distance = 0
        return ground_distance
        
    def land(acf, h=None, dT=None, surface="grass"):
        foo =1


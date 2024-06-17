import numpy as np
import aircraft
from helper import density
import matplotlib.pyplot as plt

g = 9.80665


class MAircraft:
    """
    Mission-Aircraft class.
    
    Wraps Aircraft object for use for mission analysis.
    """
    
    def __init__(self, elevation, acf=aircraft.Aircraft()):
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
        self.h = elevation # current geopotential altitude [m]
        
        self.payload = self.acf.W_pl_no_pilot # current payload [N]
        self.fuel = self.acf.W_MF # current fuel weight [N]
        self.bat_cap = self.acf.max_bat_cap # current battery capacity [Wh]
                
        self.ground_distance = 0 # distance flown [m]
        self.flight_time = 0 # time since start takeoff [s]
        
        
        #
        # lists to keep track of changing parameters
        #
        
        self.t_lst  = [0] # time since start takeoff [s]
        self.v_lst  = [self.V] # speed [m/s]
        self.h_lst  = [self.h] # altitude [m]
        self.s_lst  = [0] # ground distance [m]
        self.f_lst  = [self.fuel] # fuel weight [N]
        self.b_lst  = [self.bat_cap] # battery capacity [Wh]
        self.CL_lst = [0] # lift coefficient [-]
        
        
        #
        # variables for printing energy consumption
        #
        
        self.last_fuel = self.fuel # to calculate fuel used per flight part
        self.last_bat_cap = self.bat_cap # to calculate fuel ued per flight part

    
    def _update_lists(self):
        """
        Updates flight log that will be plotted with plot_flight().

        Returns
        -------
        None.

        """
        # print warnings if neccessary
        if self.b_lst[-1] > 0 and self.bat_cap < 0:
            print("***** WARNING: battery capacity below 0 ******")
        if self.f_lst[-1] > 0 and self.fuel < 0:
            print("***** WARNING: fuel level below 0 ******")
            
        # append current values to lists
        self.t_lst.append(self.flight_time)
        self.v_lst.append(self.V)
        self.h_lst.append(self.h)
        self.s_lst.append(self.ground_distance)
        self.f_lst.append(self.fuel)
        self.b_lst.append(self.bat_cap)
        if self.V > 0:
            self.CL_lst.append(self.CL)
        else:
            self.CL_lst.append(0)
            
        
    def plot_flight(self):
        """
        Plots the flight logs.

        Returns
        -------
        None.

        """
        fig, axs = plt.subplots(6,1)
        
        fig.set_figheight(13)
        fig.set_figwidth(10)
        
        axs[0].plot(self.t_lst, self.h_lst)
        
        axs[1].plot(self.t_lst, self.v_lst)
        
        axs[2].plot(self.t_lst, self.f_lst)
        
        axs[3].plot(self.t_lst, self.b_lst)
        
        axs[4].plot(self.t_lst, self.s_lst)
        
        axs[5].plot(self.t_lst[1:-1], self.CL_lst[1:-1])
        
        plt.show()
        
        
    def _print_update(self, label):
        """
        Prints fuel used for last flight segment.

        Parameters
        ----------
        label : string
            A description of the flight segment to be printed with info.

        Returns
        -------
        None.

        """
        # get and print difference
        fuel_used = self.last_fuel - self.fuel
        bat_cap_used = self.last_bat_cap - self.bat_cap
        print(f"** {fuel_used:8.1f} N fuel and {bat_cap_used:8.1f} Wh energy used during {label}")
        
        # update
        self.last_fuel = self.fuel
        self.last_bat_cap = self.bat_cap
        
    def print_energy_level(self):
        """
        Prints how much fuel is used and how much is left over.

        Returns
        -------
        None.

        """
        f_1 = self.f_lst[0]
        b_1 = self.b_lst[0]
        print(f"{self.fuel:8.1f} N fuel left over ({self.fuel/f_1*100:.1f}%)")
        print(f"{self.bat_cap:8.1f} Wh battery left over ({self.bat_cap/b_1*100:.1f}%)")

        
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
    
    
    def fly_accelerate_const_alt(self, V_target, electric=False, time_step=.5):
        """
        This function will either accelerate or decelerate the aircraft.
        depending on whether V_target is greater or smaller than the current
        velocity.
        
        Acceleration is performed with max const. throttle.
        
        Deceleration is performed with zero throttle.

        Parameters
        ----------
        V_target : float
            True air speed at end of de/acceleration.
        electric : boolean, optional
            Whether this is performed with the electric motor. The default is
            False.
        time_step : float, optional
            Time step for integration in seconds. The default is .5.

        Returns
        -------
        None.

        """
        # check if we need to change speed at all
        if np.isclose(self.V, V_target):
            return
        
        # check whether we need to accelerate or decelerate
        if self.V < V_target:
            accelerate = True
        else:
            accelerate = False
        
        #
        # integration
        #
        while (accelerate and self.V < V_target) or (not accelerate and self.V > V_target):
            if accelerate:
                # max constant thrust
                T = self.acf.T(self.V, self.h, self.dT, electric=electric)
                P_shaft = self.acf.P_shaft(self.h, self.dT, electric=electric)
                if electric:
                    # update battery capacity
                    self.bat_cap -= self.acf.bat_cap_rate(P_shaft) * time_step
                else:
                    # update fuel level
                    self.fuel -= self.acf.fuel_rate(P_shaft) * time_step
            else:
                # zero throttle
                T = 0
                
            # calculate deceleration or acceleration, assumed constant during time_step
            a = g * (T - self.D) / self.W

            # update aircraft parameters
            self.ground_distance += (self.V + a*time_step/2) * time_step # avg speed over dt
            self.V               += a * time_step
            self.flight_time     += time_step
            self._update_lists()
            
        # done, print energy usage
        self._print_update("de/acceleration")
        
        
    def fly_const_CL_const_climbrate(self, RC, h_target, electric=False, time_step=30):
        """
        Make aircraft climb at constant climb rate while maintaining current
        lift coefficient.

        Parameters
        ----------
        RC : float
            Climb rate [m/s], may not be negative.
        h_target : float
            Target altitude [m].
        electric : boolean, optional
            Whether electric power is used. The default is False.
        time_step : float, optional
            Time step for integration [s]. The default is 30.

        Raises
        ------
        Exception
            Raises exception if h_target is below current h.

        Returns
        -------
        None.
        
        """
        # check if we are indeed having to climb
        if h_target < self.h:
            raise Exception("This function cannot be used for descent")
        elif np.isclose(self.h, h_target):
            return # already at target altitude
        
        # get current CL, this will be kept constant
        CL = self.CL
        
        #
        # integrate
        #
        while self.h < h_target:
            # find the acceleration required to maintain current CL
            rho2 = density(self.h+RC*time_step, self.dT) # density at next integration step
            V2 = np.sqrt(self.W / (0.5 * rho2 * self.S * CL)) # speed at next integration step
            a = (V2 - self.V) / time_step
            
            # calculate shaft power required
            # with small angle approximation P = RC*W + P_r
            # and P_r = V * (D + W/g*a)
            P = RC * self.W + self.V * (self.D + self.W/g*a)
            P_shaft = P / self.acf.prop_eff(self.V, self.h, self.dT)
            
            if electric:
                # update battery capacity
                self.bat_cap -= self.acf.bat_cap_rate(P_shaft) * time_step
            else:
                # update fuel level
                self.fuel -= self.acf.fuel_rate(P_shaft) * time_step

            # determine ground speed
            V_avg = self.V + (V2 - self.V)/2 # avg speed over dt assuming const a
            V_g   = np.sqrt(V_avg**2 - RC**2) # ground speed (so we are taking into account flight path angle here)
            
            # update aircraft parameters
            self.ground_distance += V_g * time_step
            self.V += a * time_step
            self.flight_time += time_step
            self.h += RC * time_step
            self._update_lists()
            
        # done, print energy consumed
        self._print_update("climb")
            
        
    def fly_const_V_descent(self, h_target, time_step=30):
        """
        Descent at variable descent rate to maintain current speed.

        Parameters
        ----------
        h_target : float
            Target altitude [m].
        time_step : float, optional
            Time step for integration in seconds. The default is 30.

        Raises
        ------
        Exception
            Raised if h_target > current h or is descent rate exceeds 3m/s.

        Returns
        -------
        None.

        """
        # check if we are going to descent
        if h_target > self.h:
            raise Exception("This function cannot be used for ascent")
        elif np.isclose(self.h, h_target):
            return # already at target altitude
        
        # store current speed so that we can keep it constant
        V = self.V
        
        #
        # integration
        #
        while self.h > h_target:
            # with small angle approximation RC = (P_a - P_r)/W
            # thus with P_a = 0 -> RC = -P_r/W
            RC = - V*self.D / self.W
            
            # TODO: uncomment
            # check if descent rate is reasonable
            #if RC < -3:
            #    raise Exception(f"Descent RoC unreasonably large: {RC:.2f} m/s")
            #print(f"RoC descent: {RC:.2f} m/s")
            
            # determine ground speed
            V_g   = np.sqrt(V**2 - RC**2) # we are taking into account flight path angle here
            
            # update aircraft parameters
            self.ground_distance += V_g * time_step
            self.flight_time += time_step
            self.h += RC * time_step
            self._update_lists()
            
        # print energy consumed
        self._print_update("descent")

            
    def fly_const_V_const_alt(self, target_distance, time_step=60, electric=False):
        """
        Fly keeping current altitude and speed constant (good for cruise). When
        electric is set to true it automatically switches to the fuel engine 
        when the battery is depleted (reaches max Depth of Discharge).

        Parameters
        ----------
        target_distance : float
            Ground distance to cover [m].
        time_step : float, optional
            Time step for integratin [s]. The default is 60.
        electric : boolean, optional
            Whether to use the electric motor. Once the battery is depleted
            it automatically switches to the fuel engine. The default is False.

        Returns
        -------
        None.

        """
        print(f"CL before: {self.CL:.2f} CL/CD: {self.CL/self.acf.CD(self.CL):.2f} opt: {self.acf.LDmax:.2f}")
        # stores distance covered during this phase        
        distance = 0
        
        # constant V means constant prop efficiency
        prop_eff = self.acf.prop_eff(self.V, self.h, self.dT)
        
        # constant alt = constant rho
        rho = self.rho
        
        #
        # integration
        #
        while distance < target_distance:
            # get required shaft power
            P_shaft = self.V * 0.5 * rho * self.V**2 * self.S * self.CL/12 / prop_eff
            
            if electric:
                # make sure DoD is not exceeded
                energy = self.acf.bat_cap_rate(P_shaft) * time_step
                DoD = 1 - (self.bat_cap-energy) / self.acf.max_bat_cap
                
                if DoD < self.acf.max_DoD_bat:
                    # update battery capacity
                    self.bat_cap -= self.acf.bat_cap_rate(P_shaft) * time_step
                else:
                    # switch to fuel
                    electric = False
                    self.fuel -= self.acf.fuel_rate(P_shaft) * time_step
            else:
                # update fuel level
                self.fuel -= self.acf.fuel_rate(P_shaft) * time_step

            # update ground distance of this segment
            distance += self.V * time_step
            
            # update aircraft parameters
            self.ground_distance += self.V * time_step
            self.flight_time += time_step
            self._update_lists()
            
        # done, print energy consumption
        self._print_update("cruise")
        print(f"CL after: {self.CL:.2f}")

        
    def fly_const_CL_const_alt(self, duration, time_step=60, electric=False):
        """
        Loiter at constant CL at constant altitude without covering more ground
        distance. The current CL is maintained.

        Parameters
        ----------
        duration : float
            Duration of this flight segment [s].
        time_step : float, optional
            Time step for integration. The default is 60.
        electric : boolean, optional
            Whether electric power is used. The default is False.

        Returns
        -------
        None.

        """
        # const alt = const rho        
        rho = self.rho
        
        # store CL to keep constant
        CL = self.CL

        # duration of current phase        
        t = 0
        
        #
        # integrate
        #
        while t < duration:
            # get required shaft power, neglecting deceleration to maintain CL
            P_shaft = self.V * self.D / self.acf.prop_eff(self.V, self.h, self.dT)
            
            if electric:
                # update battery capacity
                self.bat_cap -= self.acf.bat_cap_rate(P_shaft) * time_step
            else:
                # update fuel level
                self.fuel -= self.acf.fuel_rate(P_shaft) * time_step

            # update speed to maintain CL
            self.V = np.sqrt(self.W / (0.5 * rho * self.S * CL))
            
            # update time
            t += time_step
            
            # update aircraft parameters
            # no ground distance covered
            self.flight_time += time_step
            self._update_lists()
            
        # done, print energy usage
        self._print_update("loiter")
        
    def take_off(self, elevation, surface, electric=False):
        """
        Make the aircraft take off. 
        
        At the end the aircraft parameters are:
            h = airfield elevation
            V = liftoff speed
            t = takeoff duration
        
        Weight and battery capacity are also updated but fuel consumption is
        not taken into account during takeoff run (tiny effects anyway).

        Parameters
        ----------
        elevation : float
            Airfield elevation [m].
        surface : string
            Runway surface ("grass" or "paved).
        electric : boolean, optional
            Whether to use electric power. The default is False.
            
        Raises
        ------
        Exception
            Raised if gross weight exceeds MTOW.

        Returns
        -------
        None.

        """         
        # check if MTOW is not exceeded
        if self.W > self.acf.W_MTO:
            raise Exception(f"Aircraft weight ({self.W:.0f} N) exceeds MTOW ({self.acf.W_MTO:.0f} N)")
            
        # get duration and final speed of takeoff phase
        t, V, _ = self.acf.takeoff_ground_run(self.W, elevation, self.dT, 0, surface, electric=electric, calc_time=True)
        
        # get max TO shaft power
        P_shaft = self.acf.P_shaft(elevation, self.dT, use_takeoff_power=True, electric=electric)
        
        
        if electric:
            # update battery capacity
            self.bat_cap -= self.acf.bat_cap_rate(P_shaft) * t
        else:
            # update fuel level
            self.fuel -= self.acf.fuel_rate(P_shaft) * t
                        
        # update aircraft parameters
        # note that there may be multiple takeoffs in a flight so we cannot
        # reset these parameters
        self.flight_time += t
        self.V = V
        self.h = elevation # neglect what little alt has been gained
        self._update_lists()
        
        # print energy used
        self._print_update("take off")
    
        
    def land(self):
        """ Lands aircraft - only sets V to zero """
        self.V = 0
        self._update_lists()

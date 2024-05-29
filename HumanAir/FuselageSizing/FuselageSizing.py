import math
import numpy as np

class FuselageSizing:
    '''
    Seat specification
    '''
    w_pax = 0.56  # seat width [m]
    l_pax = 1.09  # seat length [m]
    h_pax = 1.40  # seat height [m]
    s = 0.01  # clearance [m]
    n_row_seat = 2  # seats per row [-]
    
    '''
    Cockpit specification
    '''
    l_cock = 1.42  # cockpit length [m]
    h_cock = 1.2  # cockpit height [m]
    h_cockeye = 1.1  # cockpit eye height [m]
    
    '''
    Battery specification
    '''
    s_battery = 0.1  # battery clearance [m]
    
    '''
    Landing gear specification
    '''
    s_aft = 0.05  # landing gear clearance [m]
    
    '''
    Other specifications
    '''
    l_enbu = 0.15  # engine/firewall buffer [m]
    h_floor = 0.1  # floor height [m]
    l_extra = 0.2  # extra space [m]
    t_fuse = 0.04  # fuselage structural thickness [m]
    s_bat_wheel = 0.05  # margin for battery and wheels
    
    def __init__(self, n_seat, w_engine, l_engine, h_engine, D_nose, h_nose, D_main, h_main, h_nose_strut, h_main_strut, l_main_lateral, l_long_nose, l_long_main, l_tailcone, h_tail, V_battery): 
        self.n_seat = n_seat  # number of total seats
        self.w_engine = w_engine  # width of engine [m]
        self.h_engine = h_engine  # height of engine [m]
        self.l_engine = l_engine  # length of engine [m]
        self.D_nose = D_nose  # diameter of nose tire [m]
        self.h_nose = h_nose  # height of nose tire [m]
        self.D_main = D_main  # diameter of main tire [m]
        self.h_main = h_main  # height of main tire [m]
        self.h_nose_strut = h_nose_strut  # height of nose strut [m]
        self.h_main_strut = h_main_strut  # height of main strut [m]
        self.l_main_lateral = l_main_lateral  # main landing gear position when opened
        self.l_long_nose = l_long_nose  # nose landing gear longitudinal position
        self.l_long_main = l_long_main  # main landing gear longitudinal position
        self.h_battery = h_main  # battery height
        self.l_tailcone = l_tailcone  # length of tailcone
        self.h_tail = h_tail  # height of tail
        self.V_battery = V_battery  # volume of battery [m^3]
        
    def n_row(self):
        return math.ceil(self.n_seat / FuselageSizing.n_row_seat)
    
    def l_nosecone(self):
        return self.l_engine + FuselageSizing.l_enbu
    
    def l_cabin(self):
        return FuselageSizing.l_pax * self.n_row()
    
    def length_main_strut(self, s_gear):
        l_lateral_strut = (self.D_main / 2) + s_gear
        h = self.h_main_strut
        l = self.l_main_lateral - l_lateral_strut
        return np.sqrt(l**2 + h**2) 
    
    def l_battery(self, w_battery):
        return self.V_battery / (self.h_battery * w_battery)
    
    def w_battery(self, l_battery):
        return self.V_battery / (self.h_battery * l_battery)
    
    def l_end_nose_land(self):
        return self.l_long_nose + self.h_nose_strut + (self.D_nose / 2) + FuselageSizing.s_battery
    
    def l_end_main_land(self, s_gear):
        return self.l_long_main - self.length_main_strut(s_gear) - (self.D_main / 2)
        
    def battery_dim(self, s_gear):
        l_end_main_land = self.l_end_main_land(s_gear)
        l_end_nose_land = self.l_end_nose_land()
        
        if l_end_main_land - l_end_nose_land < 0:
            print('Landing gear folds backward')
            l_end_main_land = self.l_long_main + self.length_main_strut(s_gear) + (self.D_main / 2)
            w_battery = s_gear * 2 + (self.D_main / 2)
            l_battery = self.l_battery(w_battery)
            
            while self.length_main_strut(s_gear) - (self.D_main / 2) - l_battery < 0:
                s_gear += 0.01 
                w_battery = s_gear * 2 + (self.D_main / 2)
                l_battery = self.l_battery(w_battery)
            
            l_fus = self.length_fus()
            
            if l_fus - self.l_tailcone < l_end_main_land:
                print('Landing gear cannot be folded backward or forward')
                
        else:
            print('Landing gear folds forward')
            w_battery = s_gear * 2
            l_battery = self.l_battery(w_battery)

        return round(l_battery, 3), round(w_battery, 3), round(s_gear, 3)
    
    '''Overall dimensions'''
    def top_width(self):
        return round((FuselageSizing.w_pax + FuselageSizing.s + FuselageSizing.t_fuse) * 2, 3)
    
    @staticmethod
    def bigger_mag(a, b):
        return a if a > b else b
        
    def bottom_width(self, s_gear):
        # s_gear is the distance between the two areas 
        l_battery, w_battery, s_gear = self.battery_dim(s_gear)
        # bottom width from seats:
        bw1 = (FuselageSizing.w_pax + FuselageSizing.s + FuselageSizing.t_fuse) * 2
        
        # bottom width from landing gear: 
        bw2 = 2 * (2 * s_gear + self.D_main + FuselageSizing.s_aft)
    
        return round(FuselageSizing.bigger_mag(bw1, bw2), 3)
    
    def height(self): 
        return round(FuselageSizing.h_pax + (2 * FuselageSizing.t_fuse) + FuselageSizing.h_floor + FuselageSizing.s_bat_wheel + self.h_main, 3)

    def length_fus(self):
        return round(self.l_tailcone + self.l_nosecone() + self.l_cabin() + FuselageSizing.l_cock, 3)
    
    def fuselage_wetted(self, s_gear):
        # assuming trapezoidal shape 
        h_fus = self.height()
        l_fus = self.length_fus()
        wb_fus = self.bottom_width(s_gear)
        wu_fus = self.top_width()
        
        # area of trapezium
        A_trap = 0.5 * (wb_fus + wu_fus) * h_fus
        
        # length of the side of trapezium
        side_c = (((wb_fus - wu_fus)/2)**2 + h_fus**2)**0.5
        l_p = wb_fus + wu_fus + 2*side_c
        a_lat = l_p*l_fus
        
        return 2*A_trap + a_lat
        

# Example usage:
n_seat = 8
w_engine = 1.48
h_engine = 0.15
l_engine = 0.90 
D_nose = 0.33 
h_nose = 0.13
D_main = 0.6
h_main = 0.19
h_nose_strut = 0.89
h_main_strut = 0.89
l_main_lateral = 1.93
l_long_nose = 0.777
l_long_main = 2.09 
l_tailcone = 5.02 
h_tail = 2.52 
V_battery = 0.2

fuselage_size = FuselageSizing(n_seat, w_engine, l_engine, h_engine, D_nose, h_nose, D_main, h_main, h_nose_strut, h_main_strut, l_main_lateral, l_long_nose, l_long_main, l_tailcone, h_tail, V_battery)

print(fuselage_size.top_width())
print(fuselage_size.bottom_width(s_gear=0.2))
print(fuselage_size.height())
print(fuselage_size.length_fus())
print('fuselage wetted area', fuselage_size.fuselage_wetted(s_gear=0.2))
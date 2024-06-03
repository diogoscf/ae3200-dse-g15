import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

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
    
    def __init__(self, n_seat, w_engine, l_engine, h_engine, s_engine, D_nose, h_nose, D_main, h_main, h_nose_strut, h_main_strut, l_main_lateral, l_long_nose, l_long_main, l_tailcone, h_tail, V_battery): 
        self.n_seat = n_seat  # number of total seats
        self.w_engine = w_engine  # width of engine [m]
        self.h_engine = h_engine  # height of engine [m]
        self.l_engine = l_engine  # length of engine [m]
        self.s_engine = s_engine # small space in front of engine [m]
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
    
    def check_back(self,s_gear):
        l_end_main_land = self.l_end_main_land(s_gear)
        l_end_nose_land = self.l_end_nose_land()
        
        return True if  l_end_main_land - l_end_nose_land < 0 else False
        
    def battery_dim(self, s_gear):        
        if self.check_back(s_gear):
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
            #print('Landing gear folds forward')
            w_battery = (self.top_width()-2*FuselageSizing.s)* 2
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
        bw2 = 2 * (s_gear + self.D_main + FuselageSizing.s_aft + FuselageSizing.s)
    
        return round(FuselageSizing.bigger_mag(bw1, bw2), 3)
    
    def height(self): 
        return round(FuselageSizing.h_pax + (2 * FuselageSizing.t_fuse) + FuselageSizing.h_floor + FuselageSizing.s_bat_wheel + self.h_main + FuselageSizing.s, 3)
    
    def y_floorheight(self):
        return FuselageSizing.t_fuse + FuselageSizing.s_bat_wheel + self.h_main

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
    
    def maximum_perimeter(self, s_gear):
        # assuming trapezoidal shape 
        h_fus = self.height()
        wb_fus = self.bottom_width(s_gear)
        wu_fus = self.top_width()
                
        # length of the side of trapezium
        side_c = (((wb_fus - wu_fus)/2)**2 + h_fus**2)**0.5
        l_p = wb_fus + wu_fus + 2*side_c        
        return l_p
    
    def plot_side_drawing(self,s_gear):
        fig, ax = plt.subplots()
        
        y_floor = self.y_floorheight()
        y_engine = y_floor + self.h_floor
        ax.axhline(y=0, color='gray', linewidth=0.3)
        
        # Define the position and size of the rectangles
        rectangles = [
            patches.Rectangle((FuselageSizing.t_fuse, y_floor), self.length_fus() - 2*FuselageSizing.t_fuse, self.h_floor, linewidth=0.5, edgecolor='b', facecolor='none', label='floor'),  # Floor
            patches.Rectangle((FuselageSizing.t_fuse + self.s_engine, y_engine), self.l_engine, self.h_engine, linewidth=0.5, edgecolor='r', facecolor='none', label='engine'),  # Engine
            patches.Rectangle((FuselageSizing.t_fuse + self.s_engine + self.l_engine, y_engine), self.l_cock, self.h_cock, linewidth=0.5, edgecolor='g', facecolor='none', label='cockpit'),  # Cockpit
            patches.Rectangle((0,0), self.length_fus(), self.height(), linewidth=0.5, edgecolor='c', facecolor='none'), # outer fuselage
            patches.Rectangle((FuselageSizing.t_fuse, FuselageSizing.t_fuse), self.length_fus()-2*FuselageSizing.t_fuse, self.height()-2*FuselageSizing.t_fuse, linewidth=0.5, edgecolor='c', facecolor='none'),
            patches.Rectangle((self.l_long_nose, FuselageSizing.t_fuse + self.h_nose/2), self.h_nose_strut, 0.01, linewidth=0.5, edgecolor='k', facecolor='none'), # nose landing strut
            patches.Rectangle((self.l_long_nose+self.h_nose_strut-self.D_nose/2, FuselageSizing.t_fuse), self.D_nose, self.h_nose, linewidth=0.5, edgecolor='m', facecolor='none') # nose landing gear
        ]
        
        # Add passengers' rectangles
        num_passengers = self.n_row()
        for i in range(num_passengers):
            passenger_x = self.t_fuse+self.s_engine + self.l_engine + self.l_cock + (self.l_pax) * i
            rect = patches.Rectangle((passenger_x, y_engine), self.l_pax, self.h_pax, linewidth=0.5, edgecolor='k', facecolor='none')
            rectangles.append(rect)
            
        if self.check_back(s_gear):
            l_battery, w_battery, s_gears = self.battery_dim(s_gear)
            main_land_strut = patches.Rectangle((self.l_long_main, FuselageSizing.t_fuse + self.h_main/2), self.length_main_strut(s_gear), 0.01, linewidth=0.5, edgecolor='k', facecolor='none')
            main_land_gear = patches.Rectangle((self.l_long_main+self.length_main_strut(s_gear)-self.D_main/2, FuselageSizing.t_fuse), self.D_main, self.h_main, linewidth=0.5, edgecolor='m', facecolor='none')
            battery = patches.Rectangle((self.l_long_main, FuselageSizing.t_fuse), l_battery, self.h_battery, linewidth=0.5, edgecolor='g', facecolor='none')
            rectangles.append(main_land_strut)
            rectangles.append(main_land_gear)
            rectangles.append(battery)
            #print('s_gears', s_gears)
            #print('w_battery', w_battery)
        
        else:
            l_battery, w_battery, s_gears = self.battery_dim(s_gear)
            main_land_strut = patches.Rectangle((self.l_long_main, FuselageSizing.t_fuse + self.h_main/2), -self.length_main_strut(s_gear), 0.01, linewidth=0.5, edgecolor='k', facecolor='none')
            main_land_gear = patches.Rectangle((self.l_long_main-self.length_main_strut(s_gear)-self.D_main/2, FuselageSizing.t_fuse), self.D_main, self.h_main, linewidth=0.5, edgecolor='m', facecolor='none')
            battery = patches.Rectangle((self.l_end_nose_land()+FuselageSizing.t_fuse, FuselageSizing.t_fuse), l_battery, self.h_battery, linewidth=0.5, edgecolor='g', facecolor='none')
            rectangles.append(main_land_strut)
            rectangles.append(main_land_gear)
            rectangles.append(battery)
            #print('s_gears', s_gears)
            #print('w_battery', w_battery)
            
        # Add the rectangles to the axes
        for rect in rectangles:
            ax.add_patch(rect)
        
        #ax.set_xlim(-0.05, self.length_fus() + 0.2)
        #ax.set_ylim(-0.05, self.height() + 0.2)
        
        plt.legend()
        plt.title('Fuselage side view')
        plt.axis('equal')
        plt.show()
    
    def plot_front_view(self, s_gear):
        fig, ax = plt.subplots()
        ax.axhline(y=0, color='gray', linewidth=0.3)

        # vertices of the trapezoid
        
        vertice1 = [(-self.bottom_width(s_gear)/2, self.h_nose_strut),
                    (self.bottom_width(s_gear)/2, self.h_nose_strut),
                    (self.top_width()/2, self.h_nose_strut + self.height()),
                    (-self.top_width()/2, self.h_nose_strut + self.height())]
        vertice2 = [(-self.bottom_width(s_gear)/2+FuselageSizing.s, self.h_nose_strut + FuselageSizing.s),
                    (self.bottom_width(s_gear)/2 -FuselageSizing.s, self.h_nose_strut + FuselageSizing.s),
                    (self.top_width()/2 - FuselageSizing.s, self.h_nose_strut + self.height()-FuselageSizing.s),
                    (-self.top_width()/2 + FuselageSizing.s, self.h_nose_strut + self.height()-FuselageSizing.s)]
        
         
        l_battery, w_battery, s_gear = self.battery_dim(s_gear)
                             
        trapezoid1 = patches.Polygon(vertice1, closed=True, edgecolor='c', facecolor='none', linewidth=0.5)
        trapezoid2 = patches.Polygon(vertice2, closed=True, edgecolor='c', facecolor='none', linewidth=0.5)
        #circle1 = patches.Circle((self.l_main_lateral, self.D_main/2), self.D_main/2, edgecolor='r', facecolor='none')
        #circle2 = patches.Circle((-self.l_main_lateral, self.D_main/2), self.D_main/2, edgecolor='r', facecolor='none')
        circle3 = patches.Circle((s_gear + self.D_main/2, self.h_nose_strut + FuselageSizing.s + self.h_main/2), 0.05, edgecolor = 'r', facecolor='none')
        circle4 = patches.Circle((-s_gear - self.D_main/2,  self.h_nose_strut + FuselageSizing.s + self.h_main/2), 0.05, edgecolor = 'r', facecolor='none')
        
        rectangles = [
            patches.Rectangle((self.l_main_lateral - self.h_main/2, 0), self.h_main, self.D_main, edgecolor = 'm', facecolor='none'), 
            patches.Rectangle((-self.l_main_lateral - self.h_main/2,0), self.h_main, self.D_main, edgecolor = 'm', facecolor='none'),
            patches.Rectangle((-self.top_width()/2, self.y_floorheight()+self.h_nose_strut), self.top_width(), self.h_floor, linewidth=0.5, edgecolor='b', facecolor='none')
            ]
        
        
        # Line between two points:
        x_val1 = [self.l_main_lateral, s_gear + self.D_main/2]
        x_val2 = [-self.l_main_lateral, -s_gear - self.D_main/2]
        y_val1 = [self.D_main/2, self.h_nose_strut + FuselageSizing.s + self.h_main/2]
        
        for rect in rectangles:
            ax.add_patch(rect)
        
        ax.plot(x_val1, y_val1, 'k-', linewidth=1)
        ax.plot(x_val2, y_val1, 'k-', linewidth=1)
        
        ax.add_patch(trapezoid1)
        ax.add_patch(trapezoid2)
        #ax.add_patch(circle1)
        #ax.add_patch(circle2)
        ax.add_patch(circle3)
        ax.add_patch(circle4)
        ax.axvline(x=0, color='gray', linestyle='--', linewidth=1)
        plt.axis('equal')
        plt.title('Fuselage front view')
        plt.show()       
             


# Example usage:
n_seat = 8
w_engine = 0.85
h_engine = 0.75
l_engine = 1.1
s_engine = 0.1
D_nose = 0.329565
h_nose = 0.127
D_main = 0.55626
h_main =  0.2286
h_nose_strut = 0.65
h_main_strut = 0.65
l_main_lateral = 1
l_long_nose = 0.40874
l_long_main = 4.21989077
l_tailcone = 5.02 
h_tail = 2.52 
V_battery = 0.3

fuselage_size = FuselageSizing(n_seat, w_engine, l_engine, h_engine, s_engine, D_nose, h_nose, D_main, h_main, h_nose_strut, h_main_strut, l_main_lateral, l_long_nose, l_long_main, l_tailcone, h_tail, V_battery)

print('top_width', fuselage_size.top_width())
print('bottom_width',fuselage_size.bottom_width(s_gear=0.2))
print('fuselage height', fuselage_size.height())
print('fuselage_length without engine', fuselage_size.length_fus()-l_engine-fuselage_size.l_enbu)
print('fuselage length', fuselage_size.length_fus())
print('maximum perimeter', fuselage_size.maximum_perimeter(s_gear=0.2))
#print('fuselage wetted area', fuselage_size.fuselage_wetted(s_gear=0.2))
print('main_strut_length',fuselage_size.length_main_strut(s_gear=0.2))
print('nose_strut_length',h_nose_strut)

fuselage_size.plot_side_drawing(s_gear=0.2)
fuselage_size.plot_front_view(s_gear=0.2)

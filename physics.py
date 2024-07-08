import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
import time as tm
import os as os # I like adding bullshit like this to my code sometimes

T = tm.localtime() # Pull time at which program is run
filename = f"Logs/Logs_{T[7]}_{T[3]}-{T[4]}-{T[5]}.txt" # Set log file name

# Define logging functions
def log(text):
    try:
        T = tm.localtime()
    except:
        print("[Error]: Library `time` imported improperly!")
    with open(filename , 'a') as log_file:
        log_file.write(f"\n{T[7]} - {T[3]}:{T[4]}:{T[5]} > {text} \n")
    print(f"{T[3]}:{T[4]}:{T[5]} > {text}")        

def backend_log(text):
    try:
         tm.localtime()
    except:
        print("[Error]: Library `time` imported improperly!")
    with open(filename , 'a') as log_file:
        log_file.write(f"\n{T[7]} - {T[3]}:{T[4]}:{T[5]} > {text} \n")

def debug(text):
    with open('debug_log.txt' , 'a') as log_file:
        log_file.write(f">>> {text}\n")

log("Imports completed")
log(f"Logs saved in file {filename}")

# Define mathematical tools

def vector_add(v1 , v2):
    v_res = np.array([v1[0] + v2[0] , v1[1] + v2[1]])
    return v_res

def vector_sub(v1 , v2):
    v_res = np.array([v2[0] - v1[0] , v2[1] - v1[1]]) # subtracts v1 from v2 (v2 -v1)
    return v_res

def vector_mag(v):
    mag = np.sqrt(v[0]**2 + v[1]**2)
    return mag

# Define fundamental constants of nature

G = 6.6743e-11 
k = 8.99e+10

# Define fundamental forces of nature

def N_Gravity(G , m1 , m2 , r1, r2): # Gravitational force extered on m1 by m2     
    displacement = r2 - r1
    disp_mag = vector_mag(displacement)
    if disp_mag < 0.01:
        log("[Incident]: Collsion detected! - by gravity function")
    
    F_mag = (-1) * G * m1 * m2 / disp_mag**2
    F_head = displacement/ disp_mag
    
    F = F_mag * F_head
    return F
    
def Electro(k , q1 , q2 , r1, r2): # Coloumb's Law of Electrostatics 
    displacement = r2 - r1
    disp_mag = vector_mag(displacement)
    if disp_mag < 0.01:
        log("[Incident]: Collsion detected! - by electrostatic function")
    
    F_mag = (1) * k * q1 * q2 / disp_mag**2
    F_head = displacement/ disp_mag
    
    F = F_mag * F_head
    return F

# skip nuclear forces

# Define force interactions and gallilean transforms
def delta_v(a , u , t):
    v = u + a*t
    return v

def delta_x(x , v , t):
    X = x + v*t
    return X

def acceleration(n):
    backend_log(f"[Debug]: Acceleration compuation for element {n} begins")
    properties = System[n]
    F_g = np.array([0.0,0.0])
    F_e = np.array([0.0,0.0])
    
    for i in range(len(System)):
        if i == n:
            backend_log(f"[Debug]: Skipping acceleration compuation of element {n} on itself ({i})")
            continue # Skip interactions with itself 
        F_g += N_Gravity(G, properties['mass'], System[i]['mass'] , properties['position'], System[i]['position'])
        F_e += Electro(G, properties['charge'], System[i]['charge'] , properties['position'], System[i]['position'])
        
    F_net = F_g + F_e
    backend_log(f"[Debug]: Acceleration Compuation for {n}: F_g = {F_g} , F_e = {F_e} , F_net = {F_net}")
    debug(f"Total force is {F_net} \n ======")
    
    a = F_net / properties['mass']
    return a

def Update(n, delta_t):
    backend_log(f"[Element_Cursor_Debug]: Logging for element number {n}")
    properties = System[n]
    X = delta_x(properties['position'] , properties['velocity'] , delta_t) # Compute updated position
    
    a = acceleration(n) # Compute acceleration
    backend_log(f"[Debug]: Acceleration computation complete for element number {n}")
    V = delta_v(a , properties['velocity'] , delta_t) # Compute updated velocity 
    properties['position'] = X # Update properties
    properties['velocity'] = V
    backend_log(f"[Debug]: Updated properties for element number {n} : X = {X} ; V = {V}")
    return X , V, a
      
    
log("Functions defined")    


# Define elements in system:
    
System_Archive = [ # Holds all the obkects
   { 
    "mass" : 1, 
    "charge" : 50000,
    "velocity" : np.array([0.0 , 0.0]),
    "position" : np.array([-1.0,0.0]),
    "color" : 'green'
    },
   {
    "mass" : 1,
    "charge" : -50000,
    "velocity" : np.array([0.0 , 0.0]),
    "position" : np.array([1.0,0.0]),
    "color" : 'aqua'    
    },
   {
    "mass" : 10,
    "charge" : 30000,
    "velocity" : np.array([0.5 , 0.0]),
    "position" : np.array([-5.0,0.0]),
    "color" : 'aqua'    
    },
]

System = [] # System list pulls seleted objects from System_Archive using loops defined below

s = 1 # Select trunkator

for i in range(len(System_Archive) - s): # Pull data from System_Archive and generate a list System 
    System.append(System_Archive[i])
t = np.linspace(0 , 5 , 5 * 100)

log("System defined")

System_data = [] # Holds positions for each object in the system for each instance of time 
System_data_helper = [] # Holds accelerations of each object in the system for each instance of time

for _ in t: # Add rows ro System_data and System_data_helper required for each time instance
    System_data.append([])
    System_data_helper.append([])
    
for i in range(len(System) - s):
    System_data[0].append(System[i]['position']) # Append initial positions into System_data
    System_data_helper[0].append(np.array([0.0 , 0.0])) # Append initial zeros into System_data_helper
log("System data list defined")

for j in range(len(t)): # Computes System for each instance of time
    backend_log(f"[Time_Debug]: Running on time instance {t[j]} with delta_t = {t[j] - t[j-1]}")
    delta_t = t[j] - t[j - 1]
    for n in range(len(System)): # Computes interactions on each object by the rest of the objects in the system
        debug(f"At time instance {t[j]}:")
        X , V , a = Update(n , delta_t ) # Pull updated values of position X and velocity V and acceleration a from Update() function
        System_data[j].append(X) # Append the updated position to the System_data for current time instance
        System_data_helper[j].append(a) # Append the updated acceleration to the System_data_helper for current time instance
    backend_log(f"[Time_Debug]: End time instance {t[j]}")
    
d = 3 # Dimensions of plot

log("Plotting diagrams...")

alpha = 2 # Scaling factor for plots should it be useful

for i in range(len(t)):
    plt.figure() # Set up plot
    plt.xlim(-d , d)
    plt.ylim(-d , d)
    plt.grid(True)
    plt.title(f"Time instance {t[i]}")
    for j in range(len(System_data[i])): # Pull data from System_data and System_data_helper to plot
        #print(j)
        x = System_data[i][j][0]
        y = System_data[i][j][1]
        
        a_x = System_data_helper[i][j][0]
        a_y = System_data_helper[i][j][1]
        
        Ex = [x , (a_x + x) * alpha]
        Why = [y , (a_y + y) * alpha]
        
        try:
            if System[j]['color'] == '$code_charge': # Assign colors based on charge         
                try:
                    if System[j]['charge'] == 0:
                        cl = 'grey'  # neutral charge
                        #print(cl)
                    elif System[j]['charge'] <= 0:
                        cl = 'orange'  # negative charge
                        #print(cl)
                    else:
                        cl = 'blue'  # positive charge"""
                        #print(cl)
                except IndexError as e:
                    log(f"[Error]: {e}, ignored and proceed with special cases") # To handle errors that happen here for some reason
                    cl = 'red' # colour code erronious points
            else:
                cl = System[j]['color'] # Pull colour setting from system object properties
        except:
            log(f"[Error]: Unknown indexing error (?)") # Handle an indexing error that happens for some god forsaken reason
            cl = 'red' # colour code erronious points
        try:
            plt.scatter(x, y, color = cl) # Plot the position
        except:
            log(f"[Error]: The color {cl} defined for object {System[j]} is not recognised")
        plt.plot(Ex, Why , color = 'red') # Plot the scaled acceleration vector
    try:
        plt.savefig(f"Results/Frame_{i}.png") # Save plots into a folder Results
    except:
        log("[Error]: Results folder doesn't exist")
        os.makedirs('Results') # Make folder if it doesn't exist
        log("Results folder made using os.makedirs")
        plt.savefig(f"Results/Frame_{i}.png") # Save plots into a folder Results
    plt.show() # Make plot visible
log("Finishing plotting figures. Figures saved to file Results")
     
    
log("End program")

import math
import numpy as np
import matplotlib.pyplot as plt


def verlet(t, theta, r_old, v_old, mass, time, max_alt, proj_on, check):                                                                                  # Verlet Algorithm for 2D projectile motion with modified Force function to simulate lift
    
    
    x_new = r_old[0] + v_old[0]*t + 0.5*(t**2)*Force(0, r_old, v_old, mass, time, max_alt, proj_on)/mass                                                  # Eq. 1 of Verlet Algorithm for new position
    y_new = r_old[1] + v_old[1]*t + 0.5*(t**2)*Force(1, r_old, v_old, mass, time, max_alt, proj_on)/mass
    
    r_new = [x_new, y_new]
    
    
    if max_alt[4] == True and check == False:
        check = check_intersection(max_alt[3], r_new)
        if check == True:
            x_mid = (max_alt[2] + r_new[0])/2
            y_mid = (max_alt[0] + r_new[1])/2
            max_alt.insert(5,x_mid)

    if max_alt[4] == False:
        max_alt = max_glide(r_new, r_old, max_alt, v_old)                                                                                                 # Finds when maximum altitude is reached and then solves for max glide range for case with constant velocity and acceleration (Used for reference)
    
   
    
    vx_new_update = v_old[0] + 0.5*t/mass*(Force(0, r_old, v_old, mass, time, max_alt, proj_on) + Force(0, r_new, v_old, mass, time, max_alt, proj_on))   # Eq. 2 of Verlet Algorithm for updated velocity
    vy_new_update = v_old[1] + 0.5*t/mass*(Force(1, r_old, v_old, mass, time, max_alt, proj_on) + Force(1, r_new, v_old, mass, time, max_alt, proj_on))
    
    v_new_update = math.sqrt(vx_new_update**2 + vy_new_update**2)                                                                                         # Solving for updated velocity magnitude and the new flight path angle
    v_new_update = [vx_new_update, vy_new_update, v_new_update] 
    
    vx_new = v_old[0] + 0.5*t/mass*(Force(0, r_old, v_old, mass, time, max_alt, proj_on) + Force(0, r_new, v_new_update, mass, time, max_alt, proj_on))   # Eq. 3 of Verlet Algorithm for final new velocity 
    vy_new = v_old[1] + 0.5*t/mass*(Force(1, r_old, v_old, mass, time, max_alt, proj_on) + Force(1, r_new, v_new_update, mass, time, max_alt, proj_on))
    
    v_new = math.sqrt(vx_new**2 + vy_new**2)                                                                                                              # Solving for new velocity magnitude and the new flight path angle
    v_new = [vx_new, vy_new, v_new]


    ## Verlet formulas for reference ##
	# r_new = r_old + v_old*t + 0.5*(t**2)*Force(r_old, v_old)/mass
	# v_new_update = v_old + 0.5*t/mass*(Force(r_old, v_old) + Force(r_new, v_old))
	# v_new = v_old + 0.5*t/mass*(Force(r_old, v_old) + Force(r_new, v_new_update))
    
    
    return r_new, v_new, max_alt, check

def Force(x, r, v, mass, time, max_alt, proj_on):                                           # Since equations of motion do not contain an x term, x will flag which equation of motion is being used (ie y direction has -mg term)

        y = r[1]
        cd = c(y)                                                                           # Approximates form drag on an object relative to altitude
        L = 3*cd                                                                            # Calculates lift from quadratic drag & from L/D = 3 => L = 3*D
        theta = math.atan2(v[1], v[0])
        
        if time <= 2:
            T = mass*9*9.81                                                                 # Thrust as a calculation of 9-g's acceleration
            if proj_on == 0:
                if x == 0:                                                                  # Checking flag for x/y direction, then calculates drag force from the Quadratic Eqs of Motion
                    f = T*math.cos(theta) - L*math.sin(theta) - cd*v[2]*v[0]
                    return f
                elif x == 1:
                    f = T*math.sin(theta) + L*math.cos(theta) - cd*v[2]*v[1] - (mass*9.81)
                    return f
                    
            elif proj_on == 1:
                if x == 0:                                                                  # Checking flag for x/y direction, then calculates drag force from the Quadratic Eqs of Motion
                    f = T*math.cos(theta) - cd*v[2]*v[0]
                    return f
                elif x == 1:
                    f = T*math.sin(theta) - cd*v[2]*v[1] - (mass*9.81)
                    return f
        else:
            if proj_on == 0:
                if x == 0:                                                                  # Checking flag for x/y direction, then calculates drag force from the Quadratic Eqs of Motion
                    f = -L*math.sin(theta) -cd*v[2]*v[0]
                    return f
                elif x == 1:
                    f =  L*math.cos(theta) -cd*v[2]*v[1] - (mass*9.81)
                    return f
            elif proj_on == 1:
                if x == 0:                                                                  
                    f = -cd*v[2]*v[0]
                    return f
                elif x == 1:
                    f = -cd*v[2]*v[1] - (mass*9.81)
                    return f

def c(y):
    area = 1                # Normalizing area of glider to value of 1
    lambda_sym = 10000      # Constant approximated from altitude tables which is used to find air density 
    gamma = 0.25            # Constant used for quadratic drag coefficient
    
    cd = gamma*(area**2)*math.exp((-y/lambda_sym))
    return cd

def max_glide(new, old, flag, v):
    
    max_alt[4] = True if new[1] < old[1] and max_alt[4] == False else False                                                 # Determines when altitude is at peak height
    
    
    theta_check = math.atan2(v[1], v[0])
    FPA_min = math.atan(1/3)

     
   
    
    if max_alt[4] == True and theta_check <= FPA_min:                                                                       # Solves for maximum range from glider equations using L/D = 3 
        R_max = old[1] * 3
        Total = old[0] + R_max
        print("Peak Altitude: ", old[1], "m")
        
        max_alt.insert(0,old[1])                                                                                            # Inserting values for plotting theoretical maximum glide range
        max_alt.insert(1,0)
        max_alt.insert(2,old[0])
        max_alt.insert(3,Total)
        max_alt.insert(4,True)
        
    return max_alt                                                                                                          # Returning flag for detecting when peak altitude is met                

def check_intersection(range, r_new):
    y = -(1/3)*r_new[0] + (1/3)*range                       # linear glide range function, used for determining intersection of projectile motion of glider

    if abs(r_new[1] - y) < 0.001:
        check = True
    else:
        check = False
    return check






## Setting flight path angle and mass ##

step = 15;   # Step size for incident angle          
theta_deg = np.arange(15, 76, step)
theta = 2*math.pi/360*theta_deg;
colors = ["b", "g", "red", "m", "c" , "y", "m", "k"] 

mass = input("Please input desired mass for vehicle [kg]: ")

try:
    mass = int(round(float(mass)))
except:
    mass = 10000 # default mass



glide_on = input("Do you want to display glider range from peak Altitude? [Y/N]: ")

if glide_on == "Y":
    glide_on = True
else:
    glide_on = False
    

proj_on = input("Do you want to display standard projectile motion? [Y/N]: ")

if proj_on == "Y":
    proj_on = True
else:
    proj_on = False

count = 0


for angle in theta:
    print("\nPlotting Results for Flight Path Angle: ", int(round(float(angle*360/(2*math.pi)))), "°")
    
    # Setting initial position and velocity
    r_new = 0; x_new = 0; y_new = 0;
    v_new = 88.29*2                                                                                     # Kinematic Equation vf = vi + a*t = 88.29*2 for setting verlet algorithm initial velocity 

    # Breaking velocity into its vector components
    vx_new = v_new*math.cos(angle)
    vy_new = v_new*math.sin(angle)

    v_new = [vx_new, vy_new, v_new]                                                                     # Creating velocity array for storing x/y components and magnitude
    r_new = [x_new, y_new]                                                                              # Creating position array for storing x/y coordinate
    r_new2 = r_new; v_new2 = v_new;         
    
    x_glide = []; y_glide = [];                                                                         # Creating empty arrays for plotting maximum glide range

    t_step = 0.0001                                                                                     # Setting time step to 0.001 seconds for accurate range, and initial time to 2 seconds for climb stage after thrust
    t = 2
    x_values = []; y_values = []; x_values.append(x_new); y_values.append(y_new);                       # Creating arrays for storing simulated data to eventually plot
    x_values2 = []; y_values2 = []; x_values2.append(x_new); y_values2.append(y_new);


    max_alt = [0,0,0,0,False,0,0]
    max_alt2 = [0,0,0,0,False,0,0]
    check = False
    range = 1                                                                                                         # Setting range as a placeholder of 1 for if statement in while loop

    while t <= 60 and r_new[1] >= 0:                                                                                  # while loop to ensure that the simulation runs for 60 seconds and ends once it hits the ground
        r_new, v_new, max_alt, check = verlet(t_step, angle, r_new, v_new, mass, t, max_alt, 0, check)                # calling custom verlet function which updates inputs as the output calculation each loop iteration
        
        if proj_on == True:
            r_new2, v_new2, max_alt2, check = verlet(t_step, angle, r_new2, v_new2, mass, t, max_alt2, 1, check)      # Verlet simulation for standard projectile motion
            x_values2.append(r_new2[0])                                                                    
            y_values2.append(r_new2[1])
        
        x_values.append(r_new[0])                                                                                     # appending data to arrays for plotting
        y_values.append(r_new[1])
        
        if r_new[1] <= 0 and range == 1:                                                                              # If statement for determining horizontal range of x(t), when y(t) <= 0 the simulated range integer is found
            range = r_new[0]
            print("Simulated Horizontal Range:", range, "m", "At height value of:", r_new[1], "m")

            
        t = t + t_step
    print("Duration: ", t, "s")    
    
    
    
    legend_txt = "Flight Path Angle of " + str(round(float(angle*360/(math.pi*2)))) + "°"
    

    x_values = np.array(x_values)
    y_values = np.array(y_values)
    idx = np.argmin(np.abs(x_values - max_alt[5]))
        
    tan_height = y_values[idx]
    tan_position = x_values[idx]

    glide_range = tan_height*3 + tan_position                                        # solves for maximum range with accounting for glider momentum (tangent of parabolic trajectory)
    
    if max_alt[5] == 0:
        print("Total Distance: ", r_new[0], "m \n")
    else:
        print("Total Distance: ", glide_range, "m \n")
        
    plt.plot([max_alt[5], glide_range], [tan_height, 0], color = colors[count])      # Plots glider trajectory from tangent of parabolic trajectory
     
    if angle >= 30*math.pi*2/360:
        x_values_new = x_values[0:idx]
        y_values_new = y_values[0:idx]
    else:
        x_values_new = x_values
        y_values_new = y_values
        
    if max_alt[5] == 0 and glide_on == True:
        x_values_new = x_values
        y_values_new = y_values
        
    elif max_alt[5] != 0 and glide_on == True:                                       # Flag for plotting glider range from peak altitude 
        
        x_glide = [max_alt[2], max_alt[3]]
        y_glide = [max_alt[0], max_alt[1]]
        
        plt.plot(x_glide, y_glide, linestyle = '--', color = 'k')                    # Plots glider trajectory from peak altitude
        
        
    if proj_on == True:
        plt.plot(x_values2, y_values2, linestyle = '--', color = 'gray')             # Plots standard projectile motion
    
    plt.plot(x_values_new, y_values_new, label = legend_txt, color = colors[count])  # Plots simulated data
    count = count + 1
    
plt.title('Trajectory of Aircraft')
plt.xlabel('X Distance Traveled (m)')
plt.ylabel('Y Height Traveled (m)')
plt.legend()
plt.show()    
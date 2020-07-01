# This file will be a combination of all the stuff we have modelled so far.
# Some variables defined in individual code have been deleted to help integrate
# each file.

# Import all necessary modules or whatever

import math
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from itertools import product
from scipy.integrate import simps
from scipy.integrate import odeint
from matplotlib import cm

# REALLY IMPORTANT !!! READ BEFORE RUNNING !!!
# Since this program contains multiple parts, I
# organized everything into a function that you
# can call. You'll need to call the respective 
# function in the console or uncomment for it 
# to work. Thank you for reading, especially
# if you commented your code well so I didn't 
# have to stress and make sense of it myself :)

# First, some universal constants/arrays so they don't
# need to be changed or assigned in every individual function

# Constants assumming water
poi = 1 # Poise mPa * s
p = 1 # density g/cm^3

# Dimensions of the channel
W = 700 # width in um (.5mm)
H = 500 # height in um (.7mm)

n1 = 100 #Grid size

delta_P = -20  #change in pressure in (please put units here if you know)
delta_L = 22 # length of path on microfluidic device in mm

# Go to end of code to call functions

# Pouiselle Flow
def PFlow():
    time = 10000 # time of the liquid in us (micro seconds)
    
    # Arrays that need to be initilized
    z = np.linspace(0, W, n1) # Width axis array
    y = np.linspace(0, H, n1) # Height axis array
    
    y, z = np.meshgrid(y, z) # mesh y and z values

    vx = np.zeros((len(y), len(z))) # Initializing another array for vx as a 2D array

    # Here is where things get tricky. So for some reason I can't get the for loop
    # to match with the array list thingy we need to manually do the start and the
    # stop values which should be the same as the z and y start and stops as
    # previously defined. Also the increment should be 1 to match the vx array
    # dimensions. Basically these values correspond to the y and z arrays but I
    # couldn't do it with the actual arrays for some reason.
    # the two summation series
    for y1 in range(0, 700, 7):
        for z1 in range(0, 700, 7):
            temp = 0

            # This calculates the value given the equation from paper.
            # n and m variables are summation
            # I only went from 0 - 10 assumming they get infitiely small when values
            # get larger.
            for n in range(0, 10, 1):
                for m in range(0, 10, 1):
                    temp += -(16 / (poi * math.pi**2)) * (delta_P / delta_L) * \
                    (1 - math.exp(-(((((2 * n + 1) * math.pi) / W)**2) + ((((2 * m + 1) * math.pi) / H)**2)) * ((poi / p) * time))) * \
                    math.sin((2 * n + 1) * math.pi * (y1 / W)) * math.sin((2 * m + 1) * math.pi * (z1 / H)) * \
                    (((2 * n + 1) * (2 * m + 1) * (((((2 * n + 1) * math.pi) / W)**2) + ((((2 * m + 1) * math.pi) / H)**2))) ** -1)

            vx[int(y1/7), int(z1/7)] = temp # Assigns value in vx array

    # I recommend making it the interactive chart which you can
    # do by going to Tools > Preferances > iPython cousol > Graphics > and changing
    # Graphics Backend from Inline to Automatic. You then have to restart spyder
    fig = plt.figure()
    ax = plt.axes(projection = "3d")

    ax.plot_surface(y, z, vx, cmap="plasma")
    ax.set_xlabel("Height (um)", labelpad=15)
    ax.set_ylabel("Width (um)", labelpad=15)
    ax.set_zlabel("Velocity of Fluid in the x direction (mm/s)", labelpad=15)
    ax.set_title("Velocity of Fluid in X Direction Within Rectangular Channel")

    plt.show()

# Volumetric Flow Rate
    
def volFlow():
    
    time = 1000 # time of the liquid in us (micro seconds)

    # This function just finds the value associated with an y and z coordinate. I 
    # return the answer multiplied by 100 because from the paper it seems like we
    # are a factor of 10 off for mm/s but the later process needs the output to be 
    # in nm/s. So if you want the mm/s x velocity output just divide the answer by
    # 10 instead. To be honest I think there's something wrong with the constants
    # we are using so if someone could figure that out too that would be great!
    def equation(y, z):
        num = 0
        for n in range(0, 10, 1):
            for m in range(0, 10, 1):
    
                num += -(16 / (poi * math.pi**2)) * (delta_P / delta_L) * \
                (1 - math.exp(-(((((2 * n + 1) * math.pi) / W)**2) + ((((2 * m + 1) * math.pi) / H)**2)) * ((poi / p) * time))) * \
                math.sin((2 * n + 1) * math.pi * (y / W)) * math.sin((2 * m + 1) * math.pi * (z / H)) * \
                (((2 * n + 1) * (2 * m + 1) * (((((2 * n + 1) * math.pi) / W)**2) + ((((2 * m + 1) * math.pi) / H)**2))) ** -1)
    
        return (num * 100)

    # Lists that need to be initialized to hold the each y and z coordinate being
    # tested and the sum_array holds the values for each y and z equation output
    y_array = []
    z_array = []
    sum_array = []

    # Now here is where I start explaining myself ugh. So to find the volume under 
    # the curve there's a way where we can split the y and z plane into rectangles
    # then we can find a point in each rectangle that corresponds to the equation.
    # We then find the volume of that rectangle by multiplying the length, width, 
    # and height, think about this with respect to the graph of the equation. It's 
    # a lot like the Riemenns sum thing we did in calc. For example if I divided
    # the y and z plane into 1 um by 1 um squares I then can find the equation 
    # output which is the height, remember the 3d graph? I then can find the small
    # rectangular prisim's volume and if you do this for the whole grid you can get
    # a pretty good estimate for the volume under the 3d curve. For this code the
    # equation value is estimated at the middle of the square. Remember there are 
    # different types of Reimmen's approximation (right, left, midpoint, trapizoid).
    # This is really similar to the midpoint approximation way. 

    y_len = 1 # Rectangle length for y 
    z_len = 1 # Rectangle length for z

    tempy = y_len / 2 # Holds the step value for y
    tempz = z_len / 2 # Holds the step value for z

    # While loops to initilize the y and z array that holds all the possible points
    # that will be tested
    while tempy < W:
        y_array.append(tempy)
        tempy += y_len
    while tempz < H:
        z_array.append(tempz)
        tempz += z_len

    # For loop that gets the the equation value for each point and since we divided
    # the grid into even rectangles, it doesn't have to correspond to the point, so 
    # I just added it to the sum_array.
    for y in y_array:
        for z in z_array:
            sum_array.append(equation(y, z))

    # This actually finds the volume and multiplies the each point of the equation
    # in sum_array with the y length and z width that was set before. Remember the 
    # the equation point acts like the height in this case.
    volume = sum([x * y_len * z_len for x in sum_array])
    mmVolume = volume / (1000**3) # Converts the nm^3/s to mm^3/s

    # Prints the volume in mm^3/s
    print("The volumetric flow rate is " + str(mmVolume) + " mm^3/s")

# Real Time Position
    
def PosRealTime():

    delta_T = 1e-5 #Time step [s]

    # Arrays that need to be initilized
    z = np.linspace(0, W, n1) # Width axis array
    y = np.linspace(0, H, n1) # Height axis array
    
    y, z = np.meshgrid(y, z) # Idk some fancy to to mesh them together??? who knows

    vx = np.zeros((n1, n1, 11, 11)) # Reinitializing array for vx as a 4D array

    perm = list(product(range(11),repeat=2)) #11x11 matrix of terms to add up

    const = -(16 / (poi * math.pi**2)) * (delta_P / delta_L) #Constant term (Depends only on constants input by the user)

    #Time Independent part of the equation, only needs to calculate it once
    def TimeInd():
        for y1 in range(n1):
            for z1 in range(n1):
                temp = np.zeros((11,11))
            
                for n in perm:
                    temp[n] = 0.1 * math.sin((2 * n[0] + 1) * math.pi * (y[y1][y1] / W)) * math.sin((2 * n[1] + 1) * math.pi * (z[z1][z1] / H)) * \
                            (((2 * n[0] + 1) * (2 * n[1] + 1) * (((((2 * n[0] + 1) * math.pi) / W)**2) + ((((2 * n[1] + 1) * math.pi) / H)**2))) ** -1)
                
                vx[y1][z1] = const * temp # Assigns value in vx array
        return vx;

    #Time dependent part of the equation, needs to be updated every time point
    def TimeDepend(t):
        upd = 1 - np.exp(part * t)
        return upd;

    part2 = np.zeros((11,11))

    #Summation of parts in the time independent part so it doesn't need to calculate it more than once (Might be overcomplicating it idk)
    def Expn():
        for k in perm:
            part2[k] = -(((((2 * k[0] + 1) * math.pi) / W)**2) + ((((2 * k[1] + 1) * math.pi) / H)**2)) * (poi / p)
        return part2;

    part = Expn()

    independ = TimeInd()

    newArr = np.zeros((n1,n1))
    distTravel = np.zeros((n1,n1))

    timePoints = 1001 #Number of time points to model


    distLapse = np.zeros((timePoints,n1,n1)) #Array holding distance travelled at every time point
    velLapse = np.zeros((timePoints,n1,n1)) #Array holding velocity at every time point

    #Just creates the figure
    fig1 = plt.figure(figsize=(14,10))
    ax = plt.axes(projection = "3d")
    
    for x in range(timePoints):
        for k in range(n1): #y loop
            for l in range(n1): #z loop
                newArr[k][l] = sum(sum(independ[k][l] * TimeDepend(x*10))) #New array calculates velocity at every point on the grid
        
        velLapse[x] = newArr #Saves new velocity
        distTravel += newArr * delta_T #Updates distance velocity
        distLapse[x] = distTravel #Saves new distance
        timepoint = x * delta_T * 1e6 #Time point currently being evaluated (in us)
        ax.cla() #Clear axis
        #Plots everything
        ax.set_title("Fluids position at %0.2f us" % timepoint)
        ax.set_xlabel("Width [um]")
        ax.set_zlabel("Height [um]")
        ax.set_ylabel("Position across the tube [mm]")
        ax.plot_surface(y, distLapse[x],z)
        plt.pause(1) #Wait a second before moving to the new plot, can adjust it if you'd like to wait longer

# Real Time Velocity
        
def VelRealTime():

    delta_T = 1e-5 #Time step [s]

    # Arrays that need to be initilized
    z = np.linspace(0, W, n1) # Width axis array
    y = np.linspace(0, H, n1) # Height axis array
    
    y, z = np.meshgrid(y, z) # Idk some fancy to to mesh them together??? who knows

    vx = np.zeros((n1, n1, 11, 11)) # Reinitializing array for vx as a 4D array

    perm = list(product(range(11),repeat=2)) #11x11 matrix of terms to add up

    const = -(16 / (poi * math.pi**2)) * (delta_P / delta_L) #Constant term (Depends only on constants input by the user)

    #Time Independent part of the equation, only needs to calculate it once
    def TimeInd():
        for y1 in range(n1):
            for z1 in range(n1):
                temp = np.zeros((11,11))
            
                for n in perm:
                    temp[n] = 0.1 * math.sin((2 * n[0] + 1) * math.pi * (y[y1][y1] / W)) * math.sin((2 * n[1] + 1) * math.pi * (z[z1][z1] / H)) * \
                            (((2 * n[0] + 1) * (2 * n[1] + 1) * (((((2 * n[0] + 1) * math.pi) / W)**2) + ((((2 * n[1] + 1) * math.pi) / H)**2))) ** -1)
                
                vx[y1][z1] = const * temp # Assigns value in vx array
        return vx;

    #Time dependent part of the equation, needs to be updated every time point
    def TimeDepend(t):
        upd = 1 - np.exp(part * t)
        return upd;

    part2 = np.zeros((11,11))

    #Summation of parts in the time independent part so it doesn't need to calculate it more than once (Might be overcomplicating it idk)
    def Expn():
        for k in perm:
            part2[k] = -(((((2 * k[0] + 1) * math.pi) / W)**2) + ((((2 * k[1] + 1) * math.pi) / H)**2)) * (poi / p)
        return part2;

    part = Expn()

    independ = TimeInd()

    newArr = np.zeros((n1,n1))
    distTravel = np.zeros((n1,n1))

    timePoints = 1001 #Number of time points to model


    distLapse = np.zeros((timePoints,n1,n1)) #Array holding distance travelled at every time point
    velLapse = np.zeros((timePoints,n1,n1)) #Array holding velocity at every time point

    #Just creates the figure
    fig1 = plt.figure(figsize=(14,10))
    ax = plt.axes(projection = "3d")
    
    for x in range(timePoints):
        for k in range(n1): #y loop
            for l in range(n1): #z loop
                newArr[k][l] = sum(sum(independ[k][l] * TimeDepend(x*10))) #New array calculates velocity at every point on the grid
        
        velLapse[x] = newArr #Saves new velocity
        distTravel += newArr * delta_T #Updates distance velocity
        distLapse[x] = distTravel #Saves new distance
        timepoint = x * delta_T * 1e6 #Time point currently being evaluated (in us)
        ax.cla() #Clear axis
        #Plots everything
        ax.set_title("Fluid velocity at %0.2f us" % timepoint)
        ax.set_xlabel("Width [um]")
        ax.set_zlabel("Height [um]")
        ax.set_ylabel("Velocity across the tube [mm/s]")
        ax.plot_surface(y, velLapse[x],z)
        plt.pause(1) #Wait a second before moving to the new plot, can adjust it if you'd like to wait longer


# RPA model (I'll add RT to this eventually). Also need to find new constants :(
        
def RPA():
    
    #Find whole paper at: https://www.sciencedirect.com/science/article/pii/S1369703X16301152#bib0065
    
    # Sources of experimental data to potentially try (will do this when I'm not big depresso):
    # https://www.sciencedirect.com/science/article/pii/S095869461930216X
    # https://link.springer.com/article/10.1007/s12161-017-0820-7
    # https://link.springer.com/article/10.1007/s00604-017-2144-0

    #Reaction rate constants (Check paper for units)
    k_1 = 1e5
    k_1f = 5e8
    k_1r = 5e3
    k_eq2a = 68
    k_eq2af = 1e8
    k_eq2ar = 1.471
    k_2b = 47e-3
    k_eq2c = 3e6
    k_eq2cf = 1e8
    k_eq2cr = 331
    k_2d = 4.6e-3
    KM_3a = 20.35e-6
    k_3f = 1e8
    k_3r = 59.37
    k_4acat = 4.22
    k_4bcat = 8.32
    k_5f = 1.2e7
    k_5r = 0.06
    k_6f = 87
    k_7 = 4.1
    k_8 = 1.13

    m = 4 #Number of gp32 binding sites on a primer
    n = 8 #Number of recombinase binding sites on a primer
    bp = 32 #Number of base pairs in primer
    B = 116 #Number of base pairs in template DNA

    G = 9.4e-7 #Concentration of gp32 protein [M]
    F = 4.8e-7 #Forward primer concentration [M]

    #Initial concentrations for reactants [M]
    R0 = 5.9e-6 # Recombinase complex
    FG0 = 0 # Forward primer/Gp32 complex
    FGR_prime0 = 0 # Forward primer/Gp32 complex intermediate
    FGR0 = 0 # Stable forward primer/Gp32 complex
    FGnR_prime0 = 0 # Forward primer/Gp32 complex w/ n recombinase molecules
    FnR0 = 0 # Forward primer complexed with stable filament of n recombinase molecules
    FnRD0 = 0 # Forward primer complexed with unstable filament of n recombinase molecules and DNA
    D0 = 0.1 #Concentration in ng/uL specified on the paper
    FD0 = 0 # Forward primer/DNA complex
    P0 = 1.34e-6 # Polymerase
    PFD0 = 0 # Polymerase/forward primer/DNA complex
    dNTPs0 = 2e-4 # Deoxynucleotide triphosphate
    PPi0 = 0 # Inorganic phosphate

    y0 = [R0,FG0,FGR_prime0,FGR0,FGnR_prime0,FnR0,FnRD0,D0,FD0,P0,PFD0,dNTPs0,PPi0] #Initial conditions array

    t = np.linspace(0,1200,1000000) #Time array

    def model(y,t):
    
        R = y[0]
        FG = y[1]
        FGR_prime = y[2]
        FGR = y[3]
        FGnR_prime = y[4]
        FnR = y[5]
        FnRD = y[6]
        D = y[7]
        FD = y[8]
        P = y[9]
        PFD = y[10]
        dNTPs = y[11]
        PPi = y[12]
    
        dRdt = - 2 * k_eq2af * R * n * FG + 2 * k_eq2ar * FGR_prime - 2 * k_eq2cf * R * FGR_prime + 2 * k_eq2cr * FGnR_prime + 2 * k_4acat * FnRD + 2 * k_4bcat * FnRD
        
        dFGdt = (1 / m) * k_1f * G * F - (1 / m) * k_1r * FG - k_eq2af * R * n * FG + (1 / n) * k_eq2ar * FGR_prime
    
        dFGR_primedt = k_eq2af * R * (n * FG) - (1 / n) * k_eq2ar * FGR_prime - k_2b * FGR_prime
    
        dFGRdt = k_2b * FGR_prime - (1 / (n-1)) * k_eq2cf * R * FGR + (1 / (n-1)) * k_eq2cr * FGnR_prime
    
        dFGnR_primedt = (1 / (n-1)) * k_eq2cf * R * FGR - (1 / (n-1)) * k_eq2cr * FGnR_prime - k_2d * FGnR_prime
    
        dFnRdt = k_2d * FGnR_prime + (1 / n) * k_3r * FnRD - (1 / n) * k_3f * FnR * D
    
        dFnRDdt = (1 / n) * k_3f * FnR * D - (1 / n) * k_3r * FnRD - (1 / n) * k_4acat * FnRD - (1 / n) * k_4bcat * FnRD
    
        dDdt = - (1 / n) + k_3f * FnR * D + (1 / n) * k_3r * FnRD + 2 * k_6f * PFD + 2 * k_7 * PPi * FD - 2 * k_8 * PPi * D
        
        dFDdt = (1 / n) * k_4acat * FnRD + (1 / n) * k_4bcat * FnRD - k_5f * P * FD + k_5r * PFD - k_7 * PPi * FD + k_8 * PPi * D
    
        dPdt = - k_5f * P * FD + k_5r * PFD + k_6f * PFD
    
        dPFDdt = k_5f * P * FD - k_5r * PFD - k_6f * PFD
        
        ddNTPsdt = - 2 * B * k_6f * PFD + 2 * B * k_8 * PPi * D + 2 * bp * k_7 * PPi * FD
    
        dPPidt = 2 * k_4acat * FnRD + 2 * B * k_6f * PFD - 2 * B * k_8 * PPi * D - 2 * bp * k_7 * PPi * FD
    
        return [dRdt, dFGdt, dFGR_primedt, dFGRdt, dFGnR_primedt, dFnRdt, dFnRDdt, dDdt, dFDdt, dPdt, dPFDdt, ddNTPsdt, dPPidt]

    soln = odeint(model,y0,t) #Solves differential equations using Euler method

    DNA = soln[:,7]

    #Plots DNA concentration over time
    plt.figure(1)
    plt.plot(t,DNA)
    plt.xlabel('Time (s)')
    plt.ylabel('DNA Concentration [M]')
    plt.title('DNA Concentration vs time')

# Now, functions you can run:

# Pouiselle Flow
#PFlow()

# Volumetric Flow Rate
#volFlow()

# Real Time Position
#PosRealTime()

# RPA
#RPA()
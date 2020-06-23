"""
Created on Sun Jun 21 16:52:15 2020

@author: Matthew Chan
"""

import math

# Constants assumming water
poi = 1 # Poise mPa * s
p = 1 # density g/cm^3

# Dimensions of the dumb thing
W = 500 # width in um   .5mm
H = 700 # height in um  .7mm

# Chnage of pressure and length thingys that I have no clue what these are
delta_P = -20  #change in pressure in
delta_L = 200 # length of path on microfluidic device in mm

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
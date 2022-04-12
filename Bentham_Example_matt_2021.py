# Christopher Joel Williams 17/08/2016
# Update by Matthew Young 05/10/21
# Import the necessary libraries
from ctypes import *
import time
import pyvisa as visa
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# Set up the VISA rescourse for mulitimeter

# Call the DLL and create pointers in mutable memory space to the config and attribute files designed to build a specified system model

# Build the specified system model and assign its attributes

# Initialise the Bentham active group and park it ready to perform tasks

# pull parameters from file crerated by the app

# set up VISA rescources for multi meter

# obtain reading from multimeter from buffer






# Set up the VISA rescourse for mulitimeter
rm = visa.ResourceManager()
print(rm.list_resources())
my_resource = rm.open_resource(
    'GPIB0::1::INSTR', read_termination='\r', send_end=True)


# Call the DLL and create pointers in mutable memory space to the config and attribute files designed to build a specified
# system model
mydll = CDLL("C:/Users/Public/Documents/UoN_SDK\SDK\Benhw\Win64/benhw64.dll")
CFGF = create_string_buffer(
    b'C:\Users\Public\Documents\UoN_SDK\SDK\Configuration Files\system.cfg', size=None)
errorrep1 = create_string_buffer(8)
ATT = create_string_buffer(
    b'C:\Users\Public\Documents\UoN_SDK\SDK\Configuration Files\system.atr', size=None)

# Build the specified system model and assign its attributes
result = mydll.BI_build_system_model(CFGF, errorrep1)
result2 = mydll.BI_load_setup(ATT)
# Initialise the Bentham active group and park it ready to perform tasks
mydll.BI_initialise()
park = mydll.BI_park()
print(park)
# Check that the system model has been built successfully
print(result)
print(repr(errorrep1.value))

# Define some important starting values
STARTWAVELENGTH = 300
STOPWAVELENGTH = 500
WAVESTEP = 2
SETTLE = 1000
# two values are set up so that one can be iterated in the look and the other
# can be passed to the C function BI_select_wavelength
PI = POINTER(c_uint)
settle = c_longdouble(1000)
k = c_double(STARTWAVELENGTH)
k1 = STARTWAVELENGTH
# Empty lists used to pass data to
K = []
s = []

# Conditional set up to check starting parameters are correct
if STARTWAVELENGTH < STOPWAVELENGTH and WAVESTEP > 0:
    # Loop to iterate over the desired range of wavelengths
    while k1 < STOPWAVELENGTH:
        # Sets the stepper motor wavelength
        result2 = mydll.BI_select_wavelength(k, byref(settle))
        # If setting the wavelength is successful the counter and wavelength variables are increased by the wavestep.
        if result2 == 0:
            k1 += WAVESTEP
            k = c_double(k1)
            SETTLE = settle.value
            # Assign the wavelengths in each iteration to a list
            K.append(k1)
            time.sleep(SETTLE/1000.0)
            # Get data from the lock in amplifier and store it in a list
            s.append(my_resource.query('OUTP? 3'))
# Translate list of data from lock in amplifier into a numpy array

length = len(s)
data = np.zeros(length)
data[0:length] = s[0:length]

# Tranform list of wavelengths into a numpy array
K = np.array(K)
# Adjust existing variables
#TrueWavelength = K 
TrueWavelength = K - 0.0098*K + 17.44

e = 1.6021773349e-19   
c = 299792458         
h = 6.6260755404e-34   

energy = ((h*c)/(TrueWavelength))/e
plt.plot(energy, data)
#plt.plot(TrueWavelength, data)
current = data/(10*(10**6))


# Save data to txt file
now = datetime.now()
dt_string = now.strftime("%d.%m.%Y %H.%M:%S")
#path = "B%test" + str(dt_string) + ".txt"
path = r"C:\Users\Labuser\Desktop\sampleData" + str(dt_string) + ".txt"

with open(path, "w") as txt_file:
    i=0
    while i<length:
        txt_file.write(str(data[i]) + " " + str(energy[i]) + "\n")
        i+1
        

# Close the Lock in amplifier followed by the resource manager
my_resource.write('REST')
my_resource.clear()
my_resource.close()
rm.close()
plt.show()

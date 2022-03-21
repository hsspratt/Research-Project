# %%
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

def nearestRefraction(x_Value_Store, y_Value_Store, Single_x_Value):
        
    x_diffs = Single_x_Value-x_Value_Store

    if np.where(x_diffs==0)[0].size > 0:
        lowest=np.where(x_diffs==0)[0]
    elif np.where(x_diffs==0)[0].size <= 0:
        lowest=max(np.where(x_diffs>0)[0])

    highest=lowest+1

    y_diff = y_Value_Store[highest]- y_Value_Store[lowest]

    if y_diff != 0:
        Lambda_Percentage = x_diffs[lowest]/(x_Value_Store[highest]-x_Value_Store[lowest])
    elif y_diff == 0:
        Lambda_Percentage = 0

    if y_diff > 0:
        RefractionTrueValue = y_Value_Store[lowest] - (y_diff)* Lambda_Percentage
    elif y_diff < 0:
        RefractionTrueValue = y_Value_Store[lowest] + (y_diff)* Lambda_Percentage
    elif y_diff == 0:
        RefractionTrueValue = y_Value_Store[lowest]
    else:
        print('Error Alert')

    return RefractionTrueValue

your_path = '/Users/harold/Documents/Academia/Nottingham Uni/Year 4/Research Project/Report/Coding/Data/OceanOpticsData/Harry/'
files = sorted(os.listdir(your_path))

data = np.zeros((1,2))

for file in files:
    if os.path.isfile(os.path.join(your_path, file)):
        array = np.loadtxt(your_path+str(file))

        max_value = np.max(array[:,1])
        max_index = np.argmax(array[:,1])

        data_temp = array[max_index,:]

        data = np.vstack((data,data_temp))


plt.figure('Intensity Graph')

plt.plot(data[6:,0], data[6:,1],label=r'Light Intensity')
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Photon Intensity')
#plt.legend((r'Theoretical Curve $ \\ \vspace{0.25cm} \alpha = A(hv-E_g)^{\frac{1}{2}}$', r'Experimental Results'),
#            shadow=False, loc=(0.38, 0.4), handlelength=1.5, fontsize=13)
            
#ax.set_aspect(1./ax.get_data_ratio())

plt.savefig("BulbLightIntensity.svg", format = 'svg', dpi=1200)

wavelength = np.linspace(data[6,0], data[-1,0],1000)
intensity = np.zeros(np.size(wavelength))


for i in range(len(wavelength)):
    print(i)
    intensity[i] = nearestRefraction(data[6:,0], data[6:,1], wavelength[i])

#plt.figure('Intensity')
#plt.plot()

plt.plot(wavelength,intensity)

 # %%


% name 
% pratical
% date

%% Setting up output

daq.reset
clear, close all;
clc;
s = daq.createSession('ni');
% Creates session object

s.addDigitalChannel('Devi','PortO/LineO:7','OutputOnly');
% Adds 8 digital output channels (numbered 0:7) on the DAQ card

%Undefined variable "daq? or class "daq.reset".
%Error in FINALCODE (line 6)
%daq.reset

%% Setting up input

s1 = daq.createSession('ni');
% Creates session object

s1.addAnaloglnputChannel('Dev1',0,'Voltage');
% Adds 1 analog input channel (numbered 0) on the DAQ card

s1.Rate = 1000; s1.NumberOfScans = 1000;
% Rate is the number of samples collected per second
% Number of scans is the total number of samples to be collected per
% trigger

intensity=[];
% Emptry matrix to be filled in as the code runs from the loop

%% Moving motor

nout = [51 102 204 153];
% This is the decimal sequence for forward motion

borf=input('Do you want to go backwards or forwards b/f ?','s');
% Obtaining information from the user in terms of string and stores the
% result in the variable 'borf'

number_wavelength=input('How many wavelengths would you like to move?')
% Obtaining the number of wavelength to be moved and stores the
% result in the variable 'number_wavelength'

% If the user inputted 'b', the wavelength will move backwards
% whereas 'f', the wavelength will move forward
% This loops changes the wavelengths in steps while obtaining the data

if borf == 'b'
    count=4;
% As it is going backward/ start from the last (4th) element of the
% decimal sequence of the stepper motor and work backwards

for n = 1 : round(number_wavelength/0.6584);
% This converts the number of wavelength into data index in order
% for the motor to change the number of wavelength
% on the monochromator

outputSingleScan(s,dec2binvec(nout(count),8));
% This will output the m bit binary representation converted
% from integer #number# to command the stepper motor to
% alter the wavelength of the monochromator

[data,timestamps] = startForeground(s1);
intensity(n)= mean(data);
% This inputs data into the empty matrix 'intensity'

plot(n,intensity(n),'*');

% This will show the graph of intensity against data index being
% plotted every iteration of the loop until the loop is
% terminated

hold on;
count=count-1;
% This is so that the decimal sequence work backwards every
% iteration

if count==0
    count=4;
end

% When the count decreases to 0, reset it to 4 again so that the
% stepper motor keeps moving backward
        pause(0.05);
    end
end

if borf == 'f'

count=1;

% As it is going forward, start from the first element of the
% decimal sequence of the stepper motor and work forward

for n=1 : round(number_wavelength/0.6584)
% This converts the number of wavelength into data index in order
% for the motor to change the number of wavelength
% on the monochromator

outputSingleScan(s,dec2binvec(nout(count),8)) ;
% This will output the m bit binary representation converted
% from integer #number# to command the stepper motor to
% alter the wavelength of the monochromator

[data,timestamps] = startForeground(s1);
intensity(n) = mean(data);
% This inputs data into the empty matrix 'intensity'

plot(n,intensity(n),'*');
% This will show the graph of intensity against data index being
% plotted every iteration of the loop until the loop is
% terminated

hold on;
count=count+1;

% This is so that the decimal sequence work forward every
% iteration

if count==5
    count=1;
end

% When the count increases to 5, reset it to 1 again so thatthe
% stepper motor keeps moving forward

        pause(0.05);
    end
end

%% Analysis of direct band gap sample

clear all; close all; clc;
% Load saved workspace data
I_dwl=dlmread('datadirectwl.txt').';
I_dw2=dlmread('datadirectw2.txt').';
I_dwol=dlmread('datadirectwol.txt').';
I_dwo2=dlmread('datadirectwo2.txt').';

% Calculates the mean of the intensity with and without the sample
% recorded by the detector

Intensity_l= (I_dwl+I_dw2)./2;
Intensity_0= (I_dwol+I_dwo2)./2;

% The wavelength over which the measurements were taken
wavelength_direct= (750:0.67567568567:950).*10.^-9;

T_1= Intensity_1./Intensity_0;
% Equation to calculate the tranmission coefficient
T_1min= min(T_1)
% Finding the minimum of the transmission coefficient
T_2= T_1+(-T_1min);

% Adding the difference between 0 and the minimum of transmission
% coefficient to the rest of the points so there will be no negative
% y value

plot((wavelength_direct), T_2);
% This plots the transmission coefficient against the wavelength

xlabel('Wavelength (m)')
ylabel('Transmission coefficient')
title('Graph of transmission coefficient against wavelength for GaAs')
% These labels the plot

% Constants

n_GaAs= 3.5227; %https://refractiveindex.info/?shelf=main&book=GaAs&page=Skauli
% The refractive index of GaAs

R_GaAs= [(n_GaAs-1)./(n_GaAs+1)].^2;
% Equation to calculate R

x_GaAs= 0.32.*10.^-3 ;
error_GaAs= 0.01.*10.^-3 ;
% The thickness of the GaAs sample and it's error

a=(1-R_GaAs).^4;
b=4.*(T_2.^2).*((R_GaAs).^2);
c=(a+b).^0.5;
d=(1-R_GaAs).^2;
numerator=c-d;
denominator=2.*T_2.*(R_GaAs.^2);
logterm=numerator./denominator;
Alpha_GaAs=(-x_GaAs.^-1).*log(logterm);
% Equation to calculate alpha of GaAs

h=6.62607004.*10.^-34;
% Planck's constant

c=3.*10.^8;
% Speed of light

E_direct= ((h.*c)./(wavelength_direct))./(1.6.*10.^-19);
% Equation to calculate the optical energy

LogAlpha_GaAs=logl0(Alpha_GaAs);
% Calculating the log of Alpha_GaAs

Alpha_squared=(Alpha_GaAs).^2 ;
% Calculating the square of Alpha_GaAs cftool

% Using curve fit tool to plot graph and analyse the data
% Semi- log graph of absorption coefficient against optical energy
% and graph of square of absorption coefficient against optical energy
% was plotted.

pbs = predint(fittedmodel,E_direct,0.95);
plot(fittedmodel,E_direct,LogAlpha_GaAs);
hold on;
ylabel('The log of Absorption Coefficient, Log(\alpha)')
xlabel('Optical Energy, E (eV)')
title('Semi-log graph of Absorption Coefficient, \alpha against Optical Energy, E for GaAs')

%% Analysis of the indirect band gap sample

clear all; close all; clc;
% Load saved workspace data

I_idwl=dlmread('dataindirectwl.txt').';
I_idw2=dlmread('dataindirectw2.txt').';
I_idw3=dlmread('dataindirectw3.txt').';
I_idwol=dlmread('dataindirectwol.txt').';
I_idwo2=dlmread('dataindirectwo2.txt').';
I_idwo3=dlmread('dataindirectwo3.txt').';


% Calculates the mean of the intensity with and without the sample
% recorded by the detector

Intensity_1 = (I_idw1+ I_idw2+ I_idw3)./3
Intensity_0 = (I__idwo1+ I_idwo2+ I_idwo3)./3;

% The wavelength over which the measurements were taken
wavelength_indirect- (400:0.65789473694: 700).*10.*-9;

T_l= Intensity_l./Intensity_O;
% Equation to calculate the tranmission coefficient
T_lmin= min(T_l)
% Finding the minimum of the transmission coefficient
T_2- T_l+(-T_lmin);
% Adding the difference between 0 and the minimum of transmission
% coefficient to the rest of the points so there will be no negative
% y value
plot((wavelengthindirect), T_2);
% This plots the transmission coefficient against the wavelength
xlabel( Wavelength (m)')
ylabel('Transmission coefficient')
title('Graph of transmission coefficient against wavelength for GaP')

Constants

n_GaP = 3.1870; %http://refractiveindex.info/?shelf=main&book=GaP&page=Jellison




plot(E_direct,pbs,'m?') ;
xlim([1.25 1.5])
ylim([2 3])

% Semi log graph of absorption coefficient against optical energy with its
% line of best fit and 95% confidence bound are plotted and the graph
% labelled.
% The same is repeated for the graph of square of absorption coefficient
% against optical energy.


% The refractive index of GaP
R_GaP- [(n_GaP-l)./(n_GaP+l)]

% Equation to calculate R
x_GaP= 0.35.*10.*-3 ;
error_GaP= 0.01.*10.*-3;

% The thickness of the GaAs sample and it's error

a=(1-R_GaP).*4;
b=4.*(T_2.*2).*((R_GaP).^2);
c=(a+b).*0.5;
d=(1-R_GaP).^2;
numerator=c-d;
denominator=2.*T_2.*(R_GaP.*2);
logterm=numerator./denominator;
Alpha_GaP=(-x_GaP.*-l).*log(logterm);
% Equation to calculate alpha of GaP

h=6.62607004.*10.^-34;
% Planck's constant

c= 3.*10.^8;
% Speed of light

E_direct=((h.*c)./(wavelength_indirect))./(1.6.*10.^-19);
% Equation to calculate the optical energy

minGaP=min(Alpha_GaP);
% Finding the minimum value of Alpha_GaP

AlphaGaP=Alpha_GaP+(-minGaP);
% This is so that there will be no negative y- value

Alpha_root=(AlphaGaP).^0.5;
% This is to calculate the square root of AlphaGaP

% cftool
% Curve fit tool used to obtain information about the plot

Eindirect = linspace(1.6,3.2,456);

% This generates 456 points between 1.6 and 3.2

alpha_a = 28.29.*Eindirect+(-46.87);

% Calculates the value of alpha corresponding to the value of Eindirect
% generated in the previous part

Alpha_root2= 28.29.*E_indirect+(-46.87);
% This generates the value of square root of absorption coefficient
% obtainend from the fit plotted in cftool in the previous part

Alpharoot2= Alpha_root2;

% This creates a variable that has the value of Alpha_root2 so that if it
% is modified the value of Alpha_root2 will not be affected

Alpharoot2square= (Alpharoot2).^2;
% This removes the square root on the alpha

Alpharoot1 = Alpha_root;
% This creates a variable that has the value of Alpha_root so that ifit
% is modified the value of Alpha_root will not be affected

Alpharootlsquare = (Alpharootl).^2;
% This removes the square root on the alpha

actualdifference= Alpharootlsquare-Alpharoot2square.';
% This calculates the difference between the alpha without square rooton
% the alpha

finalroot= sqrt(actualdifference);
% This square root the difference between the alpha found in the previous
% part

ActualAlpharootfinal= finalroot;

% Renaming the previous part with ActualAlpharootfinal tcftool
% Curve fit tool is used to find the line of best fit for the graph of
% square root of the absorption coefficient against the optical energy and
% for the graph of square root of the difference between the alpha against
% optical energy

pbs = predint(fittedmodel, E_indirect,0.95);

plot(fittedmodel, E_indirect, ActualAlpharootfinal);
hold on;

set(0,'DefaultTextlnterpreter','tex')

ylabel('The difference between the square root of Absorption Coefficient, (\alpha-\alpha_a)*{1/2} (m*{-1/2})')

xlabel('Optical Energy, E (eV)')

title('Graph of the difference between the square root of Absorption Coefficient against Optical Energy for GaP')
plot(E_indirect,pbs,' m?');

% Graph of square root of the absorption coefficient against theoptical
% energy and the graph of square root of the difference between the
% alpha against optical energy together with their line of best fit and 95%
% confidence bound are plotted and the graph labelled.





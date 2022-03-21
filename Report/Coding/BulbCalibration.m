%% Bulb Calibration

% clear all previous variables
clear all
clc

% imports necissary functions
import OpticalAnalysisFunctions.nearestValue

% figures no display in LaTeX
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% sort through files

%myDir = uigetdir; %gets directory
%myDir = '/Users/harold/Documents/Academia/Nottingham Uni/Year 4/Research Project/Report/Coding/Data/OceanOpticsData/Harry/';
%myDir = '/Users/harold/Documents/Academia/Nottingham Uni/Year 4/Research Project/Report/Coding/Data/OceanOpticsData/08022022/';
myDir = '/Users/harold/Documents/Academia/Nottingham Uni/Year 4/Research Project/Report/Coding/Data/OceanOpticsData/13042022/'
myFiles = dir(fullfile(myDir,'*.txt')); %gets all wav files in struct
experimental_intensitydata = zeros(1,2);

for k = 1:length(myFiles)
    baseFileName = myFiles(k).name;
    fullFileName = fullfile(myDir, baseFileName);

    fprintf('\b\b\b\b\b\b\b\b\b\b\b');
    fprintf(baseFileName)
    num = importdata(fullFileName);
    
    [max_value, max_index] = max(num(:,2));
    data_temp = num(max_index,:);
    experimental_intensitydata = vertcat(experimental_intensitydata,data_temp);
end

fprintf('\n')

%% calculating intensity spectra

% look at intensity data and establish where to start and end from
% this is becuase on the very small values the peaks are in high lambdas

clean_intensitydata = experimental_intensitydata(11:end,:);

% datasize = length(clean_intensitydata);

wavelengthFileName = '/Users/harold/Documents/Academia/Nottingham Uni/Year 4/Research Project/Report/Coding/1st Resistor - 20ms - 4 readings/_17.02.22 - 16.38.22.txt'
num = importdata(wavelengthFileName);
wavelengths = num(:,2);
voltages = num(:,1);
% wavelength = linspace(min(intensitydata(:,1)), max(intensitydata(:,1)), datasize);

wavelength_MAX = max(clean_intensitydata(:,1))
wavelength_MIN = min(clean_intensitydata(:,1))

% wavelength_MAXINDEX = find(wavelengths==wavelength_MAX)
% wavelength_MININDEX = find(wavelengths==wavelength_MIN)

correct_range_wavelength = find((wavelengths>wavelength_MIN) & (wavelengths<wavelength_MAX))
correct_wavelength = wavelengths(min(correct_range_wavelength):max(correct_range_wavelength))

% n = 5;
clean_intensitytest = zeros(size(correct_wavelength));

for i=1:max(size(correct_wavelength))
    clean_intensitytest(i) = OpticalAnalysisFunctions.nearestValue(experimental_intensitydata(:,1), correct_wavelength(i,1), experimental_intensitydata(:,2));
    i
end

% wavelength = wavelength(1:end-n);
correct_voltage = voltages(min(correct_range_wavelength):max(correct_range_wavelength));
normalised_intensity = clean_intensitytest.*(max(clean_intensitytest)./clean_intensitytest);
normalised_voltage = correct_voltage.*(max(clean_intensitytest)./clean_intensitytest);

%% plotting

figure('name','Bulb Intensity')

plot(clean_intensitydata(:,1),clean_intensitydata(:,2))
hold on
plot(correct_wavelength(1:10:end),clean_intensitytest(1:10:end),'*')
hold on 
plot(correct_wavelength,normalised_intensity)
xlabel('Wavelength (nm)')
ylabel('Photon Intensity')
legend('Experimental Data','Interpolated Data','Normalised Data',...
    'Location','best')

figure('name','Voltage Normalised')
plot(correct_wavelength,normalised_voltage)

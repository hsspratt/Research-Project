%% Henini Automated 

% figures no display in LaTeX
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% imports modules

import StaightLineFit.*

%% opening and running of app - commented

% app = AutomatedApp;
% 
% try
%     while app.RunScriptButton.Value == 0
%     pause(0.1)
%     end
% catch
%     warning('Invalid or deleted object.')
% end
% 
% disp('Program execution resumed')

%% import parameters data

parameters_file = '/Users/harold/Documents/Parameters.csv';
parameters_data = readtable(parameters_file);
parameters_values = parameters_data.Var2;

s_lambda = parameters_values(1);
f_lambda = parameters_values(2);
Step = parameters_values(3);
Readings = parameters_values(4);
t_eV = parameters_values(5);
t_lambda = parameters_values(6);
t_calculation = parameters_values(7);

%% run py file - commented

% py_file = '/Users/harold/testmatlab.py';
% pyrunfile(py_file)

%% Refractive Index 

import OpticalAnalysisFunctions.nearestValue
import OpticalAnalysisFunctions.CalculateRefractiveIndex

% loads the refractive index info

wavelength = linspace(350,1089,740)'

type = 'GaAs'

CalculateRefractiveIndex(wavelength, type)


%% Experimental Data

% Data for volatage without sample -- (V0_As)
% Data for volatage with sample    -- (V1_As)
% Wavelength data for volatages    -- (wavelength_As)



%% import experimental data

% GaAs_Data = '/Users/harold/Library/CloudStorage/OneDrive-TheUniversityofNottingham/OceanOpticsData/Automated Data/_21.03.22  13.00.57  20 ms 40/21.03.22  13.00.57  .txt';
% num = importdata(GaAs_Data);
% 
% wavelengths = num(:,2);
% std_data = num(:,3);

file_0 = "/Users/harold/Library/CloudStorage/OneDrive-TheUniversityofNottingham/OceanOpticsData/Automated Data/09.05.22  13.14.01  20 ms 15.0/Transmission_Data.txt"
num_0 = importdata(file_0);

voltages_0 = num_0(1:1089-349,1);
wavelengths_0 = num_0(1:1089-349,2);

file_1 = "/Users/harold/Library/CloudStorage/OneDrive-TheUniversityofNottingham/OceanOpticsData/Automated Data/09.05.22  15.30.37  20 ms 15.0/TransmissionData.txt"
num_1 = importdata(file_1);

voltages_1 = num_1(:,1);
wavelengths_1 = num_1(:,2);

% correct for systematic errors in wavelengths

import OpticalAnalysisFunctions.WavelengthsSystematicCorrection

wavelengths_0 = WavelengthsSystematicCorrection(wavelengths_0);
wavelengths_1 = WavelengthsSystematicCorrection(wavelengths_1);

% plot(wavelengths_0, voltages_0-min(voltages_0))
% hold on
% plot(wavelengths_1,voltages_1-min(voltages_1))

voltages_0 = voltages_0-min(voltages_0)
voltages_1 = voltages_1-min(voltages_1)

% correct for voltages

voltages = smooth(num(:,1));           % smooth
V_min = min(voltages).*0.9999;
voltages = voltages - abs(V_min);      % nomalise voltages to min 0

%% constants for experiment



T = smooth(voltages_1)./smooth(voltages_0);


h = 6.62607004*10^(-34); % planks constant
c = 299792458; % speed of light
x = 0.417.*10.^(-3);                               % thickness of sample
x_err = 0.001;                                     % error on thickness of sample

energy_ev = h*c./(wavelength.*10^(-9)*1.6*10^-19); % calculates the energy in eV for the respective wavelengths

R_1 = R_As
R_2 = (3.5860 - 1)^2 / (3.5860 + 1)^2
R_3 = 0.33

alpha_1 = -(x.^(-1)).*log((((1 - R_1).^4 + 4.*(T.^2).*(R_1.^2)).^0.5 - (1 - R_1).^2)./(2.*T.*(R_1.^2)));
alpha_2 = -(x.^(-1)).*log((((1 - R_2).^4 + 4.*(T.^2).*(R_2.^2)).^0.5 - (1 - R_2).^2)./(2.*T.*(R_2.^2)));
alpha_3 = -(x.^(-1)).*log((((1 - R_3).^4 + 4.*(T.^2).*(R_3.^2)).^0.5 - (1 - R_3).^2)./(2.*T.*(R_3.^2)));

plot(energy_ev, alpha_1)
hold on
plot(energy_ev, alpha_2)
hold on 
plot(energy_ev, alpha_3)


%offsetalpha = alpha + abs(-500.2303);

squarealpha = alpha.^2;

square_alpha = smooth(squarealpha);

%sqrtalpha = sqrt(abs(alpha));

%sqrtalpha = smooth(sqrtalpha);

alpha_check =  (1/x).*log(((1-R).^2)./T);

a = alpha_check - alpha;
plot(wavelength,a)

plot(eV_energy,alpha_check)

Joules_energy = (h*c)./(wavelength.*10.^(-9));
eV_energy = Joules_energy./(1.602176634*10^(-19));

Log_alpha_GaAs = log10(offsetalpha)

plot(eV_energy,Log_alpha_GaAs)

plot(wavelength,V0abs)
hold on
plot(wavelength,V1abs,'r')
hold off


figure(1)
plot(wavelength,T,'b*')
xlabel('Wavelength (nm)','Interpreter','latex')
ylabel('Transmission Coefficient','Interpreter','latex')
hold on
% plot(fittedmodel)
% plot(y1,x1,'b-')

hold off

figure(2)
plot(wavelength,V1abs,'b-')
hold on
plot(wavelength,V0abs,'r-')
xlabel('Wavelength (nm)','Interpreter','latex')
ylabel('Voltage','Interpreter','latex')

hold off

figure(3)
plot(wavelength,alpha)

hold off

figure(4)
plot(wavelength,square_alpha)

hold off

figure(5)
plot(eV_energy,alpha)
set(gca, 'YScale', 'log')

hold off

figure(6)
plot(eV_energy,square_alpha)


g = zeros(size(wavelengths));

%% conditional plotting - inital is no

if t_eV == 1

    figure('Name', 'Transmission Graph $(eV)$');
    plot(wavelengths, voltages)
end

if t_lambda == 1

    figure('Name', 'Transmission Graph $(\lambda)$');
    plot(energy_ev, voltages)
end


%% calculation - initial is yes

if t_calculation == 1

    for i=1:size(wavelengths)-1
        g(i) = (voltages(i+1)-voltages(i))/(wavelengths(i+1)-wavelengths(i));
    end
    
    g_cutoff = (max(voltages)-min(voltages))/(0.20*(max(wavelengths)-min(wavelengths)));
    
    g_straighline = g(find(g>g_cutoff));
    
    w_straightline = wavelengths(find(g>g_cutoff));
    v_straightline = voltages(find(g>g_cutoff));
    std_straighline = std_data(find(g>g_cutoff));
    
    fit = StraightLineFit(w_straightline, v_straightline, std_straighline, wavelengths, voltages);

end

%% extra jake

file = "/Users/harold/Library/CloudStorage/OneDrive-TheUniversityofNottingham/OceanOpticsData/Automated Data/09.05.22  15.30.37  20 ms 15.0/Transmission_Data.txt"
num = importdata(file);

voltages = num(:,1)
wavelengths = linspace(350,1089,740)'
plot(wavelengths, voltages)

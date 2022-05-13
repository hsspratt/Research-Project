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

R_As = CalculateRefractiveIndex(wavelength, type)


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

% voltages = smooth(num(:,1));           % smooth
% V_min = min(voltages).*0.9999;
% voltages = voltages - abs(V_min);      % nomalise voltages to min 0

%% constants for experiment

import OpticalAnalysisFunctions.CutExcessData
import OpticalAnalysisFunctions.DetectLongestStarightLine
import OpticalAnalysisFunctions.DetectStraightLine


T = smooth(voltages_1)./smooth(voltages_0);


h = 6.62607004*10^(-34); % planks constant
c = 299792458; % speed of light
x = 0.417.*10.^(-3);                               % thickness of sample
x_err = 0.001;                                     % error on thickness of sample

energy_ev = h*c./(wavelength.*10^(-9)*1.6*10^-19); % calculates the energy in eV for the respective wavelengths

R_1 = R_As;
R_2 = (3.5860 - 1)^2 / (3.5860 + 1)^2;
R_3 = 0.33;

alpha_1 = -(x.^(-1)).*log((((1 - R_1).^4 + 4.*(T.^2).*(R_1.^2)).^0.5 - (1 - R_1).^2)./(2.*T.*(R_1.^2)));
alpha_2 = -(x.^(-1)).*log((((1 - R_2).^4 + 4.*(T.^2).*(R_2.^2)).^0.5 - (1 - R_2).^2)./(2.*T.*(R_2.^2)));
alpha_3 = -(x.^(-1)).*log((((1 - R_3).^4 + 4.*(T.^2).*(R_3.^2)).^0.5 - (1 - R_3).^2)./(2.*T.*(R_3.^2)));

alpha_1 = alpha_1 - min(alpha_1);
alpha_2 = alpha_2 - min(alpha_2);
alpha_3 = alpha_3 - min(alpha_3);

plot(energy_ev, alpha_1)
hold on
plot(energy_ev, alpha_2)
hold on 
plot(energy_ev, alpha_3)



% max_alpha_index = find(max(alpha_3) == alpha_3);
% min_alpha_index = find(min(alpha_3) == alpha_3);
% 
% index_diff = (max_alpha_index + min_alpha_index)*0.5;
% 
% if rem(index_diff,2) ~= 0
%     index_diff = index_diff + 0.5;
% end
% 
% mid_alpha = (max(alpha_3)-min(alpha_3))/2;
% difference_points = abs(mid_alpha - alpha_3);
% index_mid_alpha = find(min(difference_points)==difference_points);
% 
% spliced_alpha = alpha_3(index_mid_alpha-index_diff:index_mid_alpha+index_diff);

gradient = zeros(1, max(size(energy_ev))-1);
gradient(1) = alpha_3(1)./(energy_ev(1));
for i=1:max(size(energy_ev))-1
    gradient(i) = (alpha_3(i+1) - alpha_3(i))./(energy_ev(i+1) - energy_ev(i));
end

squarealpha = smooth(alpha.^2);

sqrtalpha = smooth(sqrt(alpha));

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

%% Sigmoid Plotting

h = 6.63*10^-34;
c = 3*10^8;
%%

a = smooth(alpha_1);

% Find the min and max alpha values and what position they lie
amax = max(a);
amin = min(a);


nmax = find(amax==a);

% Splice excess data beyond 5% of the emmision peak
perc = 0.5;
nmaxrange = round(nmax-150);


l = length(a);
ai = a((l-nmaxrange):l,:);

% Finding the energy co-ordinate that is the midpoint between alpha max and
% alpha min
ahalf = (amax - amin)/2;
apos = abs(ahalf-ai);
aposmin = min(apos);
napos = find(aposmin==apos);

la = length(ai);
na = 1089 - la+1;
Ewave = (na:1089);
E = h*c./(Ewave.*10^(-9)*1.6*10^-19);
E0 = E(napos);


aimin = min(ai);

naimin = find(aimin == ai);
naimax = find(amax == ai);

ai = ai'

save('ExperimentalData.mat', 'amin', 'amax', 'E0', 'E', 'ai')

app = Sigmoid;


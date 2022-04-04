%% Henini Automated 

% figures no display in LaTeX
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% imports modules

import StaightLineFit.*

%% opening and running of app

app = AutomatedApp;

try
    while app.RunScriptButton.Value == 0
    pause(0.1)
    end
catch
    warning('Invalid or deleted object.')
end

disp('Program execution resumed')

%% import parameters data

% parameters_file = '/Users/harold/Documents/Academia/Nottingham Uni/Year 4/Research Project/Report/Coding/Parameters.csv';
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

%% run py file

py_file = '/Users/harold/testmatlab.py';
pyrunfile(py_file)

%% import experimental data

GaAs_Data = '/Users/harold/Library/CloudStorage/OneDrive-TheUniversityofNottingham/OceanOpticsData/Automated Data/_21.03.22  13.00.57  20 ms 40/21.03.22  13.00.57  .txt';
num = importdata(GaAs_Data);

wavelengths = num(:,2);
std_data = num(:,3);

% correct for systematic errors in wavelengths

m = -0.0091;
c = 18.1;
x = wavelengths;
y = x - (m*x + c);
wavelengths = y;

% correct for voltages

voltages = smooth(num(:,1));           % smooth
V_min = min(voltages).*0.9999;
voltages = voltages - abs(V_min);      % nomalise voltages to min 0

%% constants for experiment

h = 6.63*10^-34; % planks constant
c = 2.99*10^8; % speed of light
energy_ev = h*c./(wavelengths.*10^-9*1.6*10^-19); % calculates the energy in eV for the respective wavelengths

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


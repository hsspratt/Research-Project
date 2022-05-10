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
%pyrunfile(py_file)

%% Refractive Index 

import OpticalAnalysisFunctions.nearestValue

% loads the refractive index info

wavelength = linspace(350,1089,740)'

RefractiveIndex_GaAs = csvimport('/Users/harold/Documents/Academia/Nottingham Uni/Year 4/Research Project/Report/Coding/Data/Refractive/Papatryfonoset-2021-0.260-1.88-Ga_As.csv'); 
RefractiveIndex_GaAs(1,:)  = [];  % removing column titles

RefractiveIndexInfo_GaAs  = zeros(max(size(RefractiveIndex_GaAs)),min(size(RefractiveIndex_GaAs)));

for c=1:min(size(RefractiveIndex_GaAs))
    for r=1:max(size(RefractiveIndex_GaAs))
        RefractiveIndexInfo_GaAs(r,c)  = RefractiveIndex_GaAs{r,c};
    end
end

L_As  = RefractiveIndexInfo_GaAs(:,1).*1000;
N_As  = RefractiveIndexInfo_GaAs(:,2);
K_As  = RefractiveIndexInfo_GaAs(:,3);

n_As = zeros(size(wavelength));
k_As = zeros(size(wavelength));
R_As = zeros(size(wavelength));

for i=1:max(size(wavelength))
    % Approximated value of R calculated - Not calculating R manually with this method
    n_As(i) = OpticalAnalysisFunctions.nearestValue(L_As, wavelength(i), N_As);
    k_As(i) = OpticalAnalysisFunctions.nearestValue(L_As, wavelength(i), K_As);
    R_As(i) = ((n_As(i)-1)+k_As(i).^2)/((n_As(i)+1)+k_As(i).^2);
end

figure( 'Name', 'Refractive Index');

plot(L_As,N_As) % real refractive index - n
hold on
plot(L_As,K_As) % complex refracive index - ik
hold on
plot(wavelength(1:150:end),n_As(1:150:end),'*') % real refractive index - n
hold on
plot(wavelength(1:150:end),k_As(1:150:end),'*') % complex refracive index - ik

title('Refractive Index $200 -- 830 nm$','Interpreter','latex');
legend('Real refractive index','Complex refractive index', ...
    'Experimental Real refractive','Experimental Complex refractive','Interpreter','latex')

% legand needs changing

xlabel( 'Wavelength $/nm$', 'Interpreter', 'latex' );
ylabel( 'Refractive Index', 'Interpreter', 'Latex' );

%% Experimental Data

% Data for volatage without sample -- (V0_As)
% Data for volatage with sample    -- (V1_As)
% Wavelength data for volatages    -- (wavelength_As)



%% import experimental data

GaAs_Data = '/Users/harold/Library/CloudStorage/OneDrive-TheUniversityofNottingham/OceanOpticsData/Automated Data/_21.03.22  13.00.57  20 ms 40/21.03.22  13.00.57  .txt';
num = importdata(GaAs_Data);

wavelengths = num(:,2);
std_data = num(:,3);

% correct for systematic errors in wavelengths

m = -0.0091;
c = 18.1;
x = wavelengths;
y = x + (m*x + c);
wavelengths = y;

% correct for voltages

voltages = smooth(num(:,1));           % smooth
V_min = min(voltages).*0.9999;
voltages = voltages - abs(V_min);      % nomalise voltages to min 0

%% constants for experiment

file_sample = "/Users/harold/Library/CloudStorage/OneDrive-TheUniversityofNottingham/OceanOpticsData/Automated Data/09.05.22  15.30.37  20 ms 15.0/Transmission_Data.txt"
num_sample = importdata(file_sample);

voltages_sample = num_sample(:,1)
wavelengths_sample = linspace(350,1089,740)'

file_nothing = "/Users/harold/Library/CloudStorage/OneDrive-TheUniversityofNottingham/OceanOpticsData/Automated Data/09.05.22  13.14.01  20 ms 15.0/Transmission_Data.txt"
num_nothing = importdata(file_nothing);

voltages_nothing = num_nothing(1:1089-349,1)
wavelengths_nothing = num_nothing(1:1089-349,2)

T = smooth(voltages_sample)./smooth(voltages_nothing);


h = 6.62607004*10^(-34); % planks constant
c = 299792458; % speed of light
x = 0.417.*10.^(-3);                               % thickness of sample
x_err = 0.001;                                     % error on thickness of sample

energy_ev = h*c./(wavelength.*10^(-9)*1.6*10^-19); % calculates the energy in eV for the respective wavelengths

%R = R_As
R = (3.5860 - 1)^2 / (3.5860 + 1)^2
%R = 0.33

alpha = -(x.^(-1)).*log((((1 - R).^4 + 4.*(T.^2).*(R.^2)).^0.5 - (1 - R).^2)./(2.*T.*(R.^2)));
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

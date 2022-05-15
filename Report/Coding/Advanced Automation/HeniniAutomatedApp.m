%% Henini Automated 

% figures no display in LaTeX
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% imports modules

import StaightLineFit.*

%% opening and running of app - commented

AutomatedApp;

waitfor(AutomatedApp)

disp('Program execution resumed')

%% import parameters data

parameters_file = 'Parameters.csv';
parameters_data = readtable(parameters_file);
parameters_values = parameters_data.Var2;

s_lambda = parameters_values(1);
f_lambda = parameters_values(2);
Readings = parameters_values(3);
Step = parameters_values(4);
Transmission = parameters_values(5);
FastMethod = parameters_values(6);
Tauc = parameters_values(7);
Sigmoid = parameters_values(8);


%% run py file - commented

%py_file = '/Users/harold/testmatlab.py';
%pyrunfile(py_file)

%% Experimental Data

% Data for volatage without sample -- (voltages_0)
% Data for volatage with sample    -- (voltages_1)
% Wavelength data for volatages    -- (wavelengths_0 & wavelengths_1)

%% import experimental data

file_0 = "TransmissionData0.txt";
num_0 = importdata(file_0);

voltages_0 = num_0(:,1); %(1:1089-349,1);
wavelengths_0 = num_0(:,2); %(1:1089-349,2);

file_1 = "TransmissionData1.txt";
num_1 = importdata(file_1);

voltages_1 = num_1(:,1);
wavelengths_1 = num_1(:,2);
std_1 = num_1(:,3);

% check for same size array of wavelengths

if ~isequal(wavelengths_0, wavelengths_1)
    disp('Error - Wavelength arrays are not the same')
end

% correct for systematic errors in wavelengths

import OpticalAnalysisFunctions.WavelengthsSystematicCorrection

wavelengths_0 = WavelengthsSystematicCorrection(wavelengths_0);
wavelengths_1 = WavelengthsSystematicCorrection(wavelengths_1);

voltages_0 = voltages_0-min(voltages_0);
voltages_1 = voltages_1-min(voltages_1);

% smooth data

voltages_0 = smooth(voltages_0);
voltages_1 = smooth(voltages_1);

%% constants for experiment

import OpticalAnalysisFunctions.CutExcessData
import OpticalAnalysisFunctions.DetectLongestStarightLine
import OpticalAnalysisFunctions.DetectStraightLine

% calculates the transmission values over all wavelengths

T = smooth(voltages_1)./smooth(voltages_0);


h = 6.62607004*10^(-34);  % planks constant
c = 299792458;            % speed of light
x = 0.417.*10.^(-3);      % thickness of sample
x_err = 0.001;            % error on thickness of sample

energy_ev = h*c./(wavelengths_0.*10^(-9)*1.6*10^-19); % calculates the energy in eV for the respective wavelengths

%% Refractive Index 

import OpticalAnalysisFunctions.nearestValue
import OpticalAnalysisFunctions.CalculateRefractiveIndex

% loads the refractive index info

type = 'GaAs';

R_As = CalculateRefractiveIndex(wavelengths_0, type);

% Three main methods for caluclating R for this experiment 1 accurate and 
% two estimates. The more accurate method has been used but the estimates
% are there for comparison

R = R_As;                              
R_2 = (3.5860 - 1)^2 / (3.5860 + 1)^2; % R=((n-1)/(n+1))^2, n=3.187
R_3 = 0.33;                            % taken in linear region after bandgap

alpha = -(x.^(-1)).*log((((1 - R).^4 + 4.*(T.^2).*(R.^2)).^0.5 - (1 - R).^2)./(2.*T.*(R.^2)));
alpha_2 = -(x.^(-1)).*log((((1 - R_2).^4 + 4.*(T.^2).*(R_2.^2)).^0.5 - (1 - R_2).^2)./(2.*T.*(R_2.^2)));
alpha_3 = -(x.^(-1)).*log((((1 - R_3).^4 + 4.*(T.^2).*(R_3.^2)).^0.5 - (1 - R_3).^2)./(2.*T.*(R_3.^2)));

alpha = alpha - min(alpha);
alpha_2 = alpha_2 - min(alpha_2);
alpha_3 = alpha_3 - min(alpha_3);

squarealpha = smooth(alpha.^2);

sqrtalpha = smooth(sqrt(alpha));

Log_alpha_GaAs = log10(alpha);

save('TaucPlotting.mat', 'alpha', 'energy_ev', 'squarealpha', 'sqrtalpha', 'Log_alpha_GaAs')

%% conditional plotting - inital is no

if Transmission == 1

    figure('Name', 'Transmission Graph $(eV)$');
    plot(energy_ev, voltages_1)
end

%% calculation - initial is yes

if FastMethod == 1

    file_fastdata = 'TransmissionDataFastRun.txt';
    num_fastdata = importdata(file_fastdata);
    
    voltages_fastdata = WavelengthsSystematicCorrection(num_fastdata(:,1));
    voltages_fastdata = smooth(voltages_fastdata);
    voltages_fastdata = voltages_fastdata - min(voltages_fastdata);
    wavelengths_fastdata = num_fastdata(:,2);
    std_fastdata = num_fastdata(:,3);
    
    energy_ev = h*c./(wavelengths_fastdata.*10^(-9)*1.6*10^-19); % calculates the energy in eV for the respective wavelengths
    
    energy = energy_ev;
    voltages = voltages_fastdata;
    std_data = std_fastdata;

    g = zeros(size(energy));
   
    for i=1:size(energy)-1
        g(i) = (voltages(i+1)-voltages(i))/(energy(i+1)-energy(i));
    end
    
    g_cutoff = -(max(voltages)-min(voltages))/(0.2*(max(energy)-min(energy)));
    
    g_straighline = g(g<g_cutoff);
    
    w_straightline = energy(g<g_cutoff);
    v_straightline = voltages(g<g_cutoff);
    std_straighline = std_data(g<g_cutoff);
    
    fit = StraightLineFit(w_straightline, v_straightline, std_straighline, energy, voltages);

    MyCoeffs = coeffvalues(fit);
    BandGap = (-MyCoeffs(2))/(MyCoeffs(1));
    var_err = confint(fit);
    c_err = abs(var_err(1,2)-var_err(2,2))/2;
    m_err = abs(var_err(1,1)-var_err(2,1))/2;
    
    ComplexError = sqrt(((MyCoeffs(2)/MyCoeffs(1)^2)^2)*m_err^2 + ((1/MyCoeffs(1))^2)*c_err^2);
    SimpleError = BandGap*sqrt((c_err/MyCoeffs(2))^2 + (m_err/MyCoeffs(1))^2);

    disp(['Error = ' , num2str(ComplexError)])

end
%% Sigmoid Plotting

if Sigmoid == 1

    a = smooth(alpha);
    
    % Find the min and max alpha values and what position they lie
    amax = max(a);
    amin = min(a);
    
    %From here down 
    
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
    
    ai = ai';
    
    save('SigmoidData.mat', 'amin', 'amax', 'E0', 'E', 'ai')
    
    SigmoidApp;

end

%% Tauc Plotting

if Tauc == 1
    TaucPlotApp
end

%% Optical Abosrption MATLAB code - Analysis

import OpticalAnalysisFunctions.nearestRefraction
import OpticalAnalysisFunctions.DetectStraightLine
import OpticalAnalysisFunctions.LongestConsecutive
import OpticalAnalysisFunctions.DetectLongestStarightLine

%% loads all the data from the experiment

% load('Data/Experimental/GaAs - Harry/Combined_As.mat')
load('Data/Experimental/GaAs - Harry/V0_As.mat') 
load('Data/Experimental/GaAs - Harry/V1_As.mat')
load('Data/Experimental/GaAs - Harry/wavelength_As.mat')

%% Refractive Index 

% loads the refractive index info

RefractiveIndex_GaAs = csvimport('Data/Refractive/RefractiveIndexGaAs.csv'); 
RefractiveIndex_GaAs(1,:)  = [];  % removing column titles

RefractiveIndexInfo_GaAs  = zeros(46,3);

for c=1:3
    for r=1:46
        RefractiveIndexInfo_GaAs(r,c)  = RefractiveIndex_GaAs{r,c};
    end
end

L_P  = RefractiveIndexInfo_GaAs(:,1).*1000;
N_P  = RefractiveIndexInfo_GaAs(:,2);
K_P  = RefractiveIndexInfo_GaAs(:,3);

figure( 'Name', 'Refractive Index');

plot(L_P,N_P) % real refractive index - n
hold on
plot(L_P,K_P) % complex refracive index - ik

title('Refractive Index $200 -- 830 nm$','Interpreter','latex');
legend('Real refractive index','Complex refractive index','Interpreter','latex')

xlabel( 'Wavelength $/nm$', 'Interpreter', 'latex' );
ylabel( 'Refractive Index', 'Interpreter', 'Latex' );

hold off

%nearestRefraction(L, Lambda, N);

%% Experimental Data

% Data for volatage without sample -- (V0_As)
% Data for volatage with sample    -- (V1_As)
% Wavelength data for volatages    -- (wavelength_As)

V0 = smooth(V0); % smooths V0 data
V1 = smooth(V1); % smooths V1 data

V0 = V0*(-1); % volatages negative when they should be +ve
V1 = V1*(-1); % volatages negative when they should be +ve

T = V1./V0;
T_min = min(T);
T = T + abs(T_min);

% CAN'T FIND STRAIGHT LINE AS ITS A BIT QUESTIONABLE IN THAT REGION

% [x,y] = OpticalAnalysisFunctions.DetectLongestStarightLine(wavelength, T)
% plot(wavelength, T)
% hold on
% plot(x, y, '*')
%
% meanInsideInterval = mean(T(2070:2255)) 
% find average transmsission value on flat region at the end of the graph

meanInsideInterval = mean(T(2070:2255)) 

R = (1-meanInsideInterval)/(1+meanInsideInterval) % Value of R found via transmission graph
R = (3.5160 - 1)^2 / (3.5160 + 1)^2               % Approximate value of R to check

x = 0.417.*10.^(-3);                               % thickness of sample
x_err = 0.001;                                     % error on thickness of sample
h = 6.62607004*10^(-34);                           % plancks constant
c = 299792458;                                     % speed of light
Joules_energy = (h*c)./(wavelength.*10.^(-9));     % calculates energy in joules
eV_energy = Joules_energy./(1.602176634*10^(-19)); % converts joules to eV

% calculate each alpha for each wavelength
alpha_GaAs = -(x.^(-1)).*log((((1 - R).^4 + 4.*(T.^2).*(R.^2)).^0.5 - (1 - R).^2)./(2.*T.*(R.^2)));

offset_As = min(alpha_GaAs);               % finds alpha offset
offsetalpha = alpha_GaAs + abs(offset_As); % correct alpha offset so no negative values
square_alpha = offsetalpha.^2;             % finds alpha squared
square_alpha = smooth(square_alpha);       % smooths values
Log_alpha_GaAs = log10(offsetalpha);       %log alpha

plot(wavelength,V0abs)
hold on
plot(wavelength,V1abs,'r')
hold off
xlabel('Wavelength (nm)','Interpreter','latex')
ylabel('Volatage','Interpreter','latex')
% plots V1 and V0 against wavelength

plot(wavelength,T,'b')
xlabel('Wavelength (nm)','Interpreter','latex')
ylabel('Transmission Coefficient','Interpreter','latex')
% plots T against wavelength

plot(eV_energy,Log_alpha_GaAs)
xlabel('Optical Energy (eV)','Interpreter','latex')
ylabel('Logaritm of Absorption Coefficient','Interpreter','latex')
% plots log alpha vs energy

plot(eV_energy,square_alpha)
xlabel('Optical Energy (eV)','Interpreter','latex')
ylabel('Absorption Coefficient Squared','Interpreter','latex')
% plots alpha squared vs energy, the graph of interest

plot(eV_energy(1:40:end),offsetalpha(1:40:end),'*')
set(gca, 'YScale', 'log')
hold on
plot(eV_energy(1:30:1070),offsetalpha(1:30:1070),'r--')
legend('experimental values', 'theoretical values')
% plots theroetical vs experimental values

%% end 


figure('Name','Transmission Coefficient vs Optical Wavelength GaAs')

plot(Lambda,T,'*')

xlabel('Wavelength /nm', 'Interpreter', 'latex')
ylabel('Transmission Coefficient', 'Interpreter', 'latex')
%title({'Transmission Coefficient vs Optical Wavelength - Gallium Phosphide'})
hold off

% btw i dont think they like titles in matlab graphs i think they prefder a
% caption but just check in the formal report guide

n_P = zeros(size(Lambda));
k_P = zeros(size(Lambda));
R_P = zeros(size(Lambda));

for i=1:size(Lambda,2)
    % approximated value of R calculated
    % Are you not calculating R manually, or just using the estimate?

    n_P(i) = OpticalAnalysisFunctions.nearestRefraction(L_P, Lambda(i), N_P);
    k_P(i) = OpticalAnalysisFunctions.nearestRefraction(L_P, Lambda(i), K_P);
    R_P(i) = ((n_P(i)-1)+k_P(i).^2)/((n_P(i)+1)+k_P(i).^2);
    % R_P(i) = ((n_P(i)-1))/((n_P(i)+1));
end


x = 0.408*10^-3; % thickness of sample
xerr = 0.001*10^-3; % error on thickness
h = 6.63*10^-34; % planks constant
c = 2.99*10^8; % speed of light
energy_ev = h*c./(Lambda.*10^-9*1.6*10^-19); % calculates the energy in eV for the respective wavelengths


% for large wavelengths then extrapolate with curve fitting tool
% rootalpha  = (-x^(-1).*log([((1-R).^4+4.*T1.^2.*R.^2).^1/2-(1-R).^2]/2.*T1.*R.^2)).^1/2;
% for small wavelengths find curve
% not really to sure whats going on with the section above as you're
% subtracting alpha from itself so alphae is just 0 in this case

% Also check this equation I've changed it a lot one thing to remeber
% divition and multiplication involving arrays has to be used by using '.'
% if you forget to put the full stop then it doesn't divide individually

% Anyway i think this equation is right, its the same as mine now but just
% check anyway just in case

R = R_P;

rootalpha = -x^(-1).*log((((1 - R).^4 + 4.*(T.^2).*(R.^2)).^0.5 - (1 - R).^2)./(2.*T.*(R.^2)));

figure('Name','Root Alpha vs Energy (ev)')
plot(energy_ev,rootalpha,'*')
hold on
plot(Alda_low_energy_alpha)
hold on
pbs = predint(Alda_low_energy_alpha,energy_ev,0.68);
plot(energy_ev,pbs,'m--')
xlabel('Energy / eV', 'Interpreter', 'latex')
ylabel('$\sqrt{\alpha}$', 'Interpreter', 'latex')
legend('GaAs data points','Line of best fit','$68\%$ prediction bounds', 'Interpreter', 'latex')
hold off

% plots rootalpha vs energy with line of best fit in low energy region
% you can fiddle round with the limit to get a good looking graph

A = 3.728/1.722 % this is the value of Eg - Ep

alphaa = (energy_ev - A).^2;
% technically alpha a over c
% calculates value of alpha a needed for alpha e

alphae = (rootalpha - alphaa).^(0.5);
% calculates value of alpha e^(0.5)

figure
plot(energy_ev,alphae,'*')
pbs = predint(Alda_high_energies_alpha,energy_ev,0.68);
hold on
plot(Alda_high_energies_alpha)
hold on
plot(energy_ev,pbs,'m--')
xlabel('Energy / eV')
ylabel('Square root of alpha_e')
legend('off')
hold off
xlabel('Optical Energy (eV)')
ylabel('Absorption Coefficient of e to the half')
% plots alpha_e vs enegy to find Eg + Ep

B = 2042/931.8 % this is the value for Eg + Eg

bandGaAs = (A + B)/2 % band GaAs energy
phonon = (B - A)/2 % phonon energy

% you dont have errors for anything atm on the graph, you could do error bars
% if you want

% you dont have many points so thats maybe why the band GaAs is slightly off
% but all in all its very close

% on things to be wary of is that you cant propagate error to get the error
% on the x interecept using the error propagation equation, if you do
% theyll mark you down.

%rootalphae=rootalpha-rootalphaa;
% % for values with equal energy subtract curves from fittedmodels from each
% % other

%hold on
%plot(E,rootalphaa);
%plot(E,rootalpha);
% % x1 intercept= Eg-Ep CURVE FITTING
%hold off
% %% 5)
% plot(E,rootalphae);
% % x2 intercept= Eg+Ep CURVE FITTING
% %% 6)
% Eg=(x1+x2)./2;
% Ep=(x2-x1)./2;
%% Huge Test

[x,y] = OpticalAnalysisFunctions.DetectLongestStarightLine(Lambda, T)
plot(Lambda, T)
hold on
plot(x, y, '*')
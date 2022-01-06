%% Optical Abosrption MATLAB code - Analysis

import OpticalAnalysisFunctions.nearestRefraction
import OpticalAnalysisFunctions.DetectStraightLine
import OpticalAnalysisFunctions.LongestConsecutive
import OpticalAnalysisFunctions.DetectLongestStarightLine

% loads all the data from the experiment

load('Data/Experimental/GaP - Alda/I02.mat')
load('Data/Experimental/GaP - Alda/I3.mat')
load('Data/Experimental/GaP - Alda/lambda.mat')
load('Data/Experimental/GaP - Alda/low_energies_alpha.mat')
load('Data/Experimental/GaP - Alda/high_energies_alpha.mat')

% loads the refractive index info

RefractiveIndex_GaP = csvimport('Data/Refractive/RefractiveIndexGaP.csv'); 
RefractiveIndex_GaP(1,:)  = [];  % removing column titles

RefractiveIndexInfo_GaP  = zeros(46,3);

for c=1:3
    for r=1:46
        RefractiveIndexInfo_GaP(r,c)  = RefractiveIndex_GaP{r,c};
    end
end

L_P  = RefractiveIndexInfo_GaP(:,1).*1000;
N_P  = RefractiveIndexInfo_GaP(:,2);
K_P  = RefractiveIndexInfo_GaP(:,3);

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

T_0 = abs(I3./I02); % whole data set of T recorded
T = T_0(349:410); % confining T and lambda to regiojn of interest
Lambda=LambdaStore(349:410); % also making sure they are the same size



figure('Name','Transmission Coefficient vs Optical Wavelength GaP')

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
legend('GaP data points','Line of best fit','$68\%$ prediction bounds', 'Interpreter', 'latex')
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

bandgap = (A + B)/2 % band gap energy
phonon = (B - A)/2 % phonon energy

% you dont have errors for anything atm on the graph, you could do error bars
% if you want

% you dont have many points so thats maybe why the band gap is slightly off
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
%% Optical Abosrption MATLAB code
load('I02.mat')
load('I3.mat')
load('lambda.mat')
load('low_energies_alpha.mat')
load('high_energies_alpha.mat')
% loads all the data for this experiment

T_0 = abs(I3./I02);
% whole data set of T recorded  
T = T_0(349:410);
Lambda=LambdaStore(349:410);
% confining T and lambda to regiojn of interest
% also making sure they are the same size

figure
plot(Lambda,T,'*')
xlabel('Wavelength/nm')
ylabel('Transmission Coefficient')
title({'Transmission Coefficient vs Optical Wavelength taken on',
    'Gallium Phosphide Sample'})
hold off

% btw i dont think they like titles in matlab graphs i think they prefder a
% caption but just check in the formal report guide

n=3.1870;
% value of n quoted from [source] 
R=((n-1)/(n+1))^2;
% approximated value of R calculated 
% Are you not calculating R manually, or just using the estimate?

x = 0.408*10^-3;
% thickness of sample
xerr = 0.001*10^-3;
% error on thickness
h = 6.63*10^-34;
% planks constant
c = 2.99*10^8;
% speed of light
energy_ev = h*c./(Lambda.*10^-9*1.6*10^-19);
% calculates the energy in eV for the respective wavelengths

 
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
rootalpha = -x^(-1).*log((((1 - R).^4 + 4.*(T.^2).*(R.^2)).^0.5 - (1 - R).^2)./(2.*T.*(R.^2)));

figure
plot(energy_ev,rootalpha,'*')
pbs = predint(Alda_low_energy_alpha,energy_ev,0.68);
hold on 
plot(Alda_low_energy_alpha)
hold on
plot(energy_ev,pbs,'m--')
xlabel('Energy / eV')
ylabel('Square root of alpha')
legend('off')
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



%%
rootalphae=rootalpha-rootalphaa;
% % for values with equal energy subtract curves from fittedmodels from each
% % other

hold on
plot(E,rootalphaa);
plot(E,rootalpha);
% % x1 intercept= Eg-Ep CURVE FITTING
hold off
% %% 5)
% plot(E,rootalphae);
% % x2 intercept= Eg+Ep CURVE FITTING
% %% 6)
% Eg=(x1+x2)./2;
% Ep=(x2-x1)./2;
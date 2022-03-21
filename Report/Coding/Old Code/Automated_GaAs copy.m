%% Optical Abosrption MATLAB code - Analysis

import OpticalAnalysisFunctions.nearestRefraction

%% loads all the data from the experiment

% load('Data/Experimental/GaAs - Harry/Combined_As.mat')
load('Data/Experimental/GaAs - Harry/V0_As.mat') 
load('Data/Experimental/GaAs - Harry/V1_As.mat')
load('Data/Experimental/GaAs - Harry/wavelength_As.mat')

%% Refractive Index 

% loads the refractive index info

RefractiveIndex_GaAs = csvimport('Data/Refractive/Papatryfonoset-2021-0.260-1.88-Ga_As.csv'); 
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
    n_As(i) = OpticalAnalysisFunctions.nearestRefraction(L_As, wavelength(i), N_As);
    k_As(i) = OpticalAnalysisFunctions.nearestRefraction(L_As, wavelength(i), K_As);
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

V0 = smooth(V0); % smooths V0 data

x = 0.417.*10.^(-3);                               % thickness of sample
x_err = 0.001;                                     % error on thickness of sample
h = 6.62607004*10^(-34);                           % plancks constant
c = 299792458;                                     % speed of light
Joules_energy = (h*c)./(wavelength.*10.^(-9));     % calculates energy in joules
eV_energy = Joules_energy./(1.602176634*10^(-19)); % converts joules to eV
eV_energy = flip(eV_energy);

% plots V1 and V0 against wavelength
figure('Name', 'Voltage with and without sample')
plot(wavelength,V0)
legend('Voltage without sample', 'Voltage with sample')
xlabel('Wavelength (nm)','Interpreter','latex')
ylabel('Volatage','Interpreter','latex')
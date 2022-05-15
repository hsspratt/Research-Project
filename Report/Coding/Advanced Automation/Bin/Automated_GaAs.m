%% Optical Abosrption MATLAB code - Analysis

% import OpticalAnalysisFunctions
import OpticalAnalysisFunctions.nearestValue
import OpticalAnalysisFunctions.DetectStraightLine
import OpticalAnalysisFunctions.LongestConsecutive
import OpticalAnalysisFunctions.DetectLongestStarightLine
import OpticalAnalysisFunctions.TransparentStraightLine
import OpticalAnalysisFunctions.CutExcessData
%import OpticalAnalysisFunctions.CutExcessDataAlpha


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
    % Approximated value of R calculated
    % Not calculating R manually with this method

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

V0 = smooth(V0); % smooths V0 data
V1 = smooth(V1); % smooths V1 data

V0 = V0*(-1); % volatages negative when they should be +ve
V1 = V1*(-1); % volatages negative when they should be +ve

T = V1./V0;
T_min = min(T).*1.005;
T = T + abs(T_min);

[Excess_Wavelength, Excess_T, wavelength, T, Index] = OpticalAnalysisFunctions.CutExcessData(wavelength, T);
% 
R_As = R_As(Index:1:end); % redifines for new range
V0   = V0(Index:1:end);   % redifines for new range
V1   = V1(Index:1:end);   % redifines for new range


% hold on
% plot(Excess_Wavelength,Excess_T,'*')

[wavelength_trans, T_trans] = OpticalAnalysisFunctions.TransparentStraightLine(wavelength, T);
mean_T = zeros(max(size(wavelength_trans)),1)+mean(T_trans);
figure('Name', 'Title')
plot(wavelength_trans,mean_T,'--')
hold on
plot(wavelength_trans,T_trans)
hold on
plot(wavelength,T)


%R = (1-mean(T))/(1+mean(T)); % Value of R found via transmission graph
% R = (1-meanInsideInterval)/(1+meanInsideInterval) % Value of R found via transmission graph
% meanInsideInterval = mean(T(2070:2255)) 

x = 0.417.*10.^(-3);                               % thickness of sample
x_err = 0.001;                                     % error on thickness of sample
h = 6.62607004*10^(-34);                           % plancks constant
c = 299792458;                                     % speed of light
Joules_energy = (h*c)./(wavelength.*10.^(-9));     % calculates energy in joules
eV_energy = Joules_energy./(1.602176634*10^(-19)); % converts joules to eV
eV_energy = flip(eV_energy);

% calculate each alpha for each wavelength
R=R_As;
alpha_GaAs = -(x.^(-1)).*log((((1 - R).^4 + 4.*(T.^2).*(R.^2)).^0.5 - (1 - R).^2)./(2.*T.*(R.^2)));
alpha_GaAs = flip(alpha_GaAs);

[Excess_eV_energy, Excess_alpha_GaAs, eV_energy, alpha_GaAs, Index] = OpticalAnalysisFunctions.CutExcessData(eV_energy, alpha_GaAs);
wavelength = wavelength(Index:1:end)
V0   = V0(Index:1:end);   % redifines for new range
V1   = V1(Index:1:end);   % redifines for new range
T = T(Index:1:end);

offset_As = min(alpha_GaAs);               % finds alpha offset
offsetalpha = alpha_GaAs + abs(offset_As); % correct alpha offset so no negative values

[xdata, ydata] = OpticalAnalysisFunctions.DetectLongestStarightLine(eV_energy, offsetalpha, 0.07);
figure('Name', 'StraightLine alpha vs eV')
plot(xdata,ydata,'*')
hold on
plot(eV_energy, offsetalpha)

%[Excess_Wavelength_alpha, Excess_alpha, wavelength_alpha, offsetalpha, Index] = OpticalAnalysisFunctions.CutExcessDataAlpha(wavelength, offsetalpha);
% 
%StraightLine = OpticalAnalysisFunctions.DetectStraightLine(flip(eV_energy), flip(offsetalpha))
%Joules_energy = (h*c)./(wavelength_alpha.*10.^(-9));     % calculates energy in joules
%eV_energy = Joules_energy./(1.602176634*10^(-19)); % converts joules to eV
xdata = wavelength;
ydata = offsetalpha;
err = 0.05;

gradient(1) = ydata(1)./(xdata(1));
for i=1:max(size(xdata))-1
    gradient(i) = (ydata(i+1) - ydata(i))./(xdata(i+1) - xdata(i));
end

gradient = smoothdata(gradient);
min_gradient = (ydata(end) - ydata(1))/(xdata(end) - xdata(1));

for i=1:max(size(gradient))
    if gradient(i)<min_gradient
        gradient(i)=0;
    end
end

StraightLine = [];
for i = 1:size(gradient,2)-1
    section = gradient(i:i+1);
    p12 = section(1);
    p23 = section(2);
    %per = p23*err;
    rolling_average = mean(gradient(1:i+1));
    per = rolling_average*err;

    if (p23-per <= p12) && (p12 <= p23+per) && (gradient(i)~=0)
       StraightLine(i) = section(1);
    end
end

if isempty(StraightLine) == 1
    disp('There is no clear straightline')
end


% gradient(1) = offsetalpha(1)./(wavelength(1));
% for i=1:max(size(wavelength))-1
%     gradient(i) = (offsetalpha(i+1) - offsetalpha(i))./(wavelength(i+1) - wavelength(i));
% end
% 
% gradient = smoothdata(gradient);
% StraightLine = [];
% for i = 1:size(gradient,2)-1
%     section = gradient(i:i+1);
%     p12 = section(1);
%     p23 = section(2);
%     per = p23*0.05;
% 
%     if (p23-per <= p12) && (p12 <= p23+per)
%        StraightLine(i) = section(1);
%     end
% end

figure('Name', 'Title')
%plot(eV_energy(1:end),offsetalpha(1:end))
plot(wavelength, offsetalpha)

square_alpha = offsetalpha.^2;             % finds alpha squared
square_alpha = smooth(square_alpha);       % smooths values
Log_alpha_GaAs = log10(offsetalpha);       %log alpha

% plots V1 and V0 against wavelength
figure('Name', 'Voltage with and without sample')
plot(wavelength,V0)
hold on
plot(wavelength,V1,'r')
hold off
legend('Voltage without sample', 'Voltage with sample')
xlabel('Wavelength (nm)','Interpreter','latex')
ylabel('Volatage','Interpreter','latex')

% plots T against wavelength
figure('Name', 'Transmission against wavelength')
plot(wavelength,T,'b')
xlabel('Wavelength (nm)','Interpreter','latex')
ylabel('Transmission Coefficient','Interpreter','latex')

% plots log alpha vs energy
figure('Name', 'Log Alpha against energy (eV)')
plot(eV_energy,real(Log_alpha_GaAs))
xlabel('Optical Energy (eV)','Interpreter','latex')
ylabel('Logaritm of Absorption Coefficient','Interpreter','latex')

% plots alpha squared vs energy, the graph of interest
figure('Name', 'Alpha Squared against energy (eV)')
plot(eV_energy,real(square_alpha))
xlabel('Optical Energy (eV)','Interpreter','latex')
ylabel('Absorption Coefficient Squared','Interpreter','latex')

% % plots theroetical vs experimental values
% figure('Name', 'Alpha Squared against energy (eV)')
% plot(eV_energy(1:40:end),offsetalpha(1:40:end),'*')
% hold on
% plot(eV_energy(1:30:1070),offsetalpha(1:30:1070),'r--')
% %set(gca, 'YScale', 'log')
% legend('experimental values', 'theoretical values')

%% end 

%hold on
%pbs = predint(Alda_low_energy_alpha,energy_ev,0.68);

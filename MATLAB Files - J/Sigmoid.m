%loads in the data
load('Three_Alphas.mat')
load('Alphas_vs_EnergyEv.mat')

% Experimental constants
h = 6.63*10^(-34);
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

deltE = 0.009;

asig = (amax + (amin - amax)./(1+exp((E-E0)./deltE)));
asig = asig';

errpl = err_alpha((l-nmaxrange):l,:);

hold on
p1 = plot(E,asig,'LineWidth',6);
errorbar(E(1:10:end),ai(1:10:end),errpl(1:10:end),'LineWidth',1.5);
p1.Color(4) = 0.5;
xlabel('Energy (eV)')
ylabel('$\alpha$ ','Interpreter','latex')
set(gca,'FontSize',16)


%%
nboltzdir = 0.3;
nboltzindir = 4.3;
Eboltz = E0 - (nboltzdir*deltE);

%%
amaxerr = err_alpha(naimax);
aminerr = err_alpha(naimin);

amaxmax = amaxerr + amax;
aminmax = amin - aminerr;

ahalfmax = (amaxmax - aminmin)/2;
aposmax = abs(ahalfmax-ai);
aposminmax = min(aposmax);
naposmax = find(aposminmax==aposmax);
E0max = E(naposmax);

amaxmin = amax - amaxerr;
aminmin = amin + aminerr;

ahalfmin = (amaxmin - aminmin)/2;
aposmin = abs(ahalfmin-ai);
aposminmin = min(aposmin);
naposmin = find(aposminmin==aposmin);
E0min = E(naposmin);

E0_err = E0max - E0min;

toterr = sqrt((E0_err*(1+0.3*deltE))^2 + ((0.0102-0.0043)*(E0 + 0.3))^2)

%%

wavemax = h*c./(E(naimax).*10^(-9)*1.6*10^-19);
wavemin = h*c./(E(naimin).*10^(-9)*1.6*10^-19);

wavemaxmax = wavemax - 3;
wavemaxmin = wavemax + 3;

waveminmax = wavemin - 3;
waveminmin = wavemin + 3;

Emaxmax = h*c./(wavemaxmax.*10^(-9)*1.6*10^-19);
Emaxmin = h*c./(wavemaxmin.*10^(-9)*1.6*10^-19);
Eminmax = h*c./(waveminmax.*10^(-9)*1.6*10^-19);
Eminmin = h*c./(waveminmin.*10^(-9)*1.6*10^-19);


gradmaxmaxco = [Emaxmin,amaxmax];
gradmaxminco = [Eminmax,aminmax];


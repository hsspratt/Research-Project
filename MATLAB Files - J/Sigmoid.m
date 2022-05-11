%gaasdatfile = importdata('Three_Alphas.mat');

%gaasdat = gaasdatfile(:,1);

%transfile = importdata('Nosample.txt');

h = 6.63*10^-34;
c = 3*10^8;
%%

a = smooth(alpha_1);

amax = max(a);
amin = min(a);

%From here down 
nmax = find(amax==a);
%nmin = find(amin==a);

perc = 0.2;
nmaxrange = round(nmax);
%nminrange = round(nmin - perc*nmin);

l = length(a);
ai = a((l-nmaxrange):l,:);

ahalf = (amax - amin)/2;
apos = abs(ahalf-ai);
aposmin = min(apos);
napos = find(aposmin==apos);

la = length(ai);
na = la;
Ewave = (350+na:la);
E = h*c./(Ewave.*10^(-9)*1.6*10^-19);
E0 = E(100);
% 543

deltE = 0.01;

asig = (amax + (amin - amax)./(1+exp((E-E0)./deltE)));
asig = asig';

hold on
plot(E,asig)
plot(E,ai)

%%
nboltz = 0.3;
Eboltz = E0 - (nboltz*deltE);

%%
lam = (350:1089);


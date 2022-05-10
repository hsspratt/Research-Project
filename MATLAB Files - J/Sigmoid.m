gaasdatfile = importdata('GaAs_Transmission.txt');

gaasdat = gaasdatfile(:,1);

transfile = importdata('Nosample.txt');

%%
adata = importdata('');
a = adata(:,1);

amax = max(a);
amin = min(a);

ahalf = (amax - amin)/2;
apos = find(a == ahalf);

lam = [350,1089,1];
E0 = E(apos);

asig = amax + (amin - amax)/(1+exp((E-E0)/deltE));


%%
nboltz = 0.3;
Eboltz = E0 - nboltz*deltE;

%%
lam = (350:1089);


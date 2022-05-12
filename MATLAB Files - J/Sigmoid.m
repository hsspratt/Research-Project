%gaasdatfile = importdata('Three_Alphas.mat');

%gaasdat = gaasdatfile(:,1);

%transfile = importdata('Nosample.txt');
load('/Users/jakehurley/Documents/GitHub/Research-Project/MATLAB Files - J/Three_Alphas.mat')
load('/Users/jakehurley/Documents/GitHub/Research-Project/MATLAB Files - J/Alphas_vs_EnergyEv.mat')

h = 6.63*10^-34;
c = 3*10^8;
%%

a = smooth(alpha_1);

amax = max(a);
amin = min(a);

%From here down 
nmax = find(amax==a);


perc = 0.5;
nmaxrange = round(nmax-150);
%nminrange = round(nmin - perc*nmin);

l = length(a);
ai = a((l-nmaxrange):l,:);

ahalf = (amax - amin)/2;
apos = abs(ahalf-ai);
aposmin = min(apos);
napos = find(aposmin==apos);

la = length(ai);
na = 1089 - la+1;
Ewave = (na:1089);
E = h*c./(Ewave.*10^(-9)*1.6*10^-19);
E0 = E(napos);
%543

%aimin = min(ai);

%naimin = find(aimin == ai);
%naimax = find(amax == ai);

%%gradient(1) = ai(1)./(E(1));
        %for i=1:max(size(E))-1
            %gradient(i) = (ai(i+1) - ai(i))./(E(i+1) - E(i));
        %end

 %gradmax = max(gradient);
 %ngradmax = find(gradmax == gradient);

%Eminslope = E(ngradmax);
%Emaxslope = E(naimax);



%deltE = (Emaxslope - Eminslope)/2;
deltE = 0.009;

asig = (amax + (amin - amax)./(1+exp((E-E0)./deltE)));
asig = asig';

hold on
p1 = plot(E,asig,'LineWidth',6);
plot(E,ai,'LineWidth',2);
p1.Color(4) = 0.5;
xlabel('Energy (eV)')
ylabel('$\alpha$ ','Interpreter','latex')
set(gca,'FontSize',16)


%%
nboltzdir = 0.3;
nboltzindir = 4.3;
Eboltz = E0 - (nboltzdir*deltE);

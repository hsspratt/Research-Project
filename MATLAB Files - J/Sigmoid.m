%gaasdatfile = importdata('Three_Alphas.mat');

%gaasdat = gaasdatfile(:,1);

%transfile = importdata('Nosample.txt');
load('Three_Alphas.mat')
load('Alphas_vs_EnergyEv.mat')

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
aminmin = amin - aminerr;

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

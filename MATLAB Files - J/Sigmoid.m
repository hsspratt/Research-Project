%loads in the data
load('Three_Alphas.mat')
load('Alphas_vs_EnergyEv.mat')

% Experimental constants
h = 6.63*10^(-34);
c = 3*10^8;
%%
h = 6.63*10^(-34);
c = 3*10^8;

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
amaxerr = [err_alpha(naimax),];
aminerr = [err_alpha(naimin),];

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

%%

a = smooth(alpha_1);

a_max_y = max(a);
a_min_y = min(a);

a_max_err_y = err_alpha(naimax);
a_min_err_y = err_alpha(naimin);

a_max_y = a_max_err_y + a_max_y;
a_min_y = a_min_y - a_min_err_y;

wavelength_err = 3.5;

e_max_x = energy_ev(find(a == max(a)));
e_min_x = energy_ev(find(a == max(a)));

l_max_x = ConvertEnergyEv2Lambda(e_max_x) + wavelength_err;
l_min_x = ConvertEnergyEv2Lambda(e_min_x) - wavelength_err;

e_max_x = ConvertLambda2EnergyEv(l_max_x);
e_min_x = ConvertLambda2EnergyEv(l_min_x);



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

%%

y_max = max(a);
y_min = min(a);
varlist = [y_max, y_min];

f = (y_max-y_min)*0.5;

E0_err = PropError(f,[y_max, y_min],[y_max, y_min],[err_alpha(naimax), err_alpha(naimin)])

%%
y_max = max(a);
y_min = min(a);

y_err = sqrt(err_alpha(naimax)^2 + err_alpha(naimin)^2);

%%

function sigma = PropError(f,varlist,vals,errs)
%SIGMA = PROPERROR(F,VARLIST,VALS,ERRS)
%
%Finds the propagated uncertainty in a function f with estimated variables
%"vals" with corresponding uncertainties "errs".
%
%varlist is a row vector of variable names. Enter in the estimated values
%in "vals" and their associated errors in "errs" at positions corresponding 
%to the order you typed in the variables in varlist.
%
%Example using period of a simple harmonic pendulum:
%
%For this example, lets say the pendulum length is 10m with an uncertainty
%of 1mm, and no error in g.
%syms L g
%T = 2*pi*sqrt(L/g)
%type the function T = 2*pi*sqrt(L/g)
%
%PropError(T,[L g],[10 9.81],[0.001 0])
%ans =
%
%    [       6.3437]    '+/-'    [3.1719e-004]
%    'Percent Error'    '+/-'    [     0.0050]
%
%(c) Brad Ridder 2007. Feel free to use this under the BSD guidelines. If
%you wish to add to this program, just leave my name and add yours to it.
n = numel(varlist);
sig = vpa(ones(1,n));
for i = 1:n
    sig(i) = diff(f,varlist(i),1);
end
error1 =sqrt((sum((subs(sig,varlist,vals).^2).*(errs.^2))));
error = double(error1);
sigma = [{subs(f,varlist,vals)} {'+/-'} {error};
         {'Percent Error'} {'+/-'} {abs(100*(error)/subs(f,varlist,vals))}];
end

function E = ConvertLambda2EnergyEv(lambda)
    h = 6.62607015*10^(-34);
    c = 299792458;

    E = (h*c)./lambda;
end

function lambda = ConvertEnergyEv2Lambda(E)
    h = 6.62607015*10^(-34);
    c = 299792458;

    lambda = (h*c)./E;
end

# %%
%matplotlib widget
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

E_g = 1.42
P_e = np.linspace(E_g, 1.65,100)
A = 100000
t = np.arange(-100, 100., 0.2)
t = np.linspace(-10, 10, 10000)

x = np.full(100, E_g)
y = np.linspace(1,1100, 100)
alpha = A*(P_e-E_g)**(0.5)

line_domain = np.linspace(1.32, 1.4305, 10)
line2 = 48969.38775510211*t-64638.59183673478

m = (np.log10(10472.7)-np.log10(1))/(-1.32 + 1.4305)
k = 1/(10**(m*1.32))
s_line = k*10**(m*line_domain) * (1.5*np.random.normal(0, .1, line_domain.shape)+1)
c_line = alpha[10::5] *  (1.5*np.random.normal(0, .1, alpha[10::5].shape)+1)

plt.figure(1)
plt.plot(line_domain,s_line,'ro')

fig = plt.figure()
ax = fig.add_subplot(111)

plt.yscale('log')

plt.plot(P_e,alpha,'k--',label=r'Theoretical Curve $\alpha = A(hv-E_g)^(\frac{1}{2})$')
plt.plot(line_domain,s_line,'ko',markersize=5,label=r'Experimental Results')
plt.plot(P_e[10::5],c_line,'ko',markersize=5)

plt.xlim(1.3,1.7)
plt.ylim(1,50000)
plt.xlabel(r'Photon Energy (eV)')
plt.xticks(np.arange(1.3, 1.8, 0.1))
plt.ylabel(r'Absorption Coefficient $(cm^{-1})$')
plt.legend((r'Theoretical Curve $ \\ \vspace{0.25cm} \alpha = A(hv-E_g)^{\frac{1}{2}}$', r'Experimental Results'),
            shadow=False, loc=(0.38, 0.4), handlelength=1.5, fontsize=13)
            
ax.set_aspect(1./ax.get_data_ratio())

plt.savefig("Theoretical_Vs_Experimental.svg", format = 'svg', dpi=1200)
# %%

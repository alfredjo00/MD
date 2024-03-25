import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.linear_model import LinearRegression
t0 = 3.335641e-19 * 10**12


sns.set('talk')
N_eq = 10000
N_a0 = 1500

state = 'liq'

data = np.loadtxt(f"data/task4_{state}.txt")

def plot_main_liq():
	fig,ax = plt.subplots(3,figsize=(10,10))

	#* 160.21766208 / 0.0001
	
	print("P in bar", np.mean(data[N_eq:,3]) * 160.21766208 / 0.0001, np.std(data[N_eq:,3]) * 160.21766208 / 0.0001)
	print("a0 mean ", np.mean(data[N_eq:,5]))
	print("Potential e", np.mean(data[N_eq:,2]), np.std(data[N_eq:,2]))
	t = data[:,0] * t0
	ax[0].plot(t, data[N_eq:,1], "-b", markerfacecolor='none', 
	    label='Potential energy')
	ax[1].plot(t, data[N_eq:,1] + data[:,2], "-r", markerfacecolor='none',
	    label='Total energy')
	ax[2].plot(t, data[N_eq:,2], "-k", markerfacecolor='none',
	    label='Kinetic energy')

	for a in ax:
		a.set_xlabel("Time (ps)")
		a.set_ylabel("Energy (eV)")
		a.legend()

	fig.tight_layout()
	fig.savefig(f"graphs/task4_energy_{state}.jpg", dpi=90)
	
	
def plot_temps_pressures():
	fig,ax = plt.subplots(2,figsize=(10,10))
	T = data[:,4] - 273.15
	P = data[:,3] * 160.21766208 / 0.0001
	M  = len(P[N_eq:]) / 100
	
	print("P mean", np.mean(P[N_eq:]), np.std(P[N_eq:])/np.sqrt(M))
	print("T mean", np.mean(T[N_eq:]), np.std(T[N_eq:])/np.sqrt(M))
	
	print("a0 mean ", np.mean(data[N_a0:N_eq,5]), np.std(data[N_a0:N_eq,5]))
	t = data[:,0] * t0
	
	ax[0].plot(t, P, "-b", markerfacecolor='none', 
	    label='Pressure')
	ax[1].plot(t, T, "-r", markerfacecolor='none',
	    label='Temperature')

	ax[0].set_ylabel('Pressure (bar)')
	ax[1].set_ylabel(r'Temperature (C$^\circ$)')
	
	ax[0].set_ylim(-15000,15000)
	ax[1].set_ylim(-50,1500)
	for a in ax:
		a.set_xlabel("Time (ps)")
		a.legend()

	fig.tight_layout()
	fig.savefig(f"graphs/task4_{state}_PT.jpg", dpi=90)
	
	
def plot_pos_liq():
	fig,ax = plt.subplots(3,figsize=(10,10))

	data = np.loadtxt(f"data/task4_{state}_pos.txt")
	t = data[:,0] * t0
	p1 = [data[:,1], data[:,2], data[:,3]]
	p2 = [data[:,4], data[:,5], data[:,6]]
	p3 = [data[:,7], data[:,8], data[:,9]]
	ps = [p1,p2,p3]
	tils = [r'$p_1$', '$p_{50}$', '$p_{128}$']
	
	for a,p in zip(ax,ps):
		a.plot(t , p[0], "-b", label='$x$')
		a.plot(t , p[1], "-r", label='$y$')
		a.plot(t , p[2], "-k", label='$z$')
	    
	for i,a in enumerate(ax):
		a.set_title(tils[i])
		a.set_xlabel("Time (ps)")
		a.set_ylabel(u"Position (\u00C5)")
		a.legend()

	fig.tight_layout()
	fig.savefig(f"graphs/task4_{state}_pos.jpg", dpi=90)

def plot_msd_liq():	
	msd_ls = data[:N_eq,-2]
	t = data[:N_eq,0]
	b = 0
	reg = LinearRegression().fit(t[b:].reshape(-1,1)* t0, msd_ls[b:].reshape(-1,1))
	
	pred = reg.predict( t[b:].reshape(-1,1)* t0)
	
	print("Coefficients: \n", reg.coef_)
	print("Prediction:", pred)
	print()
	print("Diff const ",  (pred[-1][0] ) / (6 * data[-1,0] * t0 ))
	fig,ax = plt.subplots(figsize=(8,6))
	ax.plot(t* t0,msd_ls, '-b', label=r'$\Delta_{MSD}(t)$')
	ax.plot(t[b:] * t0, pred, '--k', label='Linear fit')
	ax.set_xlabel("Time (ps)")
	ax.legend()
	ax.set_ylabel(u"Mean squared displacement (\u00C5$^2$)")
	
	fig.tight_layout()
	fig.savefig(f"graphs/{state}_msd.jpg", dpi=90)
	
def plot_vcf_liq():	
	sns.set('talk')	
	b = 100
	t = data[:b,0] * t0
	vcf = data[:b,-1] / data[0,-1]
	fig,ax = plt.subplots(figsize=(8,6))
	ax.plot(t,vcf, '-k')
	
	ax.set_xlabel("Time (ps)")
	
	ax.set_ylabel(u"Velocity correlation function (a.u.)")
	
	fig.tight_layout()
	fig.savefig(f"graphs/{state}_vcf.jpg", dpi=90)
	
def plot_spec():	
	sns.set('talk')
	fft_data = np.loadtxt("data/fft_data.txt")
	
	b = 100
	t = fft_data[:b,0]
	ifft = fft_data[:b,1] / fft_data[0,1]
	
	f = fft_data[:,2]
	spec = np.sqrt(fft_data[:,3])
	
	fig,ax = plt.subplots(2,figsize=(12,7))
	ax[0].set_title(r"Velocity correlation function - $\mathcal{F}^{-1} \left[  \hat{\Phi}(f)^2  \right]$")
	ax[1].set_title(r"FT of velocity correlation function - $\hat{\Phi}(f)$")
		
	ax[0].plot(t, ifft, '-k')
	ax[1].plot(f, spec, '-k')
	
	ax[0].set_yscale("log")
	ax[1].set_yscale("log")
	
	ax[0].set_xlabel("Time (ps)")
	ax[0].set_ylabel(r"$\Phi(t)$ (a.u)")
	
	ax[1].set_xlabel(r"$f$ (THz)")
	ax[1].set_ylabel(r"$ \mathcal{F}[\Phi(t)] $")
	
	fig.tight_layout()
	fig.savefig("graphs/liq_fft.jpg", dpi=100)

	
	
if __name__ == "__main__":
	plot_temps_pressures()
	plot_pos_liq()
	plot_msd_liq()
	plot_vcf_liq()
	plot_spec()
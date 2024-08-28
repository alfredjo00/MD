import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
t0 = 3.335641e-19 * 10**12

def plot_energy_task3():
	sns.set_theme('talk')

	data = np.loadtxt("data/task4_solid.txt")
	fig,ax = plt.subplots(3,figsize=(10,10))

	#* 160.21766208 / 0.0001
	
	print("P in bar", np.mean(data[1000:,3]) * 160.21766208 / 0.0001)
	print("a0 mean ", np.mean(data[1000:,5]))
	t = data[:,0] * t0
	ax[0].plot(t, data[:,1], "-b", markerfacecolor='none', 
	    label='Potential energy')
	ax[1].plot(t, data[:,1] + data[:,3], "-r", markerfacecolor='none',
	    label='Total energy')
	ax[2].plot(t, data[:,3], "-k", markerfacecolor='none',
	    label='Kinetic energy')

	for a in ax:
		a.set_xlabel("Time (ps)")
		a.set_ylabel("Energy (eV)")
		a.legend()

	fig.tight_layout()
	fig.savefig("graphs/task3.jpg", dpi=90)
	
	
def plot_positions_task3():
	sns.set_theme('talk')

	data = np.loadtxt("data/task4_solid_pos.txt")
	fig,ax = plt.subplots(3,figsize=(10,10))


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
	fig.savefig("graphs/task3_positions.jpg", dpi=90)

def plot_msd():	
	sns.set_theme('talk')

	data = np.loadtxt("data/task4_solid_pos.txt")
	
	t = data[:,0] * t0
	msd = data[:,-1]
	fig,ax = plt.subplots(figsize=(8,6))
	ax.plot(t,msd, '-k')
	ax.set_xlabel("Time")
	
	ax.set_ylabel(u"Mean squared displacement (\u00C5)")
	
	fig.tight_layout()
	fig.savefig("graphs/task3_msd.jpg", dpi=90)
	
	

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def plot_energy_kinetic():
    sns.set_theme('talk')

    data = np.loadtxt("data/task2_0.txt")
    fig,ax = plt.subplots(3,figsize=(10,10))

    b = 5
    ax[0].plot(data[::b,0], data[::b,1], "-b", markerfacecolor='none', 
        label='Potential energy')
    ax[1].plot(data[::b,0], data[::b,1] + data[::b,3], "-r", markerfacecolor='none',
        label='Total energy')
    ax[2].plot(data[::b,0], data[::b,3], "-k", markerfacecolor='none',
        label='Kinetic energy')

    for a in ax:
        a.set_xlabel("Time (ps)")
        a.set_ylabel("Energy (eV)")
        a.legend()

    fig.tight_layout()
    fig.savefig("graphs/task2_0.jpg", dpi=90)
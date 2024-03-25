import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set('talk')

data = np.loadtxt("data/volume_Al-energy_pot.txt")
fig,ax = plt.subplots(figsize=(6,6))
x,y = data[:,0], data[:,1]
ax.plot(x, y, ":ob", label='Data')
ax.set_xlabel(u"Volume (\u00C5$^3$)")
ax.set_ylabel("Energy (eV / unit cell)")

z = np.polyfit(x, y, deg=2)
print(z)
p = np.poly1d(z)
f_z = [p(t) for t in x]
ax.plot(x, f_z, '--k' ,label='Quadratic fit')
i = np.argmin(f_z)
# Lowest energy point
print(x[i], x[i] ** (1/3))

ax.legend()
fig.tight_layout()
fig.savefig("graphs/task1.jpg", dpi=90)
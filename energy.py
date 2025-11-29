import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("energy_1000_leapfrog.txt")

time = data[:,0]
kinetic_energy = data[:,1]
potential_energy = data[:,2]
whole_energy = data[:,1] + data[:,1]

fig =plt.figure(figsize=(10,10))
plt.plot(time*10000000,kinetic_energy,linewidth=0,marker=".")
plt.xlabel("time")
plt.ylabel("kinetic energy")

plt.title("time-energy")
plt.axis("equal")
plt.grid(True)

#plt.xlim(0,1)
#plt.ylim(0,200000)

plt.savefig("Kenergy_1000_leapfrog.png")   

plt.show()
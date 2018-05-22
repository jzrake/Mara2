import numpy as np
import matplotlib.pyplot as plt
import math, os

x = np.linspace(0,20,100)
y = np.sin(x * 2.0 * math.pi / 10.0)



plt.plot(x,y)
plt.xlabel("Distance in x (m)")
plt.ylabel("Vertical Displacement in y (m)")
plt.title("Periodic wave in two regions")
plt.savefig("Wave.png")
plt.show()

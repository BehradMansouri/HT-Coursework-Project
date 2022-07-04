import matplotlib
from matplotlib import pyplot as plt
import numpy as np

# M = delta_tau / (delta_eta ** 2)
# M = 0.0001 / (0.05 ** 2)
M = 0.04
eta = np.arange(0, 1.05, 0.05, dtype=float)  # discretizing the space array from 0 to 1
tau = np.arange(0, 1.0001, 0.0001, dtype=float)  # discretizing the time array from 0 to 1
theta = np.zeros((10001, 21, 6))
one = []
two = []
three = []
four = []
five = []
six = []
Qone = np.zeros(10001)
Qtwo = np.zeros(10001)
Qthree = np.zeros(10001)
Qfour = np.zeros(10001)
Qfive = np.zeros(10001)
Qsix = np.zeros(10001)


# Initial conditions
for q in range(0, 5+1):
    for b in range(0, 20+1):
        theta[0, b, q] = 1
for v in range(0, 5+1):
    for w in range(1, 10000+1):
        theta[w, 0, v] = 0
        theta[w, 20, v] = 0

P = [-0.9, -0.5, 0, 1, 3, 5]

for j in range(0, 9999+1):    # for different times
    for k in range(1, 19+1):      # for different points
        for i in range(0, 5+1):       # for different P values
            theta[j + 1, k, i] = (theta[j, k, i] * (1 - 2 * M)) + (M * (theta[j, k + 1, i] + theta[j, k - 1, i])) + (P[i] * M * theta[j, k, i] * (theta[j, k + 1, i] + theta[j, k - 1, i] - (2 * theta[j, k, i]))) + (0.25 * P[i] * M * ((theta[j, k + 1, i] ** 2) + (theta[j, k - 1, i] ** 2) - (2 * theta[j, k + 1, i] * theta[j, k - 1, i])))
            # formulas for points k=eta=0 and k=20 or eta=1 are not needed since their temperatures are always constant

for v in range(0, 10000+1):
    one.append(theta[v, 10, 0])
    two.append(theta[v, 10, 1])
    three.append(theta[v, 10, 2])
    four.append(theta[v, 10, 3])
    five.append(theta[v, 10, 4])
    six.append(theta[v, 10, 5])
    for w in range(1, 19 + 1):
        Qone[v] += (theta[v, w, 0]) * 0.05
        Qtwo[v] += (theta[v, w, 1]) * 0.05
        Qthree[v] += (theta[v, w, 2]) * 0.05
        Qfour[v] += (theta[v, w, 3]) * 0.05
        Qfive[v] += (theta[v, w, 4]) * 0.05
        Qsix[v] += (theta[v, w, 5]) * 0.05
    Qone[v] += ((theta[v, 0, 0]) + (theta[v, 20, 0])) * 0.025
    Qtwo[v] += ((theta[v, 0, 1]) + (theta[v, 20, 1])) * 0.025
    Qthree[v] += ((theta[v, 0, 2]) + (theta[v, 20, 2])) * 0.025
    Qfour[v] += ((theta[v, 0, 3]) + (theta[v, 20, 3])) * 0.025
    Qfive[v] += ((theta[v, 0, 4]) + (theta[v, 20, 4])) * 0.025
    Qsix[v] += ((theta[v, 0, 5]) + (theta[v, 20, 5])) * 0.025

fig = plt.figure(figsize=(16, 9))
plt.subplot(221)
plt.plot(tau, one)
plt.plot(tau, two)
plt.plot(tau, three)
plt.plot(tau, four)
plt.plot(tau, five)
plt.plot(tau, six)
plt.title('Fig. 1: Centerline temperature')
matplotlib.pyplot.ylabel('(Tc-Tb)/(Ti-Tb)')
plt.grid(True)

plt.subplot(222)
plt.plot(tau, Qone, label='P=-0.9')
plt.plot(tau, Qtwo, label='P=-0.5')
plt.plot(tau, Qthree, label='P=0')
plt.plot(tau, Qfour, label='P=1')
plt.plot(tau, Qfive, label='P=3')
plt.plot(tau, Qsix, label='P=5')
plt.title('Fig. 2: Dimensionless internal energy')
matplotlib.pyplot.ylabel('Q\'\'/Q0\'\'')
plt.legend()
plt.grid(True)

plt.subplot(223)
plt.plot(tau, one)
plt.plot(tau, two)
plt.plot(tau, three)
plt.plot(tau, four)
plt.plot(tau, five)
plt.plot(tau, six)
plt.yscale('log')
plt.title('Fig. 3: Centerline temperature - logarithmic')
matplotlib.pyplot.xlabel('Dimensionless time')
matplotlib.pyplot.ylabel('(Tc-Tb)/(Ti-Tb)')
plt.grid(True)

plt.subplot(224)
plt.plot(tau, Qone)
plt.plot(tau, Qtwo)
plt.plot(tau, Qthree)
plt.plot(tau, Qfour)
plt.plot(tau, Qfive)
plt.plot(tau, Qsix)
plt.yscale('log')
plt.title('Fig. 4: Dimensionless internal energy - logarithmic')
matplotlib.pyplot.xlabel('Dimensionless time')
matplotlib.pyplot.ylabel('Q\'\'/Q0\'\'')
plt.grid(True)
plt.show()

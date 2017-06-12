import numpy as np
import matplotlib.pyplot as plt

from collections import deque

K = 1
n = 20
rho = 0.5
lmb = 1
dt = 0.1

print("== Question 1 ==")

def simulate(rho=rho, T=1000, full=False):
    mu = lmb / rho

    N = int(T / dt)

    X = np.zeros(N)
    tEvent = np.random.exponential(scale=1/lmb)

    for k in range(1, N):
        X[k] = X[k - 1]
        t = k * dt

        if t >= tEvent:
            if X[k] == 0:
                X[k] += 1
            else:
                direction = np.random.uniform(0, 1)

                if direction < lmb / (lmb + mu):
                    X[k] += 1
                else:
                    X[k] -= 1

            if X[k] == 0:
                tEvent += np.random.exponential(scale=1/lmb)
            else:
                tEvent += np.random.exponential(scale=1/(lmb + mu))

    # t = np.arange(T)
    # plt.step(t, X)
    # plt.show()

    if full:
        return X
    else:
        return X[-1]


plt.step(np.arange(int(100 / dt)), simulate(rho=0.5, T=100, full=True))
plt.title(r"$\rho = 0.5$")
plt.show()

plt.step(np.arange(int(100 / dt)), simulate(rho=1, T=100, full=True))
plt.title(r"$\rho = 1$")
plt.show()

plt.step(np.arange(int(100 / dt)), simulate(rho=2, T=100, full=True))
plt.title(r"$\rho = 2$")
plt.show()

print("== Question 2 ==")

N_sim = 1250
trajectoires = np.zeros(N_sim)

for k in range(N_sim):
    trajectoires[k] = simulate()

print("Moyenne et variance théoriques : ", rho / (1 - rho), rho / (1 - rho) ** 2)
print("Moyenne et variance empiriques : ", np.mean(trajectoires), np.std(trajectoires) ** 2)

pi_1 = np.zeros(n)
for k in range(n):
    pi_1[k] = np.mean(trajectoires == k)

moyenne_pi_1 = np.dot(pi_1, np.arange(n))
variance_pi_1 = np.dot(pi_1, np.arange(n) ** 2) - moyenne_pi_1 ** 2
print("Moyenne et variance par théorème ergodique (premier point) : ", moyenne_pi_1, variance_pi_1)
print("PI : ", pi_1)

X = simulate(T=100000, full=True)
pi_2 = np.zeros(n)
for k in range(n):
    pi_2[k] = np.mean(X == k)

moyenne_pi_2 = np.dot(pi_2, np.arange(n))
variance_pi_2 = np.dot(pi_2, np.arange(n) ** 2) - moyenne_pi_2 ** 2
print("Moyenne et variance par théorème ergodique (second point) : ", moyenne_pi_2, variance_pi_2)
print("PI : ", pi_2)

pi_th = rho ** np.arange(n) * (1 - rho)
print("PI théorique : ", pi_th)

print("== Question 2 / Convergence ==")

dts = [2, 1, 0.75, 0.5, 0.25, 0.1, 0.05, 0.01]
moyennes = []
variances = []
for deltaT in dts:
    dt = deltaT
    X = simulate(T=100000, full=True)
    pi_2 = np.zeros(n)
    for k in range(n):
        pi_2[k] = np.mean(X == k)

    moyenne_pi_2 = np.dot(pi_2, np.arange(n))
    variance_pi_2 = np.dot(pi_2, np.arange(n) ** 2) - moyenne_pi_2 ** 2

    moyennes.append(moyenne_pi_2)
    variances.append(variance_pi_2)

plt.semilogx(dts, moyennes, label="Moyenne ergodique")
plt.axhline(y=rho / (1 - rho), label="Moyenne théorique", color="g")
plt.semilogx(dts, variances, label="Variance ergodique")
plt.axhline(y=rho / (1 - rho) ** 2, label="Variance théorique", color="r")
plt.legend()
plt.show()


print("== Question 3 ==")

plt.hist(trajectoires, normed=True, label="Histogramme de $X_t$")
plt.step(np.arange(n), pi_1, where='post', label=r"$\pi_1$")
plt.step(np.arange(n), pi_2, where='post', label=r"$\pi_2$")
plt.step(np.arange(n), pi_th, where='post', label=r"$\pi_{th}$")
plt.legend()
plt.show()

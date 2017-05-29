import numpy as np
import matplotlib.pyplot as plt

K = 1
n = 50
rho = 0.5
lmb = 1
dt = 1

print("== Question 1 ==")

def simulate(rho=rho, T=5000, full=False):
    mu = lmb / rho

    X = np.zeros(T)
    tEvent = np.random.exponential(scale=1/(lmb + mu))

    for t in range(1, T):
        X[t] = X[t - 1]

        # if t >= tM and X[t] > 0:
        #     X[t] = X[t] - 1
        #     tM += np.random.exponential(scale=1/mu)
        #
        # if t >= tL:
        #     if X[t] == 0:
        #         tM += np.random.exponential(scale=1/mu)
        #
        #     X[t] += 1
        #     tL += np.random.exponential(scale=1/lmb)

        if t >= tEvent:
            tEvent += np.random.exponential(scale=1/(lmb + mu))

            direction = np.random.uniform(0, 1)

            if direction <= lmb / (lmb + mu) or X[t] == 0:
                X[t] += 1
            else:
                X[t] -= 1

    # t = np.arange(T)
    # plt.step(t, X)
    # plt.show()

    if full:
        return X
    else:
        return X[-1]

plt.step(np.arange(100), simulate(rho=0.5, T=100, full=True))
plt.title(r"$\rho = 0.5$")
plt.show()

plt.step(np.arange(100), simulate(rho=1, T=100, full=True))
plt.title(r"$\rho = 1$")
plt.show()

plt.step(np.arange(100), simulate(rho=2, T=100, full=True))
plt.title(r"$\rho = 2$")
plt.show()

print("== Question 2 ==")

N_sim = 1000
trajectoires = np.zeros(N_sim)

for k in range(N_sim):
    trajectoires[k] = simulate()

print("Moyenne et variance théoriques : ", rho / (1 - rho), rho / (1 - rho) ** 2)
print("Moyenne et variance empiriques : ", np.mean(trajectoires), np.std(trajectoires) ** 2)

X = simulate(full=True)

pi = np.zeros(n)
for k in range(n):
    pi[k] = np.mean(X == k)

moyenne_pi = np.dot(pi, np.arange(n))
variance_pi = np.dot(pi, np.arange(n) ** 2) - moyenne_pi ** 2
print("Moyenne et variance par théorème ergodique : ", moyenne_pi, variance_pi)
print("Pi : ", pi)

print("== Question 3 ==")

plt.hist(trajectoires, normed=True, label="Histogramme de $X_t$")
plt.hist(pi, normed=True, label=r"$\pi$")
plt.legend()
plt.show()

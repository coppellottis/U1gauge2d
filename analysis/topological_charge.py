import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullLocator
import emcee
from scipy.stats import norm
from scipy.integrate import quad
from scipy.special import i0

def tau_int(rho,c):
    tau = .5+np.cumsum(rho[1:])
    M = np.arange(1, len(rho))

    idx = np.where(M >= c*tau)[0]
    if len(idx) == 0:
        print(f"Nessuna finestra autoconsistente trovata")
        return 0,0

    M_star = M[idx[0]]
    tau_star = tau[idx[0]]

    return tau_star, M_star

def phi2(beta):
    integrand = lambda phi: phi**2 * np.exp(beta*np.cos(phi))
    integral, _ = quad(integrand, -np.pi, np.pi)

    return integral / (2*np.pi*i0(beta))

def chi(beta):
    return 1/(4*(np.pi**2))*phi2(beta)

def block_obs(data, b_size, function, iterations = 300):
    b_number = len(data) // b_size
    blocks = data[:b_size*b_number].reshape(b_number, b_size)

    rng = np.random.default_rng()

    # Un ricampionamento alla volta: in memoria c'e' sempre un solo campione
    # (len(data) elementi) invece di tutti i `iterations` insieme (iterations x len(data)).
    estimates = np.empty(iterations)
    for i in range(iterations):
        idx = rng.integers(0, b_number, size=b_number)
        sample = blocks[idx].reshape(-1)
        estimates[i] = np.mean(function(sample))   # function puo' restituire un array (es. Q**2): ne prendo la media

    return np.mean(estimates), np.std(estimates, ddof=1)

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 20,

    "axes.labelsize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,

    "axes.linewidth": 1.0,

    "xtick.direction": "in",
    "ytick.direction": "in",

    "xtick.top": True,
    "ytick.right": True,

    "xtick.major.size": 6,
    "ytick.major.size": 6,
    "xtick.minor.size": 3,
    "ytick.minor.size": 3,

    "legend.fontsize": 20,
})

fig, ax = plt.subplots(4,1, figsize=(8, 16))

c = np.array([4,3,3,2])
L = np.array([16, 24, 28, 32])
beta = np.array([2,4,6,8])

arr = []
arr2 = []

for i, l in enumerate(L):
    Q = pd.read_csv(f"data/sigma/Q_W_L{l}_beta{beta[i]}.csv").iloc[:, 0].to_numpy()
    mean = Q.mean()
    rho = emcee.autocorr.function_1d(Q)
    tau, _ = tau_int(rho, c[i])
    err = Q.std(ddof=1) * np.sqrt(2*tau/len(Q))

    mean2, err2 = block_obs(data=Q, b_size=int(5*np.ceil(tau)), function=lambda x: x**2, iterations=100)

    QQ = np.arange(-7,8)
    hist = [np.sum(Q == i)/np.size(Q,axis=0) for i in QQ]
    ax[i].bar(QQ, hist, width=0.8, alpha=.5,label=fr'$\beta={beta[i]}$')

    x = np.linspace(-7, 7, 500)
    sigma = np.sqrt(l**2/(4*np.pi**2)*phi2(beta[i]))
    y = norm.pdf(x, loc=0, scale=sigma)
    ax[i].plot(x,y, '-.', color='red')
    #sigma2 = np.sqrt(l**2/(4*np.pi**2*beta[i]))
    #y2 = norm.pdf(x, loc=0, scale=sigma2)
    #ax[i].plot(x,y2, '-.', color='orange')

    ax[i].set_xlabel(fr"$Q$")
    ax[i].set_ylabel(fr"$P(Q)$")

    ax[i].legend(loc="best")

    arr.append([mean, err])
    arr2.append([mean2, err2])

    print(fr"[L={l}] -> <Q>: {mean}+-{err}, <Q^2>: {mean2}+-{err2}, th (beta={beta[i]}): {sigma**2}")

arr = np.asarray(arr)
arr2 = np.asarray(arr2)

x = np.linspace(1, 10, 1000)
chi_vec = np.vectorize(chi)
y = chi_vec(x)

fig2, ax2 = plt.subplots(figsize=(12,8))
ax2.errorbar(beta[:3], arr2[:3,0]/(L[:3]**2), yerr=arr2[:3,1],fmt='o',capsize=3)
ax2.plot(x,y,'r-.', label=r"$\hat\chi(\beta)_{\text{th}}$")

ax2.set_xlabel(r"$\beta$")
ax2.set_ylabel(r"$\hat\chi=a^2\chi$")

fig2.tight_layout()

fig.tight_layout()
fig.savefig("results/topo_test.pdf",bbox_inches="tight")
fig2.tight_layout()
fig2.savefig("results/topo2.pdf",bbox_inches="tight")
plt.show()

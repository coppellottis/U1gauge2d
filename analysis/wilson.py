import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullLocator
import emcee

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

T = np.arange(1, 9)
R = np.arange(1, 9)

arr = []
arr2 = []

for t in T:
    for r in range(1, t+1):

        W = pd.read_csv(f"data/W_R{r}_T{t}.csv").iloc[:, 0].to_numpy()
        mean = W.mean()
        rho = emcee.autocorr.function_1d(W)
        tau, _ = tau_int(rho, 2)
        err = W.std(ddof=1) * np.sqrt(2*tau/len(W))

        arr.append([r*t, mean, err])

        W2 = pd.read_csv(f"data/Wc_R{r}_T{t}.csv").iloc[:, 0].to_numpy()
        mean2 = W2.mean()
        rho2 = emcee.autocorr.function_1d(W2)
        tau2, _ = tau_int(rho2, 2)
        err2 = W2.std(ddof=1) * np.sqrt(2*tau2/len(W2))

        arr2.append([r*t, mean2, err2])

arr = np.asarray(arr)
arr2 = np.asarray(arr2)

fig, ax = plt.subplots(figsize=(12, 8))

ax.errorbar(
    arr[:, 0], arr[:, 1], arr[:, 2],
    fmt='o', linestyle='none', capsize=4,
    color='blue', label=r"$N=20^2$", alpha=.5
)

ax.errorbar(
    arr2[:, 0], arr2[:, 1], arr2[:, 2],
    fmt='o', linestyle='none', capsize=4,
    color='red', label=r"$N=20^2$, multihit", alpha=.5
)

# Inset
axins = ax.inset_axes([.25, .10, .2, .35])

axins.errorbar(
    arr[:, 0], arr[:, 1], arr[:, 2],
    fmt='o', linestyle='none', capsize=4,
    color='blue', alpha=.5
)

axins.errorbar(
    arr2[:, 0], arr2[:, 1], arr2[:, 2],
    fmt='o', linestyle='none', capsize=4,
    color='red', alpha=.5
)

axins.set_yscale("log")
axins.set_ylim(2.705e-2, 2.75e-2)
axins.set_xlim(9.95, 10.05)

# FIX: esattamente 3 etichette sull'asse y
axins.yaxis.set_major_locator(
    FixedLocator([2.71e-2, 2.7275e-2, 2.745e-2])
)
axins.yaxis.set_minor_locator(NullLocator())

# Zoom connectors
ax.indicate_inset_zoom(axins, edgecolor="black")

# Grafico principale
ax.set_yscale("log")
ax.set_xlim(0, 22)
ax.set_ylim(2e-4, 1e0)

ax.legend(loc="best")
ax.grid(alpha=.2)

ax.set_xlabel(r"$(r\times t)/a^2$")
ax.set_ylabel(r"$\langle W(r/a,t/a)\rangle$")

fig.tight_layout()

fig.savefig("results/wilson_log_N20_b2.pdf", bbox_inches="tight")
plt.show()

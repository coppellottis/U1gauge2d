import numpy as np
import emcee
import pandas as pd
import  matplotlib.pyplot as plt

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


W = pd.read_csv("data/sigma/W_L40_beta12_R3_T3.csv").iloc[:,0].to_numpy()

N = len(W)

rho = emcee.autocorr.function_1d(W)
tau, _ = tau_int(rho, 1000)

print(fr"tau: {tau}, N/tau={(N/tau)}")


fig, ax = plt.subplots(figsize=(12,8))


rho = emcee.autocorr.function_1d(np.array(W))
tau = 0.5 + np.cumsum(rho[1:])
M = np.arange(1, len(rho))
tau_star, M_star = tau_int(rho,150)


ax.plot(M, tau,'-.',label=fr"$W(6,6)$")
ax.scatter(M_star, tau_star,color=ax.lines[-1].get_color(), s=80, zorder=5)

# ###########################################

# ax.set_xlabel(fr"$M$")
# ax.set_ylabel(r"$\tau_{int}^{Q}(M)/\tau_{int}^{Q}(M^\star)$")

# ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

ax.grid(alpha=0.3)

ax.legend(loc="best")

ax.set_xlim(0,1e2)
ax.set_ylim(0,2)


fig.tight_layout()


fig.savefig("results/tau_int_test_1x1.pdf")
plt.show()
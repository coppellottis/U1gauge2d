import numpy as np
import pandas as pd
import  matplotlib.pyplot as plt

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
    
        W = pd.read_csv(f"data/W_R{r}_T{t}.csv").iloc[:,0].to_numpy()
        mean = W.mean()
        err = W.std()

        arr.append([r*t, mean, err])

        W2 = pd.read_csv(f"data/Wc_R{r}_T{t}.csv").iloc[:,0].to_numpy()
        mean2 = W2.mean()
        err2 = W2.std()

        arr2.append([r*t, mean2, err2])

arr = np.asarray(arr)
arr2 = np.asarray(arr2)

fig, ax = plt.subplots(figsize=(12,8))

ax.errorbar(arr[:,0], arr[:,1], arr[:,2], fmt='o', linestyle='none', capsize=4, color='blue', label=fr"$N=20^2$")
ax.errorbar(arr2[:,0], arr2[:,1], arr2[:,2], fmt='o', linestyle='none', capsize=4, color='red', label=fr"$N=20^2$, multihit")

#ax.set_yscale("log")

ax.legend(loc="best")
ax.grid(alpha=.2)

ax.set_xlabel(fr"$(r\times t)/a^2$")
ax.set_ylabel(r"$\langle \mathcal{W}(r/a, t/a)\rangle$")

fig.tight_layout()

fig.savefig("results/wilson.pdf", bbox_inches="tight")
plt.show()

import numpy as np
from scipy.special import ive
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

def sigma_th(x) :
    y = -np.log(ive(np.abs(1), x) / (ive(0, x)))
    return y


print(sigma_th(2))

beta = np.array([2, 4, 6, 8])
sigma = np.array([0.3598527343675133, 0.14685705128809903, 0.0916402075790826, 0.06675908561098748])
err_sigma = np.array([0.0001834604533196477, 8.104047325430019e-05, 7.978733019301913e-05, 0.00018915250672332836])

fig, ax = plt.subplots(figsize=(12,8))
ax.errorbar(beta, sigma, yerr=err_sigma, fmt='o', capsize=3, alpha=.5, label=r"$\hat\sigma_{\text{exp}}$")
x = np.linspace(1, 10,10000)
y = sigma_th(x)
ax.plot(x,y,'k-.', label=r"$\hat\sigma_{\text{th}}$", alpha=.3)

ax.legend(loc="best")
ax.grid(alpha=.3)

ax.set_xlabel(r"$\beta$")
ax.set_ylabel(r"$\hat\sigma$")

fig.tight_layout()

fig.savefig("results/continuum.pdf",bbox_inches="tight")

plt.show()
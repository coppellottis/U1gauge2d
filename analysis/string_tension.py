import numpy as np
import pandas as pd
import  matplotlib.pyplot as plt
import emcee
from scipy.optimize import curve_fit
from scipy.special import ive

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

def format_measurement(mean, error, sig_digits=2):
    """
    Formatta una misura nel formato mean(error).
    sig_digits: numero di cifre significative dell'errore.
    """
    import math

    if error <= 0:
        raise ValueError("L'errore deve essere positivo")

    # Posizione della prima cifra significativa dell'errore
    exponent = math.floor(math.log10(abs(error)))

    # Decimali necessari per l'errore
    decimals = sig_digits - 1 - exponent

    # Arrotondamento dell'errore
    rounded_error = round(error, decimals)

    # Arrotondamento della media alla stessa precisione
    rounded_mean = round(mean, decimals)

    # Formattazione
    mean_str = f"{rounded_mean:.{max(0, decimals)}f}"

    # Errore espresso nelle ultime cifre della media
    error_digits = int(round(rounded_error * 10**decimals))

    return f"{mean_str}({error_digits})"

def lin(x, a, b):
    return a*x+b

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

fig, ax = plt.subplots(figsize=(12,8))
beta = 2

axins = ax.inset_axes([.65, .10, .2, .35])

T = np.arange(1, 6)
R = np.arange(1, 6)

for t in T:
    arr = []
    for r in R:
        if(l == 40 && r>4)
    
        W = pd.read_csv(f"data/Wc_R{r}_T{t}.csv").iloc[:,0].to_numpy()

        mean = W.mean()
        rho = emcee.autocorr.function_1d(W)
        tau, _ = tau_int(rho, 2)
        err = W.std(ddof=1) * np.sqrt(2*tau/len(W))

        V = -np.log(mean)/t
        V_err = err / (t * abs(mean))

        arr.append([r, V, V_err])

    arr = np.asarray(arr)
    ax.errorbar(arr[:, 0],arr[:, 1], yerr=arr[:, 2],fmt='o',capsize=3,label=fr"$t/a={t}$",alpha=0.6)
    shift = np.linspace(-0.0005, 0.0005, len(T))
    x_shifted = arr[:, 0] + shift[t-1]
    axins.errorbar(x_shifted,arr[:, 1], yerr=arr[:, 2],fmt='o',capsize=3, alpha=.6)
    # FIT in R per T=5
    if(t == 4):
        popt, pcov = curve_fit(lin, arr[:,0], arr[:,1], sigma=arr[:,2], absolute_sigma=True)
        a, b = popt
        err_a, err_b = np.sqrt(np.diag(pcov))

        x = np.linspace(0, max(arr[:,0]))
        y = lin(x, a, b)

        formatted = format_measurement(a, err_a)
        expected = -np.log( ive(np.abs(1), beta) / ive(0, beta))

        # Chi quadro e gradi di libertà
        residuals = arr[:,1] - lin(arr[:,0], a, b)
        chi2 = np.sum((residuals / arr[:,2])**2)
        dof = len(arr) - len(popt)
        chi2_r = chi2 / dof

        ax.plot(x, y, 'r-', label=fr"Fit: $\hat\sigma={formatted}$")

        print(f"Exp. b = {formatted}")
        print(f"Th. b  = {expected}")
        print(f"dof    = {dof}")
        print(f"chi²   = {chi2:.4f}")
        print(f"chi²_r = {chi2_r:.4f}")


axins.set_xlim(0.999,1.001)
axins.set_ylim(0.3590, 0.3605)
axins.yaxis.tick_right()
axins.tick_params(axis='x', bottom=False, top=False, labelbottom=False, labeltop=False)
ax.indicate_inset_zoom(axins, edgecolor="black")

ax.set_xlabel(r"$r/a$")
ax.set_ylabel(r"$aV(r)$")

ax.grid(alpha=.2)

ax.legend(loc="best")


fig.savefig("results/string_tension20x20-c.pdf", bbox_inches="tight")
plt.show()

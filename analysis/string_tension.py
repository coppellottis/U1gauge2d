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

for r in R:
    for t in R:
    
        W = pd.read_csv(f"data/Wc_R{r}_T{t}.csv").iloc[:,0].to_numpy()
        mean = W.mean()
        err = W.std()
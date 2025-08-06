#!python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Load data
p = np.load("pgrid.npy")
T = np.load("Tgrid.npy")

seedings = [
    "ind_O3_1000.npy", "ind_O_1000.npy",
    "ind_O3_10000.npy", "ind_O_10000.npy",
    "ind_none.npy"
]
labels = [
    "1000 PPM of O3", "1000 PPM of O",
    "10,000 PPM of O3", "10,000 PPM of O",
    "No O or O3"
]

# Load data
datas4 = [np.load(f) for f in seedings[:4]]
Z5 = np.load(seedings[4])  # denominator for normalization

# Compute vmin/vmax after normalization
normalized = [Z5 / Z   for Z in datas4]
vmin = 1.0
vmax = [5,10,5,10]
fig1, axes = plt.subplots(2, 2, figsize=(10, 10), sharex=True, sharey=True, constrained_layout=True)

# Plot each
mappables = []
for ax, Z, lab, vmaxi in zip(axes.flatten(), normalized, labels[:4], vmax):
    m = ax.contourf(T, p, Z, vmin=vmin, vmax=vmaxi, levels=50)
    ax.set_title(lab)
    ax.set_xlabel("Temperature")
    ax.set_ylabel("Pressure")
    ax.set_title(f"1 to {vmaxi}")


plt.show()

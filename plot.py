import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

DPI = 250

def FixPlot(lx: float, ly: float):
    from matplotlib import rcParams, cycler
    rcParams['font.family'] = 'serif'
    rcParams['font.serif'] = ['Computer Modern']
    rcParams['text.usetex'] = True
    rcParams['font.size'] = 28
    rcParams['axes.linewidth'] = 1.1
    rcParams['axes.labelpad'] = 10.0
    plot_color_cycle = cycler('color', ['000000', 'FE0000', '0000FE', '008001', 'FD8000', '8c564b',
                                        'e377c2', '7f7f7f', 'bcbd22', '17becf'])
    rcParams['axes.prop_cycle'] = plot_color_cycle
    rcParams['axes.xmargin'] = 0
    rcParams['axes.ymargin'] = 0
    rcParams['legend.fancybox'] = False
    rcParams['legend.framealpha'] = 1.0
    rcParams['legend.edgecolor'] = "black"
    rcParams['legend.fontsize'] = 22
    rcParams['xtick.labelsize'] = 22
    rcParams['ytick.labelsize'] = 22

    rcParams['ytick.right'] = True
    rcParams['xtick.top'] = True

    rcParams['xtick.direction'] = "in"
    rcParams['ytick.direction'] = "in"
    rcParams['axes.formatter.useoffset'] = False

    rcParams.update({"figure.figsize": (lx, ly),
                    "figure.subplot.left": 0.177, "figure.subplot.right": 0.946,
                     "figure.subplot.bottom": 0.156, "figure.subplot.top": 0.965,
                     "axes.autolimit_mode": "round_numbers",
                     "xtick.major.size": 7,
                     "xtick.minor.size": 3.5,
                     "xtick.major.width": 1.1,
                     "xtick.minor.width": 1.1,
                     "xtick.major.pad": 5,
                     "xtick.minor.visible": True,
                     "ytick.major.size": 7,
                     "ytick.minor.size": 3.5,
                     "ytick.major.width": 1.1,
                     "ytick.minor.width": 1.1,
                     "ytick.major.pad": 5,
                     "ytick.minor.visible": True,
                     "lines.markersize": 10,
                     "lines.markeredgewidth": 0.8, 
                     "mathtext.fontset": "cm"})

def FixPlot_(lx: float, ly: float):
    from matplotlib import rcParams, cycler
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
    rcParams['text.usetex'] = False
    rcParams['font.size'] = 28
    rcParams['axes.linewidth'] = 1.1
    rcParams['axes.labelpad'] = 10.0
    plot_color_cycle = cycler('color', ['000000', 'FE0000', '0000FE', '008001', 'FD8000', '8c564b',
                                        'e377c2', '7f7f7f', 'bcbd22', '17becf'])
    rcParams['axes.prop_cycle'] = plot_color_cycle
    rcParams['axes.xmargin'] = 0
    rcParams['axes.ymargin'] = 0
    rcParams['legend.fancybox'] = False
    rcParams['legend.framealpha'] = 1.0
    rcParams['legend.edgecolor'] = "black"
    rcParams['legend.fontsize'] = 28
    rcParams['xtick.labelsize'] = 22
    rcParams['ytick.labelsize'] = 22

    rcParams['ytick.right'] = True
    rcParams['xtick.top'] = True

    rcParams['xtick.direction'] = "in"
    rcParams['ytick.direction'] = "in"

    rcParams.update({"figure.figsize": (lx, ly),
                    "figure.subplot.left": 0.177, "figure.subplot.right": 0.946,
                     "figure.subplot.bottom": 0.156, "figure.subplot.top": 0.965,
                     "axes.autolimit_mode": "round_numbers",
                     "xtick.major.size": 7,
                     "xtick.minor.size": 3.5,
                     "xtick.major.width": 1.1,
                     "xtick.minor.width": 1.1,
                     "xtick.major.pad": 5,
                     "xtick.minor.visible": True,
                     "ytick.major.size": 7,
                     "ytick.minor.size": 3.5,
                     "ytick.major.width": 1.1,
                     "ytick.minor.width": 1.1,
                     "ytick.major.pad": 5,
                     "ytick.minor.visible": True,
                     "lines.markersize": 10,
                     "lines.markeredgewidth": 0.8, 
                     "mathtext.fontset": "dejavusans"}) #"cm"

FixPlot_(8, 8)
data = pd.read_table("./data.dat", header=None, sep="\t")

fig, ax = plt.subplots()
ax.scatter(data[0], data[1], s=5.0)
ax.set_xlim([min(data[0]), max(data[0])])
ax.set_ylim([min(data[1]), max(data[1])])
ax.set_xlabel("T")
ax.set_ylabel("E")
fig.savefig("./energy_vs_T.png", dpi=DPI, facecolor="white", bbox_inches="tight")

fig, ax = plt.subplots()
ax.scatter(data[0], np.abs(data[2]), s=5.0)
ax.set_xlim([min(data[0]), max(data[0])])
ax.set_ylim([min(np.abs(data[2])), max(np.abs(data[2]))])
ax.set_xlabel("T")
ax.set_ylabel("|M|")
fig.savefig("./magnetization_vs_T.png", dpi=DPI, facecolor="white", bbox_inches="tight")

fig, ax = plt.subplots()
ax.scatter(data[0], data[3], s=5.0)
ax.set_xlim([min(data[0]), max(data[0])])
ax.set_ylim([min(data[3]), max(data[3])])
ax.set_xlabel("T")
ax.set_ylabel("C")
fig.savefig("./specific_heat_vs_T.png", dpi=DPI, facecolor="white", bbox_inches="tight")

fig, ax = plt.subplots()
ax.scatter(data[0], data[4], s=5.0)
ax.set_xlim([min(data[0]), max(data[0])])
ax.set_ylim([min(data[4]), max(data[4])])
ax.set_xlabel("T")
ax.set_ylabel("$\\chi$")
fig.savefig("./magnetic_sus_vs_T.png", dpi=DPI, facecolor="white", bbox_inches="tight")

fig, ax = plt.subplots()
d2MdT2 = np.diff(data[4], n=2)
d2MdT2 = np.concatenate([d2MdT2, [0, 0]])
ax.scatter(data[0], d2MdT2, s=5.0)
ax.set_xlabel("T")
ax.set_ylabel("d$^2$M/dT$^2$")
fig.savefig("./d2MdT2.png", dpi=DPI, facecolor="white", bbox_inches="tight")

fig, ax = plt.subplots()
d2EdT2 = np.diff(data[1], n=2)
d2EdT2 = np.concatenate([d2EdT2, [0, 0]])
ax.scatter(data[0], d2EdT2, s=5.0)
ax.set_xlabel("T")
ax.set_ylabel("d$^2$E/dT$^2$")
fig.savefig("./d2EdT2.png", dpi=DPI, facecolor="white", bbox_inches="tight")
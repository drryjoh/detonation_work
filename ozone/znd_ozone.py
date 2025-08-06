#!/usr/bin/env python3
from sdtoolbox.postshock import CJspeed, PostShock_fr
from sdtoolbox.znd import zndsolve
from sdtoolbox.cv import cvsolve
from sdtoolbox.utilities import CJspeed_plot, znd_plot
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from inspect import signature
from values import *
from ppplot import *



def get_znd_solution(P1, T1, X):
    q = X
    mech_yaml = 'ffcm_h2_o3.yaml'

    cj_speed,R2,plot_data = CJspeed(P1,T1,q,mech_yaml,fullOutput=True)
    gas1 = ct.Solution(mech_yaml)
    gas1.TPX = T1,P1,q

    gas = PostShock_fr(cj_speed, P1, T1, q, mech_yaml)
    print(gas.TP)
    print(cj_speed)

    out = zndsolve(gas,gas1,cj_speed,t_end=2e-4,max_step=0.0001,relTol=1e-8,advanced_output=True)


    T = out['T']
    x = out['distance']
    P = out['P']
    U = out['U']

    # Include pre-shock state

    # breakpoint()

    P_stag = P*0.

    # mm
    x *= 1.e3
    # kPa
    P /= 1.e3
    P_stag /= 1.e3

    # SPR
    SPR = P_stag/np.max(P_stag)

    return cj_speed, x, T, P, U, P_stag, SPR

def induction_location(x, T):
    dTdx = np.diff(T)/np.diff(x)
    idx  = np.argmin(np.abs(T - (T[0]+100)))
    return [x[idx], T[idx]]

save_dir = ''
prepare_plot(linewidth=1.5, fontsize=24, markersize=12)

try:
    plt.style.use('seaborn-deep')
except OSError:
    plt.style.use('seaborn-v0_8')

label_01 = 'No ozone'
label_02 = '1000 ppm'
label_03 = '10,000 ppm'

style_01 = '-'
style_02 = '--'
style_03 = '-.'

P1 = 100e3
T1 = 700

# No ozone
#O3:0.054949494949495
n_znds = 10
ozones = np.linspace(0, 1, n_znds)**2 * 0.054949494949495 * 2  
cjs = []
cjs_O = []
ls_O3 = []
ls_O  = []
for ozone in ozones:
    cj, x, T, p, u, p_stage, spr = get_znd_solution(P1, T1, f"H2:0.68 O2:1 N2:3.76 O3:{ozone}")
    [x_ind, T_ind] = induction_location(x, T)
    ls_O3.append(x_ind)
    cjs.append(cj)

    cj, x, T, p, u, p_stage, spr = get_znd_solution(P1, T1, f"H2:0.68 O2:1 N2:3.76 O:{ozone}")
    [x_ind, T_ind] = induction_location(x, T)
    ls_O.append(x_ind)
    cjs_O.append(cj)


# Convert to arrays if not already
ls_O3 = np.array(ls_O3)
ls_O = np.array(ls_O)

cjs = np.array(cjs)
cjs_O = np.array(cjs_O)

ozone_ppm = ozones / 0.054949494949495 * 10000


#reduction_O3 = ls_O3[0] / ls_O3[1:]
#reduction_O = ls_O[0] / ls_O[1:]
reduction_O3 = ls_O3
reduction_O = ls_O
fig, ax1 = plt.subplots()
plot_O = False
ax1.plot(ozone_ppm, reduction_O3, 'o-', label="Induction Length")
if plot_O:
    ax1.plot(ozone_ppm, reduction_O, 'o-', label="Induction Length, O")
#ax1.plot(ozone_ppm, cjs[1:]/cjs[0], 's-', label="CJ Velocity")
#if plot_O:
#    ax1.plot(ozone_ppm, cjs_O[1:]/cjs_O[0], 's-', label="CJ Velocity, O")
#ax1.plot(ozone_ppm, np.ones_like(reduction_O3), '--')
if plot_O:
    name = ""
else:
    name = "O3 "
ax1.set_xlabel(f"{name}Concentration (PPM)")
if plot_O:
    name  = "/O"
else:
    name = ""
ax1.set_ylabel(f"Induction Length (mm)")

ax1.set_xticks([0,1000, 10000, 20000])
plt.legend()
plt.tight_layout()
plt.show()

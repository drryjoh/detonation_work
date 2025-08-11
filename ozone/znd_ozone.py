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



def get_znd_solution(P1, T1, X, name  = 'ffcm_h2_o3.yaml'):
    q = X
    mech_yaml = name

    cj_speed,R2,plot_data = CJspeed(P1,T1,q,mech_yaml,fullOutput=True)
    gas1 = ct.Solution(mech_yaml)
    gas1.TPX = T1,P1,q

    gas = PostShock_fr(cj_speed, P1, T1, q, mech_yaml)

    out = zndsolve(gas,gas1,cj_speed,t_end=2e-6,max_step=1e-3,relTol=1e-8,advanced_output=True)


    T = out['T']
    x = out['distance']
    P = out['P']
    U = out['U']

    P_stag = P*0.
    x *= 1.e3
    P /= 1.e3
    P_stag /= 1.e3
    SPR = P_stag/np.max(P_stag)

    return cj_speed, x, T, P, U, P_stag, SPR

def induction_location(x, T):
    idx  = np.argmin(np.abs(T - (T[0] + 100)))
    return [x[idx], T[idx]]

def get_induction_time_it(gas, T, p, X_string, time_resolution = 1e-9, n_steps = 10000):
    gas.TPX = T, p, X_string
    reactor = ct.IdealGasMoleReactor(gas)
    network = ct.ReactorNet([reactor])
    time = np.linspace(0, n_steps * time_resolution, n_steps)
    
    T_start =  T
    for t in time:
        network.advance(t)
        if reactor.T-T_start>20: 
            return [True, t, reactor.T]

    return [False, 0, 0]

def induction_location_2(T, p, X_string, u, name = 'ffcm_h2_o3.yaml'):
    gas = ct.Solution(name)
    p  *= 1.e3
    [found, time, temperature] = get_induction_time_it(gas, T, p, X_string)

    if found:
        return [u[0] * time * 1000, temperature] #in mm
    else:
        exit("Did not find inflection point")

    

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
n_znds = 30
ozones = np.linspace(0, 1, n_znds)**1.5 * 0.054949494949495

cjs = []
cjs_full = []
cjs_O = []
cjs_lean = []

ls_O3 = []
ls_O3_full = []
ls_O  = []
ls_lean = []

plot_O = False
plot_lean = False
plot_princeton = True
diluent  = "N2:3.76"

for k, ozone in enumerate(ozones):
    print(f" Calculating {k+1}th znd solution, with O3 = {ozone}")
    cj, x, T, p, u, p_stage, spr = get_znd_solution(P1, T1, f"H2:0.68 O2:1 {diluent} O3:{ozone}")
    [x_ind, T_ind] = induction_location_2(T[0], p[0], f"H2:0.68 O2:1 {diluent} O3:{ozone}", u)

    ls_O3.append(x_ind)
    cjs.append(cj)
    if plot_O:
        cj, x, T, p, u, p_stage, spr = get_znd_solution(P1, T1, f"H2:0.68 O2:1 {diluent} O:{ozone}")
        [x_ind, T_ind] = induction_location_2(T[0], p[0], f"H2:0.68 O2:1 {diluent} O:{ozone}", u)
        ls_O.append(x_ind)
        cjs_O.append(cj)
    if plot_lean:
        cj, x, T, p, u, p_stage, spr = get_znd_solution(P1, T1, f"H2:0.68 O2:{1.0+ozone*2/3} {diluent} O:0")
        [x_ind, T_ind] = induction_location_2(T[0], p[0], f"H2:0.68 O2:{1.0+ozone*2/3} {diluent} O:0", u)
        ls_lean.append(x_ind)
        cjs_lean.append(cj)
    if plot_princeton:
        cj, x, T, p, u, p_stage, spr = get_znd_solution(P1, T1, f"H2:0.68 O2:1 {diluent} O3:{ozone}", name = "ffcm_h2_o3_full.yaml")
        [x_ind, T_ind] = induction_location_2(T[0], p[0], f"H2:0.68 O2:1 {diluent} O3:{ozone}",u, name = "ffcm_h2_o3_full.yaml")
        ls_O3_full.append(x_ind)
        cjs_full.append(cj)       


# Convert to arrays if not already
ls_O3 = np.array(ls_O3)
cjs = np.array(cjs)

ls_O3_full = np.array(ls_O3_full)
cjs_full = np.array(cjs)

ls_O = np.array(ls_O)
cjs_O = np.array(cjs_O)

ozone_ppm = ozones / 0.054949494949495 * 10000

#plot figure
fig, ax1 = plt.subplots()
colors_ax1 = plt.cm.tab10.colors[:2]
colors_ax2 = plt.cm.tab10.colors[2:4]

if plot_O:
    l1 = ax1.plot(ozone_ppm, ls_O3, '-s', label="Induction Length, O3", color=colors_ax1[0])
    l2 = ax1.plot(ozone_ppm, ls_O, '-d', label="Induction Length, O", color=colors_ax1[1])
else:
    l1 = ax1.plot(ozone_ppm, ls_O3/ls_O3[0], '-', label="Induction Length", color=colors_ax1[0])
    l2 = []

ax2 = ax1.twinx()

if plot_O:
    l3 = ax2.plot(ozone_ppm, cjs, "-o", label="CJ Velocity, O3", color=colors_ax2[0])
    l4 = ax2.plot(ozone_ppm, cjs_O, "-^", label="CJ Velocity, O", color=colors_ax2[1])
else:
    l3 = ax2.plot(ozone_ppm, cjs, "-", label="CJ Velocity", color=colors_ax2[0])
    l4 = []

# Labels
ax1.set_xlabel(f"{'O3 ' if not plot_O else ''}Concentration (PPM)")
ax1.set_ylabel("Induction Length (mm)")
ax2.set_ylabel("CJ Velocity")

if plot_lean:
    l6 = ax2.plot(ozone_ppm, cjs_lean, "--k", label="CJ Velocity, O2 increase")
    
if plot_princeton:
    l7 = ax1.plot(ozone_ppm, ls_O3_full/ls_O3_full[0], 'ok', label="Induction Length, full")

lines = l1 + l2 + l3 + l4
if plot_lean:
    lines += l6
if plot_princeton:
    lines += l7
labels = [line.get_label() for line in lines]
ax1.legend(lines, labels, loc='upper center', ncol=2)

ax1.set_xticks([0, 1000, 10000])
ax1.set_ylim([0,1])
ax2.set_ylim([cjs[0] - cjs[0]*0.002, cjs[0] + cjs[0]*0.1])
ax2.set_yticks([cjs[0], cjs[0] + cjs[0]*0.05,cjs[0] + cjs[0]*0.1])
ax2.set_yticklabels(["$u_{cj,0}$", "$1.05\\times u_{cj,0}$","$1.1\\times u_{cj,0}$"])
plt.tight_layout()
print(f"CJ Velocity PCT Change of O3 compared to none: {cjs[-1]/cjs[0]*100-100}%")
print(f"Induction distance of O3 compared to none: {ls_O3[-1]/ls_O3[0]*100-100}%")

if plot_O:
    print(f"Induction distance of O compared to none:  {ls_O[-1]/ls_O3[0]*100-100}%")
    print(f"CJ Velocity PCT Change of O compared to none:  {cjs_O[-1]/cjs_O[0]*100-100}%")
plt.show()

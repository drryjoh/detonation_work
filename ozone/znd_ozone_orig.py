"""
Shock and Detonation Toolbox Demo Program

Computes ZND and CV models of detonation with the shock front
traveling at the CJ speed.  Evaluates various measures of the reaction
zone thickness and exothermic pulse width, effective activation energy
and Ng stability parameter.

################################################################################
Theory, numerical methods and applications are described in the following report:

    Numerical Solution Methods for Shock and Detonation Jump Conditions, S.
    Browne, J. Ziegler, and J. E. Shepherd, GALCIT Report FM2006.006 - R3,
    California Institute of Technology Revised August, 2017

Please cite this report and the website if you use these routines.

Please refer to LICENCE.txt or the above report for copyright and disclaimers.

http://shepherd.caltech.edu/EDL/public/sdt/SD_Toolbox/

################################################################################
Updated August 2018
Tested with:
    Python 3.5 and 3.6, Cantera 2.3 and 2.4
Under these operating systems:
    Windows 8.1, Windows 10, Linux (Debian 9)
"""
from sdtoolbox.postshock import CJspeed, PostShock_fr
from sdtoolbox.znd import zndsolve
from sdtoolbox.cv import cvsolve
from sdtoolbox.utilities import CJspeed_plot, znd_plot
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from inspect import signature

def inv(a): return 1/a
def sqr(a): return a*a
def magsqr(a):
    b = []
    for ai in a:
        b.append(np.sum(ai*ai))
    return np.array(b)

def sqrt(a): return np.sqrt(a)

def R(): return 8.314459800e+03
def universal_gas_constant(): return R()

def reference_length():  return 1;

def reference_density():  return 1;

def reference_pressure():  return 101325.0;

def reference_temperature():  return 1000.0;

def reference_viscosity():  return 2e-5;

def reference_diffusion():  return 3e-5;

def reference_conductivity():  return 0.04;

def reference_specific_heat_volume(): return 715.0

def reference_concentration():  return reference_pressure()/reference_temperature()/R();

def reference_reynolds():  return reference_density() * reference_velocity() * reference_length() * inv(reference_viscosity());

def reference_prandtl():  return reference_specific_heat_pressure() * reference_viscosity() * inv(reference_conductivity());

def reference_schmidts():  return reference_viscosity() * inv(reference_density() * reference_diffusion());

def reference_specific_heat_pressure():  return reference_energy() * inv(reference_temperature());

def reference_energy():  return reference_velocity_sqr();

def reference_time():  return reference_length() * inv(reference_velocity());

def reference_velocity():  return sqrt((reference_pressure()) * inv(reference_density()));

def reference_velocity_sqr():  return sqr(reference_velocity());

def reference_molecular_weight():  return reference_density() * inv(reference_concentration());

def dimensionalize_energy(energy):
    return energy*reference_energy()

def dimensionalize_concentrations(species):
    return species*reference_concentration()
def nondimensionalize_concentrations(species):
    return species/reference_concentration()

def dimensionalize_momentum(momentum):
    return momentum*reference_velocity()*reference_density()

def internal_energy(energy, momentum, density):
    return energy - 0.5*magsqr(momentum)/density

def scale(a, b):
    return a*b
def multiply(a, b, c):
    return a*b*c
def tuple_weight(a, b):
    return a*b

def inv(a): return 1./a

# AR HE N2 H2 H O O2 OH H2O HO2 H2O2 O3
def molecular_weights_dim():
    return np.array([39.950000000000003,
                      4.0026020000000004,
                      28.013999999999999,
                      2.016,
                      1.008,
                      15.999000000000001,
                      31.998000000000001,
                      17.007000000000001,
                      18.015000000000001,
                      33.006,
                      34.014000000000003,
                      47.997])

def inv_molecular_weights_dim(): return inv(molecular_weights_dim())

def contract(a, b): return np.dot(a, b)

def average_molecular_weight_massfractions_dim(massfractions):
    # breakpoint()
    return inv(contract(massfractions,
                        inv_molecular_weights_dim()))

def concentrations_from_mass_fractions(massfractions, p_dim, temperature_dim):
    return scale(inv(universal_gas_constant()),
                  scale(multiply(average_molecular_weight_massfractions_dim(massfractions),
                                 p_dim,
                                 inv(temperature_dim)),
                        tuple_weight(inv(molecular_weights_dim()),
                                     massfractions)))

def concentrations_from_mole_fractions(X, p_dim, temperature_dim):
    return scale(inv(scale(universal_gas_constant(),
                           temperature_dim)),
                 scale(p_dim,
                       X))

def prepare_plot(reset=False, defaults=False, close_all=True, fontsize=12.,
        font={'family':'serif', 'serif': ['DejaVu Sans']}, linewidth=1.5,
        markersize=4.0, axis=None, cmap='cividis', equal_AR=False):
    '''
    This function sets parameters for plotting.

    Inputs:
    -------
        reset: if True, will reset to default parameters before setting
            input arguments
        defaults: if True, will reset to default parameters and then
            immediately return (input arguments not set)
        close_all: if True, will close all current figures
        fontsize: font size
        font: font
        linewidth: line width
        markersize: size of markers
        axis: axis limits [2*ndims]
        cmap: colormap
        equal_AR: if True, will set equal aspect ratio (only affects 2D)
    '''
    if reset or defaults:
        mpl.rcdefaults() # return to default settings
        if defaults:
            return

    # Use tex syntax
    plt.rc('text', usetex=True)

    # Set parameters
    if close_all:
        plt.close("all")
    mpl.rcParams['font.size'] = fontsize
    plt.rc('font',**font)
    # font={'family':'serif', 'serif': ['computer modern roman']}
    mpl.rcParams['lines.linewidth'] = linewidth
    mpl.rcParams['lines.markersize'] = markersize
    # Note: cividis is more colorblind-friendly than viridis (which is already pretty colorblind-friendly)
    mpl.rcParams['image.cmap'] = cmap
    # Color cycle
    # plt.style.use('tableau-colorblind10')
    plt.style.use('seaborn-colorblind')

    if axis is not None:
        plt.axis(axis)
    if equal_AR:
        plt.gca().set_aspect('equal', adjustable='box')


def save_figure(file_name='fig', file_type='pdf', crop_level=1):
    '''
    This function saves the current figure to disk.

    Inputs:
    -------
        file_name: name of file to save
        file_type: type of file to save
        crop_level: 0 for no cropping, 1 for some cropping, 2 for a lot
            of cropping
    '''
    file_name = file_name + '.' + file_type
    if crop_level == 0:
        # Don't crop
        plt.savefig(file_name)
    elif crop_level == 1:
        # Crop a little
        plt.savefig(file_name, bbox_inches='tight')
    elif crop_level == 2:
        # Crop a lot
        plt.savefig(file_name, bbox_inches='tight', pad_inches=0.0)
    else:
        raise ValueError


get_SPR = True
if get_SPR:
    from jenre import *
    api = API('post')

def get_znd_solution(P1, T1, X):
    # X: mole ratios (unnormalized mole fractions)
    q = X
    # mech = 'hydrogen_34.cti'
    mech_yaml = 'ffcm_h2_o3.yaml'
    # fname = 'hydrogen_34'

    # P1 = 6670;
    # #10132.5;
    # T1 = 300
    # q = 'H2:2 O2:1 AR:7'

    # Find CJ speed and related data, make CJ diagnostic plots
    cj_speed,R2,plot_data = CJspeed(P1,T1,q,mech_yaml,fullOutput=True)
    # cj_speed = 1617.5
    # sig = signature(CJspeed)
    # print(sig)
    #CJspeed_plot(plot_data,cj_speed)

    # Set up gas object
    gas1 = ct.Solution(mech_yaml)
    gas1.TPX = T1,P1,q

    # Find post shock state for given speed
    # cj_speed *= 1.5 # overdriven
    gas = PostShock_fr(cj_speed, P1, T1, q, mech_yaml)
    print(gas.TP)
    print(cj_speed)
    # Solve ZND ODEs, make ZND plots

    out = zndsolve(gas,gas1,cj_speed,t_end=2e-4,max_step=0.0001,relTol=1e-8,advanced_output=True)
    # gas = PostShock_fr(cj_speed, P1, T1, q, mech_yaml)
    # print(gas.TP)
    #znd_plot(out,maxx=0.1,
    #         major_species=['H2', 'O2', 'H2O'],
    #         minor_species=['H', 'O', 'OH', 'H2O2', 'HO2'])
    # plt.plot(out['distance'],out['T'],'ok',mfc='white')
    # plt.show()
    # np.save('Y.npy',out['species'])
    # np.save('T.npy',out['T'])
    # np.save('P.npy',out['P'])
    # np.save('x.npy',out['distance'])
    # np.save('U.npy',out['U'])
    # np.save('U1.npy',out['U1'])
    # print(out['U1'])
    # print(np.max(out['T']))

    T = out['T']
    x = out['distance']
    P = out['P']
    U = out['U']

    # Include pre-shock state
    x = np.concatenate([[-1.e-3, 0.],x])
    T = np.concatenate([[T1, T1],T])
    P = np.concatenate([[P1, P1],P])
    U = np.concatenate([[cj_speed, cj_speed],U])

    # breakpoint()

    P_stag = P*0.

    if get_SPR:
        Y = out['species']

        for i, xi in enumerate(x):
            Pi = P[i]
            Ti = T[i]
            Ui = U[i]
            if i < 2:
                # left of shock
                C = concentrations_from_mole_fractions(gas1.X, Pi, Ti)
            else:
                Yi = Y[:,i-2]
                # breakpoint()
                C = concentrations_from_mass_fractions(Yi, Pi, Ti)

            P_stag[i] = api.scalar_function('PressureStagnation','(ConcentrationsVelocityXPressureTemperature)', np.append(C, [Ui, Pi, Ti]))

    # mm
    x *= 1.e3
    # kPa
    P /= 1.e3
    P_stag /= 1.e3

    # SPR
    SPR = P_stag/np.max(P_stag)

    return cj_speed, x, T, P, U, P_stag, SPR


save_dir = ''
prepare_plot(linewidth=1.5, fontsize=24, markersize=12)

plt.style.use('seaborn-deep') # https://matplotlib.org/3.4.3/gallery/style_sheets/style_sheets_reference.html

label_01 = 'No ozone'
label_02 = '1000 ppm'
label_03 = '10,000 ppm'

style_01 = '-'
style_02 = '--'
style_03 = '-.'

P1 = 100e3
T1 = 700

# No ozone
cj_speed_01, x_01, T_01, P_01, U_01, P_stag_01, SPR_01 = get_znd_solution(P1, T1, 'H2:0.68 O2:1 N2:3.76')

# 1000 ppm
cj_speed_02, x_02, T_02, P_02, U_02, P_stag_02, SPR_02 = get_znd_solution(P1, T1, 'H2:0.68 O2:1 N2:3.76 O3:0.005445445445445')

# 10,000 ppm
cj_speed_03, x_03, T_03, P_03, U_03, P_stag_03, SPR_03 = get_znd_solution(P1, T1, 'H2:0.68 O2:1 N2:3.76 O3:0.054949494949495')

x_start = np.min(x_01)
x_end = 10 # mm

# T
plt.figure()
plt.plot(x_01, T_01, style_01, label=label_01)
plt.plot(x_02, T_02, style_02, label=label_02)
plt.plot(x_03, T_03, style_03, label=label_03)
# plt.plot(x_04, T_04, style_04, label=label_04)
plt.xlabel('$x_1$ (mm)')
plt.ylabel('$T$ (K)')
plt.xlim(x_start, x_end)
plt.legend(loc='best')
save_figure(save_dir + 'ZND_T')

# P
plt.figure()
plt.plot(x_01, P_01, style_01, label=label_01)
plt.plot(x_02, P_02, style_02, label=label_02)
plt.plot(x_03, P_03, style_03, label=label_03)
# plt.plot(x_04, T_04, style_04, label=label_04)
plt.xlabel('$x_1$ (mm)')
plt.ylabel('$P$ (kPa)')
plt.xlim(x_start, x_end)
plt.legend(loc='best')
save_figure(save_dir + 'ZND_P')

# SPR
if get_SPR:
    plt.figure()
    plt.plot(x_01, SPR_01, style_01, label=label_01)
    plt.plot(x_02, SPR_02, style_02, label=label_02)
    plt.plot(x_03, SPR_03, style_03, label=label_03)
    # plt.plot(x_04, T_04, style_04, label=label_04)
    plt.xlabel('$x_1$ (mm)')
    plt.ylabel('SPR')
    plt.xlim(0., x_end)
    plt.ylim(0.361, 0.44)
    plt.legend(loc='best')
    save_figure(save_dir + 'ZND_SPR')

# # P_stag
# plt.figure()
# plt.plot(x_01, P_stag_01, style_01, label=label_01)
# plt.plot(x_02, P_stag_02, style_02, label=label_02)
# plt.plot(x_03, P_stag_03, style_03, label=label_03)
# # plt.plot(x_04, T_04, style_04, label=label_04)
# plt.xlabel('$x_1$ (mm)')
# plt.ylabel('$P_s$ (kPa)')
# plt.xlim(x_start, x_end)
# plt.legend(loc='best')
# save_figure(save_dir + 'ZND_U')

# U
# plt.figure()
# plt.plot(x_01, U_01, style_01, label=label_01)
# plt.plot(x_02, U_02, style_02, label=label_02)
# plt.plot(x_03, U_03, style_03, label=label_03)
# # plt.plot(x_04, T_04, style_04, label=label_04)
# plt.xlabel('$x_1$ (mm)')
# plt.ylabel('$v$ (m/s)')
# # plt.xlim(x_start, x_end)
# plt.legend(loc='best')
# save_figure(save_dir + 'ZND_U')

plt.show()





# breakpoint()
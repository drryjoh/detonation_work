#!python3
"""
No ozone
            MoleFractions
              H2          0.68
              O2          1
              N2          3.76
 
1000 ppm ozone by mole
            MoleFractions
            H2          0.68
            O2          1
            N2          3.76
            O3          0.005445445445445
 
 
10,000 ppm ozone by mole
            MoleFractions
            H2          0.68
            O2          1
            N2          3.76
            O3          0.054949494949495
"""
import cantera as ct
import numpy as np

def get_induction_time_it(gas, test_conditions, time_resolution = 1e-9, n_steps = 10000):
    
    gas.TPX = (
    test_conditions["temperature"],
    test_conditions["pressure"],
    test_conditions["molefractions"]
    )

    reactor = ct.IdealGasMoleReactor(gas)
    network = ct.ReactorNet([reactor])
    time = np.linspace(0, n_steps * time_resolution, n_steps)
    
    T_start =  test_conditions["temperature"]
    for t in time:
        network.advance(t)
        if reactor.T-T_start>50:
            return [True, t, reactor.T]

    return [False, 0, 0]

def get_induction_time(gas, test_conditions, time_resolution = 1e-9, n_steps = 10000):
    for i in range(5):
        [found, time, temperature] = get_induction_time_it(gas, test_conditions, time_resolution = time_resolution, n_steps = 10000)
        if found:
            return [time, temperature]
        else:
            time_resolution = 10 * time_resolution
    sys.exit(f"no inductiontime found T = {test_conditions['temperature']}, p = {test_conditions['pressure']}")


gas = ct.Solution("ffcm_h2_o3.yaml")
T_vn = 1550
kpa = 1000
p_vn = 200 * kpa #kpa
ozone = [0.005445445445445, 0.054949494949495]
labels = ["1000", "10,000"] # ppm
n_grid = 10

induction_times = []
Ts = np.linspace(T_vn - 200, T_vn + 200, n_grid)
ps = np.linspace(p_vn - 50 * kpa, p_vn + 200 * kpa, n_grid)
[T, P] = np.meshgrid(Ts, ps)
ind_O3_1000 = T * 0
ind_O_1000  = T * 0
ind_O3_10000 = T * 0
ind_O_10000  = T * 0
ind_none  = T * 0

for i, T_i in enumerate(Ts):
    for j, p_j in enumerate(ps): 
        test_conditions = {
                "temperature": T_i,  # K
                "pressure": p_j,  # Pa
                "molefractions": {
                    "H2": 0.68,
                    "O2": 1,
                    "N2": 3.76
                }
            }
        [induction_time_none, T_none] = get_induction_time(gas, test_conditions, time_resolution = 1e-9, n_steps = 10000)
        ind_none[i,j] = induction_time_none

        #ozone/O effects
        for k, ozone_i in enumerate(ozone):
            test_conditions = {
                "temperature": T_i,  # K
                "pressure": p_j,  # Pa
                "molefractions": {
                    "H2": 0.68,
                    "O2": 1,
                    "N2": 3.76,
                    "O3": ozone_i
                }
            }
            [induction_time_ozone, T_O3] = get_induction_time(gas, test_conditions, time_resolution = 1e-9, n_steps = 10000)

            test_conditions = {
                "temperature": T_i,  # K
                "pressure": p_j,  # Pa
                "molefractions": {
                    "H2": 0.68,
                    "O2": 1,
                    "N2": 3.76,
                    "O": ozone_i
                }
            }
            [induction_time_o, T_o] = get_induction_time(gas, test_conditions, time_resolution = 1e-9, n_steps = 10000)
            if labels[k] == "1000":
                ind_O3_1000[i,j] = induction_time_ozone
                ind_O_1000[i,j]  = induction_time_o
            else:
                ind_O3_10000[i,j] = induction_time_ozone
                ind_O_10000[i,j]  = induction_time_o
np.save("ind_O3_1000.npy", ind_O3_1000)
np.save("ind_O_1000.npy", ind_O_1000)

np.save("ind_O3_10000.npy", ind_O3_10000)
np.save("ind_O_10000.npy", ind_O_10000)

np.save("ind_none.npy", ind_none)

np.save("Tgrid.npy", T)
np.save("pgrid.npy", P)

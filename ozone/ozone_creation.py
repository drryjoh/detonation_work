#!python3
import cantera as ct
import numpy as np
from ambiance import Atmosphere
gas = ct.Solution("ffcm_h2_o3.yaml")
T_spark = 3000
M = 2.0
h = np.array(10000)  # m
atm = Atmosphere(h)
T_alt = atm.temperature  # K
p_alt = atm.pressure     # Pa
n_steps = 100
time_resolution  = 1e-9

test_conditions = {
    "temperature": T_alt,  # K
    "pressure": p_alt,  # Pa
    "molefractions": {
        "O2": 1,
        "N2": 3.76
    }
}


gas.TPX = test_conditions["temperature"], test_conditions["pressure"], test_conditions["molefractions"]

v = M * np.sqrt(gas.T * gas.cp/gas.cv * 287)
print(v)

test_conditions = {
    "temperature": T_spark,  # K
    "pressure": p_alt,  # Pa
    "molefractions": {
        "O2": 1,
        "N2": 3.76
    }
}
gas.TPX = test_conditions["temperature"], test_conditions["pressure"], test_conditions["molefractions"]

reactor = ct.IdealGasMoleReactor(gas)
network = ct.ReactorNet([reactor])
time = np.linspace(0, n_steps * time_resolution, n_steps)

T_start =  test_conditions["temperature"]
XO3 = []
XO = []
X = gas.X
for t in time:
    network.advance(t)
    X = reactor.thermo.X
    XO3.append(X[gas.species_index("O3")])
    XO.append(X[gas.species_index("O")])
print(f"after {time[-1]} seconds:")
print(f"O: {XO[-1]}")
print(f"O3: {XO3[-1]}")


gas.TPX = 300, 50000, X
reactor = ct.IdealGasMoleReactor(gas)
network = ct.ReactorNet([reactor])
time = np.linspace(0, 1e-4, n_steps)
for t in time:
    network.advance(t)
    X = reactor.thermo.X
    XO3.append(X[gas.species_index("O3")])
    XO.append(X[gas.species_index("O")])
print(f"after {time[-1]} seconds:")
print(f"O: {XO[-1]}")
print(f"O3: {XO3[-1]}")




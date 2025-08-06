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


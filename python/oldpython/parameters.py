import eos_importer
import importlib
importlib.reload(eos_importer)

### Parameters for the integration of the profiles (scalar field, pressure, ...)

# min and max radius for the integration of the differential equations
minIntegrationRadius = 1e-10
maxIntegrationRadius = 200

### parameters specific to the model we are considering

mu = 1 # scalar field mass parameter
lam = 0 # scalar field quartic self-interaction parameter

# EoS_v = eos_importer.EoS_polytrope()
EoS_v = eos_importer.EoS("DD2")
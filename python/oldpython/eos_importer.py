import numpy as np
from numba import int32, float32, char    # import the types
from numba.experimental import jitclass

class EoS:
    neutron_mass = 939.565379 # MeV

    def __init__(self, name):
        self.name = name
        self.load()
    
    def load(self):
        self.raw_data = np.loadtxt(self.name + "_" + "eos.table").T

        # the following is in compose units
        self.baryonNumberDensity = self.raw_data[1] # fm^-3
        self.pressure = self.raw_data[3] # MeV fm^-3
        self.totalEnergyPerBaryon = self.raw_data[4] # MeV
        self.totalEnergyDensity = self.raw_data[4] # MeV fm^-3
        self.restMassDensity = self.raw_data[1] * self.neutron_mass # MeV fm^-3

        # add the 0 value for better convergence
        # TODO: add check if the data is already provided like this.
        # NO: this actually doesn't seem to work, because this will generate something that is not well behaved towards the origin
        # self.baryonNumberDensity = np.insert(self.baryonNumberDensity, 0, 0.0) # fm^-3
        # self.pressure = np.insert(self.pressure, 0, 0.0) # MeV fm^-3
        # self.totalEnergyPerBaryon = np.insert(self.totalEnergyPerBaryon, 0, 0.0) # MeV
        # self.totalEnergyDensity = np.insert(self.totalEnergyDensity, 0, 0.0) # MeV fm^-3
        # self.restMassDensity = np.insert(self.restMassDensity, 0, 0.0) # MeV fm^-3

        # the following is in code units (c = G = Msun = 1). Code units are depicted with cu
        self.pressure_cu = self.pressure * 2.88e-6
        # self.totalEnergyPerBaryon_cu = self.raw_data[4]
        self.totalEnergyDensity_cu = self.totalEnergyDensity * 2.88e-6
        self.restMassDensity_cu = self.restMassDensity * 2.886e-6
        self.amtEntries = len(self.pressure)

        self.totalEnergyDensity_vs_pressure_cu = lambda x : np.interp(x, self.pressure_cu, self.totalEnergyDensity_cu)
        self.restMassDensity_vs_pressure_cu = lambda x : np.interp(x, self.pressure_cu, self.restMassDensity_cu)

        self.minPressure = min(self.pressure_cu) # TODO: think of something better on where to take the fermion radius. Maybe I could interpolate starting from rho = 0 -> p = 0?
        # self.minPressure = 0

    def get_restMassDensity(self, pressure):        
        return self.restMassDensity_vs_pressure_cu(pressure)
        
    def get_totalEnergyDensity(self, pressure):        
        return self.totalEnergyDensity_vs_pressure_cu(pressure)
    
    def get_totalInternalEnergy(self, pressure):
        totEnergy = self.get_totalEnergyDensity(pressure)
        restMassDensity = self.get_restMassDensity(pressure)
        return totEnergy / restMassDensity - 1.0

spec = [
    # ('name', char[:]),               # a simple scalar field
    ('K', float32),          # an array field
    ('gamma', float32),          # an array field
    ('minPressure', float32),          # an array field
]

# @jitclass(spec)
class EoS_polytrope:
    def __init__(self):
        # self.name = "polytropic"
        self.K = 100.0
        self.gamma = 2.0
        self.minPressure = 0.0
        pass

    def get_pressure(self, restMassDensity):
        return self.K * restMassDensity**self.gamma

    def get_restMassDensity(self, pressure):

        if(pressure < self.minPressure):
            return 0.0
        return np.power(pressure / self.K, 1.0 / self.gamma)
    
    def get_totalInternalEnergy(self, pressure):
        rho = self.get_restMassDensity(pressure)
        return self.K * rho**(self.gamma - 1.0) / (self.gamma - 1.0)

    def get_totalEnergyDensity(self, pressure):
        print("NOT IMPLEMENTED LOLOLLOL")
        return 0.0
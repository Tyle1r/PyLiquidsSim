# rocket_engine/liquid_rocket_engine.py

from .combustion_chamber import CombustionChamber
from .fuel_system import FuelSystem
from .oxidizer_system import OxidizerSystem
from .nozzle import Nozzle

class LiquidRocketEngine:
    def __init__(self, combustion_chamber, fuel_system, oxidizer_system, nozzle):
        self.combustion_chamber = combustion_chamber
        self.fuel_system = fuel_system
        self.oxidizer_system = oxidizer_system
        self.nozzle = nozzle

    def start(self):
        self.fuel_system.pressurize()
        self.oxidizer_system.pressurize()
        self.combustion_chamber.ignite()
        thrust = self.nozzle.calculate_thrust(self.combustion_chamber)
        return thrust

    def shutdown(self):
        self.combustion_chamber.extinguish()
        self.fuel_system.depressurize()
        self.oxidizer_system.depressurize()

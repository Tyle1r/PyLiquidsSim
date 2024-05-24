# rocket_engine/combustion_chamber.py

class CombustionChamber:
    def __init__(self, volume, pressure, temperature):
        self.volume = volume
        self.pressure = pressure
        self.temperature = temperature
        self.is_ignited = False

    def ignite(self):
        if not self.is_ignited:
            self.is_ignited = True
            print("Combustion chamber ignited.")

    def extinguish(self):
        if self.is_ignited:
            self.is_ignited = False
            print("Combustion chamber extinguished.")

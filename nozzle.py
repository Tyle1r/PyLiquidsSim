# rocket_engine/nozzle.py
# from .nozzle_utils import calculate_throat_area
class Nozzle:
    def __init__(self, throat_area, exit_area, efficiency):
        self.throat_area = throat_area
        self.exit_area = exit_area
        self.efficiency = efficiency

    def calculate_thrust(self, combustion_chamber):
        if combustion_chamber.is_ignited:
            thrust = (combustion_chamber.pressure * self.exit_area * self.efficiency)
            print(f"Thrust generated: {thrust} N")
            return thrust
        else:
            print("No thrust generated, combustion chamber is not ignited.")
            return 0

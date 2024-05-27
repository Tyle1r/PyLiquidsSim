# rocket_engine/fuel_system.py

class FuelSystem:
    def __init__(self, fuel_type, flow_rate, pressure):
        self.fuel_type = fuel_type
        self.flow_rate = flow_rate
        self.pressure = pressure

    def pressurize(self):
        print("Fuel system pressurized.")

    def depressurize(self):
        print("Fuel system depressurized.")

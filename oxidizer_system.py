# rocket_engine/oxidizer_system.py

class OxidizerSystem:
    def __init__(self, oxidizer_type, flow_rate, pressure):
        self.oxidizer_type = oxidizer_type
        self.flow_rate = flow_rate
        self.pressure = pressure

    def pressurize(self):
        print("Oxidizer system pressurized.")

    def depressurize(self):
        print("Oxidizer system depressurized.")

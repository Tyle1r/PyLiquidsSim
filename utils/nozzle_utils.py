# rocket_engine/nozzle_utils.py

def calculate_throat_area(exit_area, exit_mach, specific_heat_ratio=1.4):
    """
    Calculate the throat area of a nozzle given the exit area and exit Mach number.
    
    Args:
    exit_area (float): Exit area of the nozzle.
    exit_mach (float): Exit Mach number.
    specific_heat_ratio (float): Specific heat ratio (default is 1.4 for air).
    
    Returns:
    float: Throat area of the nozzle.
    """
    return exit_area / ((1 + 0.5 * (specific_heat_ratio - 1) * exit_mach ** 2) ** (1 / (specific_heat_ratio - 1)))

# Add other utility functions as needed

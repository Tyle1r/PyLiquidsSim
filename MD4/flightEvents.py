def flightEvents(t,y,alt):
    # Checks events based off of a wanted altitude 
    value = y[1] - alt
    direction = 0
    isterminal = 1
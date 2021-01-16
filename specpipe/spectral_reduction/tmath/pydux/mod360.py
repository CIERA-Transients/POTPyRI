def mod360(angle):
    import numpy as np
    angle=np.fmod(angle,360.0)
#    if angle < 0:
#        angle=angle+360.0
    return angle

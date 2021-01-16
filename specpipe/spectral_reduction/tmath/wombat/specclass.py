class Spec():
    import numpy as np
    def __init__(self,
                 wave=np.zeros(0),  #wavelength
                 flux=np.zeros(0),  # flux
                 var=np.zeros(0),  # variance, not sigma
                 obname='',  # object name
                 header=[''],  # header
            ):
        self.wave=wave
        self.flux=flux
        self.var=var
        self.obname=obname
        self.header=header

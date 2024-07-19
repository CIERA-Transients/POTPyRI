import numpy as np
import random
import matplotlib.pyplot as plt
import dustmaps.sfd

def uv_extinction(w):
    Rv=3.1
    x=1/w 
    Fax=-0.04473*(x-5.9)**2-0.009779*(x-5.9)**3
    Fbx=0.2130*(x-5.9)**2+0.1207*(x-5.9)**3
    A = []
    if x<=1/5.9:
        a=1.752-0.316*x-0.104/((x-4.67)**2+0.341)+Fax
        b=-3.090+1.825*x+1.206/((x-4.62)**2+0.263)+Fbx
    else:
        a=1.752-0.316*x-0.104/((x-4.67)**2+0.341)
        b=-3.090+1.825*x+1.206/((x-4.62)**2+0.263)
    return a+b/Rv

def opt_extinction(w):
    Rv=3.1
    y=1/w-1.82
    a=1+0.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4+0.01979*y**5-0.77530*y**6+0.32999*y**7
    b=1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-0.62251*y**5+5.30260*y**6-2.09002*y**7
    return a+b/Rv

def nir_extinction(w):
    Rv=3.1
    x=1/w
    a=0.574*x**1.61
    b=-0.527*x**1.61
    return a+b/Rv

def filter_wavelength(fil):
    if fil == 'u':
        wave = 0.36
    if fil == 'g':
        wave = 0.47
    if fil == 'r':
        wave = 0.62
    if fil == 'i':
        wave = 0.75
    if fil == 'z':
        wave = 0.89
    if fil == 'Y':
        wave = 0.99
    if fil == 'J':
        wave = 1.25
    if fil == 'H':
        wave = 1.66
    if fil == 'K' or 'Ks':
        wave = 2.19
    return wave

def calculate_mag_extinction(coords,fil):
    dustmaps.sfd.fetch()
    sfd = dustmaps.sfd.SFDQuery()
    AV = sfd(coords)*3.1/1.14
    wave = filter_wavelength(fil)
    if wave<=0.303:
        A = uv_extinction(wave)
    elif wave>0.303 and wave<=0.91:
        A = opt_extinction(wave)
    else:
        A = nir_extinction(wave)
    return A*AV
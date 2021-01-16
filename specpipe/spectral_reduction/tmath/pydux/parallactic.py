def parallactic(lat,dec,ha):
    import numpy as np
    degrad=180.0/np.pi
    


    # rewrote this.  Alex's paper gives the sine formula, but you can
    # use the tangent to get the quadrant correct without all this
    # nonsense

    # b=(sin(plat)*sin(pdec)+cos(plat)*cos(pdec)*cos(pha))
    # eta=(sin(pha)*cos(plat))/(sqrt(1-b*b))
    
    # if (dec < lat):
    #     etadeg=asin(eta)
    #     etadeg=etadeg*180.0/pi
    # else:
    #     colat=pi/2.0 - plat
    #     codec=pi/2.0 - pdec
    #     hacrit=1.0-(cos(colat)*cos(colat))/(cos(codec)*cos(codec))
    #     hacrit=sqrt(hacrit)/sin(colat)
    #     hacrit=asin(hacrit)
    #     if (abs(pha) > abs(hacrit)):
    #         etadeg=aasin(eta)*180.0/pi
    #     else: 
    #         if (pha > 0):
    #             etadeg=180-asin(eta)*180.0/pi
    #         if (pha < 0):
    #             etadeg=-180-asin(eta)*180.0/pi

    etadeg=np.arctan2(np.sin(ha),(np.cos(dec)*np.tan(lat)-np.cos(ha)*np.sin(dec)))
    etadeg=etadeg*degrad
    if (etadeg < 0):
        etadeg=etadeg+360.0
    
    return etadeg

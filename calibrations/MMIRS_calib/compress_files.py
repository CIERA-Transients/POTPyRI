from astropy.io import fits
import glob
import os

for file in glob.glob('*.gz'):

    hdulist = fits.HDUList()
    prihdu = fits.PrimaryHDU()
    imghdu = fits.ImageHDU()

    hdu = fits.open(file)

    imghdu.data = hdu[1].data
    imghdu.name = hdu[1].name

    hdulist.append(prihdu)
    hdulist.append(imghdu)

    newfilename = file.replace('.gz','')
    hdulist.writeto(newfilename, overwrite=True, output_verify='silentfix')

    os.system(f'fpack {newfilename}')


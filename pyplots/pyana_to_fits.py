import pyana
from astropy.io import fits
import sys

file_in = sys.argv[1]+'.f0'
file_out = sys.argv[1]+'.fits'

temp = pyana.fzread(file_in)
cube = temp["data"]

print cube.shape

hdu = fits.PrimaryHDU(cube)
hdulist = fits.HDUList(hdu)
hdulist.writeto(file_out)








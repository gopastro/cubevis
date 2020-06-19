from sculpt.idealpy.fits import sxpar
import math
import numpy

def xyad_cosdec(header, x, y):
    """
    Given header and x, y pixels returns ra and dec coordinates
    corresponding to that x and y. This corrects for the
    cosine(declination) term for large fields-of-view.

    @param header:  A pyfits style header object
    @param x: x pixel location (assumes RA axis)
    @param y: y pixel location (assumes Dec axis)
    @return : A tuple of (ra, dec) in decimal degrees
       """
    naxis = sxpar(header, 'NAXIS')
    for i in range(1, naxis+1):
        if "RA" in sxpar(header, 'CTYPE%d' % i):
            ra_axis = i
        if "DEC" in sxpar(header, 'CTYPE%d' % i):
            dec_axis = i            
    crval1 = sxpar(header,'CRVAL%d' % ra_axis)
    crval2 = sxpar(header,'CRVAL%d' % dec_axis)
    cdelt1 = sxpar(header,'CDELT%d' % ra_axis)
    cdelt2 = sxpar(header,'CDELT%d' % dec_axis)
    crpix1 = sxpar(header,'CRPIX%d' % ra_axis)
    crpix2 = sxpar(header,'CRPIX%d' % dec_axis)
    dec = crval2 + (y - crpix2)*cdelt2
    ra = crval1 + (x - crpix1)*cdelt1/numpy.cos(numpy.radians(dec)) 

    return (ra, dec)



def xyad(header, x, y):
    """
    Given header and x, y pixels returns world coordinates
    corresponding to that x and y. 

    @param header:  A pyfits style header object
    @param x: x pixel location (assumes RA axis)
    @param y: y pixel location (assumes Dec axis)
    @return : A tuple of world coordinates in decimal degrees
       """
    #naxis = sxpar(header, 'NAXIS')
    # for i in range(1, naxis+1):
    #     if "RA" in sxpar(header, 'CTYPE%d' % i):
    #         ra_axis = i
    #     if "DEC" in sxpar(header, 'CTYPE%d' % i):
    #         dec_axis = i            
    ra_axis = 2
    dec_axis = 3
    crval1 = sxpar(header,'CRVAL%d' % ra_axis)
    crval2 = sxpar(header,'CRVAL%d' % dec_axis)
    cdelt1 = sxpar(header,'CDELT%d' % ra_axis)
    cdelt2 = sxpar(header,'CDELT%d' % dec_axis)
    crpix1 = sxpar(header,'CRPIX%d' % ra_axis)
    crpix2 = sxpar(header,'CRPIX%d' % dec_axis)
    dec = crval2 + (y - crpix2)*cdelt2
    ra = crval1 + (x - crpix1)*cdelt1

    return (ra, dec)


from sculpt.idealpy.fits import sxpar
import math

def ad_cdec_xy(header, a, d):
    """
    Given header and ra, dec returns x and y pixel locations
    for that ra and dec. This is similar to L{adxy} except
    that it keeps track of the cosine(declination) term
    for large fields-of-view.

    @param header:  A pyfits style header object
    @param a: RA in decimal degrees
    @param d: DEC in decimal degrees
    @return : A tuple of (x, y) pixel locations
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
    x = crpix1+math.cos(math.radians(d))*(a-crval1)/cdelt1
    y = crpix2+(d-crval2)/cdelt2

    return (x, y)


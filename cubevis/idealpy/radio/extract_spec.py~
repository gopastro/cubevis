import numpy
from astropy.io import fits as pyfits
import copy
import types
import astropy.io

from sculpt.utils import SculptArgumentError
from sculpt.idealpy.fits import sxpar, sxdelpar, sxaddpar, getaxes, sxaddhist
from sculpt.idealpy.radio.smooth_image import gauss_kern

def extract_spec(hdu, x0, y0, header=None, gauss_width=2):
    """
    Given a hdu in vlm format, and a x0,y0 location (in pixel
    coordinates), obtains a gaussian-convolved spectra
    around x0, y0 with width gauss_width.
    Returns spectrum as a pyfits hdu unit as well.
    @param hdu: input pyfits style HDU (header data unit) or just the numpy
        numpy format data cube from pyfits HDU data attribute.
        If numpy format data cube is passed, then the header parameter
        should also be passed in. The cube is expected to be in vlm format.
    @type hdu: pyfits hdu type or numpy nd-array
    @param x0: x pixel location where spectrum will be centered
    @type x0: float or int
    @param y0: y pixel location where spectrum will be centered
    @type y0: float or int    
    @param header: pyfits header object if needed
    @gauss_width: the width of the gaussian kernel to use in weighting
        neighboring pixels with.
    @type gauss_width: int
    @return A HDU instance with 1d spectrum map in the data attribute of the HDU, 
        and the corresponding header in header attribute of the HDU.
    
    """
    if isinstance(hdu, astropy.io.fits.hdu.image.PrimaryHDU):
        #get data and header from the hdu
        data = hdu.data
        header = hdu.header
    elif isinstance(hdu, numpy.ndarray):
        if header is None or not isinstance(header, astropy.io.fits.header.Header):
            raise SculptArgumentError('header', "Since you passed in data that is a numpy array, set header to a pyfits header type")
        data = hdu
    else:
        raise SculptArgumentError('hdu', "can only be one of pyfits.PrimaryHDU type or numpy ndarray")
    hdr = header.copy()
    dt = data.copy()
    if gauss_width <= 0 or type(gauss_width) != types.IntType:
        raise SculptArgumentError('gauss_width', "should be positive and non-zero integer")
    gk = gauss_kern(gauss_width)
    kernx, kerny = gk.shape
    ysize, xsize, vsize = dt.shape
    if x0 < 0 or x0 >= xsize:
        raise SculptArgumentError('x0', "x0=%s is out of bounds of x-limits: (0, %d)" % (x0, xsize))
    if y0 < 0 or y0 >= ysize:
        raise SculptArgumentError('y0', "y0=%s is out of bounds of y-limits: (0, %d)" % (y0, ysize))
    #print gk.shape
    xmin = int(round(x0))-kernx/2
    xmax = xmin+kernx
    if xmin<0:
        xmin = 0
    if xmax >= xsize:
        xmax = xsize-1
    ymin = int(round(y0))-kerny/2
    ymax = ymin+kerny
    if ymin<0:
        ymin = 0
    if ymax >= ysize:
        ymax = ysize-1
    specdata = numpy.zeros(vsize, dtype='float')
    #print "Gk = " , gk
    #print xmin, xmax, ymin, ymax
    #print gk.sum().shape
    if xmax-xmin != kernx or ymax-ymin != kerny:
        specdata = dt[y0, x0, :]
    else:
        for i in range(vsize):
            specdata[i] = (gk*dt[ymin:ymax, xmin:xmax, i]).sum()/gk.sum()

    sxaddpar(hdr, "NAXIS", 1)
    sxdelpar(hdr, "NAXIS2")
    sxdelpar(hdr, "NAXIS3")
    for i in range(2,4):
        for name in ('CTYPE', 'CRVAL', 'CDELT', 'CRPIX'):
            sxdelpar(hdr, "%s%i" % (name, i))
    sxaddhist(hdr, "Extracted spectrum at (%s, %s) with gauss_width=%s" % (x0,y0,gauss_width))
    return pyfits.PrimaryHDU(specdata, header=hdr)


    

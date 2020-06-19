import numpy
from astropy.io import fits as pyfits
import copy
import types
import math
import astropy.io

from cubevis.utils import CubeVisArgumentError
from cubevis.idealpy.fits import sxpar, sxdelpar, sxaddpar, getaxes, sxaddhist, xyad
from cubevis.idealpy.radio.extract_spec import extract_spec

def extract_posvel(hdu, p1, p2, header=None, gauss_width=2):
    """
    Given a hdu in vlm format, and two positions (in pixel
    coordinates), obtains a position velocity cut along that line
    Returns p-v image as a pyfits hdu unit.
    @param hdu: input pyfits style HDU (header data unit) or just the numpy
        numpy format data cube from pyfits HDU data attribute.
        If numpy format data cube is passed, then the header parameter
        should also be passed in. The cube is expected to be in vlm format.
    @type hdu: pyfits hdu type or numpy nd-array
    @param p1: (x, y) pixel location of position 1. The x,y position can be a
        two-tuple, a list of two elements, or numpy array of two elements
    @type p1: tuple, list or array of floats or ints
    @param p2: (x, y) pixel location of position 2. The x,y position can be a
        two-tuple, a list of two elements, or numpy array of two elements
    @type p2: tuple, list or array of floats or ints
    @param header: pyfits header object if needed
    @gauss_width: the width of the gaussian kernel to use in weighting
        neighboring pixels with when extracting spectra
    @type gauss_width: int
    @return A HDU instance with 2d position-velocity image in the data
        attribute of the HDU, and the corresponding header in header
        attribute of the HDU.
    """
    if isinstance(hdu, astropy.io.fits.hdu.image.PrimaryHDU):
        #get data and header from the hdu
        data = hdu.data
        header = hdu.header
    elif isinstance(hdu, numpy.ndarray):
        if header is None or not isinstance(header, astropy.io.fits.header.Header):
            raise CubeVisArgumentError('header', "Since you passed in data that is a numpy array, set header to a pyfits header type")
        data = hdu
    else:
        raise CubeVisArgumentError('hdu', "can only be one of pyfits.PrimaryHDU type or numpy ndarray")
    hdr = header.copy()
    dt = data.copy()
    hdu = pyfits.PrimaryHDU(dt, header=hdr)
    if gauss_width <= 0 or type(gauss_width) != types.IntType:
        raise CubeVisArgumentError('gauss_width', "should be positive and non-zero integer")
    try:
        if len(p1) != 2 and len(p2) != 2:
            raise CubeVisArgumentError('positions', "Positions p1 and p2 should be 2-element tuples, lists or arrays")
    except:
        raise CubeVisArgumentError('positions', "Positions p1 and p2 should be 2-element tuples, lists or arrays")
    x1, y1 = p1
    x2, y2 = p2
    m = (y1-y2)/(x1-x2)  #slope of line
    d = math.sqrt((x1-x2)**2 + (y1-y2)**2.)
    posvel = numpy.zeros((int(d), dt.shape[2]), dtype='float')
    dist = numpy.arange(int(d))
    if x1 <= x2:
        x = x1 + numpy.sqrt(dist**2./(1+m**2.)) #2nd term is deltax
    else:
        x = x1 - numpy.sqrt(dist**2./(1+m**2.)) #2nd term is deltax
    y = y1 + m*(x-x1)
    for i in range(int(d)):
        posvel[i,:] = extract_spec(hdu, x[i], y[i],
                                   gauss_width=gauss_width).data

    xmid = (x1+x2)/2.
    ymid = (y1+y2)/2.
    r, d = xyad(header, xmid, ymid)
    sxaddpar(hdr, "CRVAL2", r)
    mid = posvel.shape[0]/2.
    sxaddpar(hdr, "CRPIX2", mid)
    sxaddpar(hdr, "NAXIS", 2)
    sxdelpar(hdr, "NAXIS3")

    for name in ('CTYPE', 'CRVAL', 'CDELT', 'CRPIX'):
            sxdelpar(hdr, "%s3" % name)
    sxaddhist(hdr, "Extracted posvel image from (%.1f, %.1f) to (%.1f, %.1f) with gauss_width=%s" % (x1, y1, x2, y2, gauss_width))
    return pyfits.PrimaryHDU(posvel, header=hdr)


    

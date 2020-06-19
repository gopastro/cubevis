import numpy
from astropy.io import fits as pyfits
import copy
import types
import math
import astropy.io

from cubevis.utils import CubeVisArgumentError
from cubevis.idealpy.fits import sxpar, sxdelpar, sxaddpar, getaxes, sxaddhist, xyad
from cubevis.idealpy.radio.extract_spec import extract_spec
from cubevis.idealpy.utils import line_rectangle_intersection

def extract_posvel_angle(hdu, p1, angle, header=None, gauss_width=2):
    """
    Given a hdu in vlm format, and one position (in pixel
    coordinates), and an angle, obtains a position velocity cut about
    that position with that angle, along that line
    Returns p-v image as a pyfits hdu unit.
    @param hdu: input pyfits style HDU (header data unit) or just the numpy
        numpy format data cube from pyfits HDU data attribute.
        If numpy format data cube is passed, then the header parameter
        should also be passed in. The cube is expected to be in vlm format.
    @type hdu: pyfits hdu type or numpy nd-array
    @param p1: (x, y) pixel location of position 1. The x,y position can be a
        two-tuple, a list of two elements, or numpy array of two elements
    @type p1: tuple, list or array of floats or ints
    @param angle: angle of PV cut with respect to X-axis in degrees
    @type angle: float
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
        if len(p1) != 2:
            raise CubeVisArgumentError('positions', "Positions p1 and p2 should be 2-element tuples, lists or arrays")
    except:
        raise CubeVisArgumentError('positions', "Positions p1 and p2 should be 2-element tuples, lists or arrays")
    a, b = p1
    #x2, y2 = p2
    m = math.tan(math.radians(angle)) #slope of line
    xmin, ymin = 0, 0
    y1, x1, v1 = hdu.data.shape
    ymax = y1-1
    xmax = x1-1
    print("xylims: ", xmin, ymin, xmax, ymax)
    # #check for x-limits first
    # xl = None
    # xh = None
    # if m > 0.0:
    #     #line sloped right of vertical
    #     yl = m*(xmin-a) + b
    # else:
    #     #line sloped left of vertical
    #     yl = m*(xmax-a) + b
    # if ymin <= yl and yl <= ymax:
    #     if m > 0.0:
    #         xl = xmin
    #     else:
    #         xl = xmax
    # if m > 0.0:
    #     yh = m*(xmax-a) + b
    # else:
    #     yh = m*(xmin-a) + b
    # print "yl, yh = ", yl, yh
    # if ymin <= yh and yh <= ymax:
    #     if m > 0.0:
    #         xh = xmax
    #     else:
    #         xh = xmin
    # if xl is None:
    #     yl = None
    #     if m > 0.0:
    #         xl = ((ymin-b)/m) + a
    #     else:
    #         xl = ((ymax-b)/m) + a
    #     if xmin <= xl and xl <= xmax:
    #         yl = ymin
    # if xh is None:
    #     yh = None
    #     if m > 0.0:
    #         xh = ((ymax-b)/m) + a
    #     else:
    #         xh = ((ymin-b)/m) + a
    #     if xmin <= xh and xh <= xmax:
    #         yh = ymax
    i1, i2 = line_rectangle_intersection(xmin, xmax, ymin, ymax, a, b, angle)
    xl,yl = i1
    xh,yh = i2
    print(xl, xh, yl, yh)
    d = math.sqrt((xl-xh)**2 + (yl-yh)**2.)
    posvel = numpy.zeros((int(d), dt.shape[2]), dtype='float')
    dist = numpy.arange(int(d))
    if xl <= xh:
        x = xl + numpy.sqrt(dist**2./(1+m**2.)) #2nd term is deltax
    else:
        x = xl - numpy.sqrt(dist**2./(1+m**2.)) #2nd term is deltax
    y = yl + m*(x-xl)
    for i in range(int(d)):
        posvel[i,:] = extract_spec(hdu, x[i], y[i],
                                   gauss_width=gauss_width).data
    r, d = xyad(header, a, b)
    sxaddpar(hdr, "CRVAL2", r)
    dmid = math.sqrt((xl-a)**2 + (yl-b)**2)
    sxaddpar(hdr, "CRPIX2", int(dmid))

    sxaddpar(hdr, "NAXIS", 2)
    sxdelpar(hdr, "NAXIS3")

    for name in ('CTYPE', 'CRVAL', 'CDELT', 'CRPIX'):
            sxdelpar(hdr, "%s3" % name)
    sxaddhist(hdr, "Extracted posvel image from (%.1f, %.1f) with angle %.2f with gauss_width=%s" % (a, b, angle, gauss_width))
    return pyfits.PrimaryHDU(posvel, header=hdr)


    

from sculpt.utils import SculptArgumentError
from sculpt.idealpy.fits import sxpar, sxaddpar, sxaddhist
import numpy
from scipy import signal
from astropy.io import fits as pyfits
import copy
import astropy.io

def cube_extract(hdu, vxy, header=None):
    """
    Given a 3-dimensionsal FITS format cube, returns a portion (sub-cube)
    of it.
    @param hdu: input pyfits style HDU (header data unit) or just the numpy
        numpy format data image from pyfits HDU data attribute.
        If numpy format data cube is passed, then the header parameter
        should also be passed in.
    @param vxy: list or numpy array with 6 parameters, first 2 are limits of
        1st axis, 2nd pair is limits of 2nd axis, and 3rd pair are limits of
        3rd axis. Using -1 as a limit implies that you can use the limiting pixel.
        For example vxy = [-1, -1, 40, 60, 40, 60] will extract the full range
        in 1st axis, but 40:60 each in 2nd and 3rd axis.
    @type vxy: list or numpy array
    @param header: pyfits header object if needed
    @return: A HDU instance with 3d output subcube in the data attribute of the HDU, 
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
    if len(vxy) != 6:
        raise SculptArgumentError('vxy', "has to be of length 6")
    hdr = header.copy()
    dt = data.copy()
    naxis1 = sxpar(hdr, 'naxis1')
    naxis2 = sxpar(hdr, 'naxis2')
    naxis3 = sxpar(hdr, 'naxis3')
    v1, v2 = vxy[0], vxy[1]
    if v1<0:
        v1 = 0
    if v2 > naxis1 or v2<0:
        v2 = naxis1
    x1, x2 = vxy[2], vxy[3]
    if x1 < 0:
        x1 = 0
    if x2 > naxis2 or x2<0:
        x2 = naxis2
    y1, y2 = vxy[4], vxy[5]
    if y1 < 0:
        y1 = 0
    if y2 > naxis3 or y2<0:
        y2 = naxis3
    dt = dt[y1:y2, x1:x2, v1:v2]
    ny, nx, nv = dt.shape
    cards = hdr.cards
    card = cards['NAXIS1']
    card.value = nv
    card = cards['NAXIS2']
    card.value = nx
    card = cards['NAXIS3']
    card.value = ny

    #calculate velocity axis
    crval1 = sxpar(hdr, 'crval1')
    cdelt1 = sxpar(hdr, 'cdelt1')
    crpix1 = sxpar(hdr, 'crpix1')
    v = crval1+(v1-crpix1)*cdelt1
    card = cards['CRVAL1']
    card.value = v

    #x-axis
    card = cards['CRPIX2']
    card.value = card.value-x1

    #y-axis
    card = cards['CRPIX3']
    card.value = card.value-y1
    return pyfits.PrimaryHDU(dt, header=hdr)

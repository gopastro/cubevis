from sculpt.idealpy.fits import sxpar, sxdelpar, sxaddpar, getaxes, sxaddhist
from sculpt.utils import SculptArgumentError

import numpy
from astropy.io import fits as pyfits
import types
import astropy.io
#import copy

def mkrmsimage (hdu, window, header=None, 
                chan=False, dontblank=False,
                kms=True,
                moment=0):
    """
    Create rms image from spectral cube. The input cube is expected
    to be in the VLM format

    @param hdu: input pyfits style HDU (header data unit) or just the numpy
        numpy format data cube from pyfits HDU data attribute.
        This cube is expected to be in vlm format. If numpy format
        data cube is passed, then the header parameter should also be
        passed in.
    @type hdu: pyfits HDU or numpy array
    @param header: pyfits style header object that corresponds to the
        hdu data variable. This header parameter is only considered if
        the hdu parameter is a numpy array type rather than a pyfits HDU
        type
    @type header: pyfits header object
    @param window: list of tuples (pairs) of channels/velocity to exclude
        from polynomial calculation in determining rms. If chan is set to
        True, these should be integers. If chan is False, this is expected to be
        velocity expressed in kms (if the input parameter kms is True). If
        chan is False, and kms is False, v1 is treated with the same
        units as the FITS header describes.
    @param chan: If True, window variables are treated as channels. If False,
        they are treated as velocity
    @type chan: Boolean
    @param dontblank: By default, mkrmsimage will use the BLANK fits
        header value and automatically zero blanked values in order
        to compute rms. If dontblank keyword is set, then
        such blanking is not performed.
    @type dontblank: Boolean
    @param kms: If True velocity units is in kms.
    @type kms: Boolean
    @return: A HDU instance with 2d output map in the data attribute of the HDU, 
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
    #check for list or tuple type of window
    if type(window) not in (types.ListType, types.TupleType):
        raise SculptArgumentError('window', 'has to be a List Type or Tuple Type')
    # calculate the x-axis (velocity) 
    crpix1 = sxpar(header,"CRPIX1")
    crval1 = sxpar(header,"CRVAL1")
    cdelt1 = sxpar(header,"CDELT1")
    if kms:
        #convert velocity to km/s
        crval1 = crval1/1000.
        cdelt1 = cdelt1/1000.
    ctype1 = sxpar(header,"CTYPE1")
    crpix2 = sxpar(header,"CRPIX2")
    crval2 = sxpar(header,"CRVAL2")
    cdelt2 = sxpar(header,"CDELT2")
    ctype2 = sxpar(header,"CTYPE2")
    crpix3 = sxpar(header,"CRPIX3")
    crval3 = sxpar(header,"CRVAL3")
    cdelt3 = sxpar(header,"CDELT3")
    ctype3 = sxpar(header,"CTYPE3")
    nv = sxpar(header,"NAXIS1")
    nx = sxpar(header,"NAXIS2")
    ny = sxpar(header,"NAXIS3")
    blank = sxpar(header,"BLANK")

    cube = data.copy()
    if not dontblank and blank is not None:
        ind = numpy.where(cube == blank)
        cube[ind] = 0.0

    #get velocity axis
    velax = getaxes(header, 1)
    indices = numpy.zeros(velax.shape, dtype=bool)
    for v1, v2 in window:
        if chan:
            indices[v1:v2] = True
        else:
            ind = numpy.logical_and(velax>= v1, velax<= v2)
            indices = numpy.logical_or(indices, ind)
    indices = numpy.logical_not(indices) # for indices to include in std calc

    rms = cube[:,:,indices].std(axis=2)
    hnew = header.copy()
    
    sxaddpar(hnew, "CRVAL1", crval2, comment="DEGREES")
    sxaddpar(hnew, "CRPIX1", crpix2)
    sxaddpar(hnew, "CDELT1", cdelt2, comment="DEGREES")
    sxaddpar(hnew, "CTYPE1", ctype2)
    sxaddpar(hnew, "CRVAL2", crval3, comment="DEGREES")
    sxaddpar(hnew, "CRPIX2", crpix3)
    sxaddpar(hnew, "CDELT2", cdelt3, comment="DEGREES")
    sxaddpar(hnew, "CTYPE2", ctype3)
    sxaddpar(hnew, "NAXIS", 2)
    sxaddpar(hnew, "NAXIS1", nx)
    sxaddpar(hnew, "NAXIS2", ny)
    sxaddpar(hnew, "NAXIS3", 1)
    sxaddpar(hnew, "NAXIS4", 1)
    if chan:
        vorc = 'CHANNEL'
    else:
        vorc = 'VELOCITY'
    sxaddhist(hnew, "WINDOW : %s; Window %s LIMITS" % (repr(window), vorc))
    #sxaddpar(hnew, "BUNIT", units, "Units")
    sxdelpar(hnew, "CRVAL3")
    sxdelpar(hnew, "CRPIX3")
    sxdelpar(hnew, "CDELT3")
    sxdelpar(hnew, "CTYPE3")
    sxdelpar(hnew, "NAXIS3")
    sxdelpar(hnew, "NAXIS4")
    hdu = pyfits.PrimaryHDU(rms, header=hnew)
    return hdu

from sculpt.utils import SculptArgumentError
import numpy
from scipy import signal
from astropy.io import fits as pyfits
import copy
import astropy.io

def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = numpy.mgrid[-size:size+1, -sizey:sizey+1]
    g = numpy.exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None, mode='same') :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im,g, mode=mode)
    return(improc)

def smooth_image(hdu, smooth=2, header=None, mode='same'):
    """Takes a 2-dimensional FITS format image, and returns a smoothed version
    of it.

    @param hdu: input pyfits style HDU (header data unit) or just the numpy
        numpy format data image from pyfits HDU data attribute.
        If numpy format data cube is passed, then the header parameter
        should also be passed in.
    @param smooth: integer smoothing factor
    @param header: pyfits header object if needed
    @return A HDU instance with 2d output map in the data attribute of the HDU, 
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
    hdr = copy.copy(header)
    ny, nx = data.shape
    dt = copy.copy(data)
    dt = blur_image(dt, smooth, mode=mode)
    shape = dt.shape
    nyy, nxx = shape
    #cards = hdr.ascardlist()
    cards = hdr.cards
    card = cards['NAXIS1']
    card.value = shape[1]
    card = cards['CDELT1']
    card.value = card.value*float(nx)/float(nxx)
    card = cards['NAXIS2']
    card.value = shape[0]
    card = cards['CDELT2']
    card.value = card.value*float(ny)/float(nyy)
    
    return pyfits.PrimaryHDU(dt, header=hdr)

import types
from cubevis.idealpy.fits import sxpar
import numpy

def getax(header, ax):
    """
    helper function for getaxes
    """
    crval = sxpar(header, 'CRVAL%d' % ax)
    crpix = sxpar(header, 'CRPIX%d' % ax)
    cdelt = sxpar(header, 'CDELT%d' % ax)
    naxis = sxpar(header, 'NAXIS%d' % ax)
    return crval + (numpy.arange(naxis) - (crpix-1)) * cdelt

def getaxes(header, axis='all', kms=True):
    """
    L{getaxes} returns the axis vector(s) listed in keyword
    axis.

    @param header: a pyfits header instance. This is derived from
        pyfits.getheader or from the hdu instance header attribute
        of a pyfits HDU
    @param axis: can be string 'all', in which case it returns all axes
        available in the header. Returned as a dictionary with
        keys equaling axis number and values equaling corresponding
        axis vector. axis can also be set to None, which is the same
        as 'all'. axis can also be set to an integer axis number or
        a list of integer axis numbers.
    @type axis: Boolean
    @param kms: If set to True (default), if it finds velocity axis
        (it finds this by looking for CTYPE starting with velo), it
        will automatically divide by 1000. to convert velocity from
        m/s to km/s.
    @type kms: Boolean
    @return: returns axis. If only a single axis is requested (based on
        the type of axis input), a single numpy array of axis is returned.
        Otherwise a dictionary with axis number as keys and numpy array of
        axis as dictionary values is returned.
      """
    allaxes = False
    if isinstance(axis, str):
        if axis.lower() == 'all':
            allaxes = True
        else:
            print("axis can be 'all' or an integer or list of integers.")
            return None
    if axis is None:
        allaxes = True
    if not allaxes:
        axes = []
        if isinstance(axis, list) or isinstance(axis, tuple):
            for ax in axis:
                axes.append(ax)
        elif isinstance(axis, int):
            axes.append(axis)
        else:
            print("axis can be 'all' or an integer or list of integers.")
            return None
    if allaxes:
        naxis = sxpar(header, 'naxis')
        axes = range(1, naxis+1)
    axdic = {}
    for ax in axes:
        axdic[ax] = getax(header, ax)
        if kms:
            if sxpar(header, 'CTYPE%d' % ax).lower().startswith('velo'):
                #velocity axis
                axdic[ax] = axdic[ax]/1000.
    if len(axes) == 1:
        return axdic[axes[0]]
    else:
        return axdic


        

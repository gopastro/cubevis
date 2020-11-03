import numpy

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return numpy.isnan(y), lambda z: z.nonzero()[0]

def interpolate_spectrum_with_nans(spec):
    nans, x= nan_helper(spec)
    spec[nans]= numpy.interp(x(nans), x(~nans), spec[~nans])
    return spec


def interpolate_cube_with_nans(hdu):
    """
    Helper to take all spectra of a cube with nans in the spectra
    and interpolate those nans with adjacent channels.
    hdu is Header data unit.
    Data is expected in vxy format.
    The data is replaced in place
    """
    data = hdu.data
    nx, ny, _ = data.shape
    for x in range(nx):
        for y in range(ny):
            data[x, y, :] = interpolate_spectrum_with_nans(data[x, y, :])
    return hdu


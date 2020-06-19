from astropy.io import fits as pyfits
from optparse import OptionParser
import os
import sys
import numpy
import copy
from sculpt.idealpy.utils.congrid import congrid
from sculpt.utils import SculptArgumentError

def transpose_data(data, origin='vxy'):
    """always transposes vxy to xyv
    or xyv to vxy. """
    if origin == 'vxy':
        #if data is vxy, data shape is (because
        #of numpy:
        # data.shape = yxv
        # need to cast it into vyx
        data1 = numpy.swapaxes(data, 0, 2)
        #now it is vxy
        data1 = numpy.swapaxes(data1, 1, 2)
        #now it is vyx
    elif origin == 'vyx':
        # if data is vyx, data shape is (because of numpy):
        # data.shape = xyv
        # need to cast into vyx
        data1 = numpy.swapaxes(data, 0, 2)
        # now it is vyx - so we are good!
    elif origin == 'xyv':
        #input is xyv
        #data.shape = vyx
        #need to cast into yxv
        data1 = numpy.swapaxes(data, 0, 2)
        #now it is xyv
        data1 = numpy.swapaxes(data1, 0, 1)
        #now it is yxv
    elif origin == 'yxv':
        # input is yxv
        # data.shape = vxy
        # need to cast into yxv
        data1 = numpy.swapaxes(data, 0, 2)
        # now it is yxv - so we are good!
    return data1

def copycard(c1, c2):
    """Copy contents of c2 to c1"""
    c1.value = c2.value
    c1.comment = c2.comment

def transpose_header(header, origin='vxy'):
    """tranpose axis1 and axis3 information. axis1->axis3,
    axis2->axis1, axis3->axis2."""
    cards = header.cards
    keepattr = {}
    for attr in ('NAXIS', 'CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CROTA'):
        keepattr[attr] = {}
        for i in range(1, 4):
            keepattr[attr][i] = copy.copy(cards['%s%d' % (attr,i)])
    if origin == 'vxy':
        for attr in ('NAXIS', 'CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CROTA'):
            #1->3
            copycard(cards['%s3' % attr], keepattr[attr][1])
            #2->1
            copycard(cards['%s1' % attr], keepattr[attr][2])
            #3->2
            copycard(cards['%s2' % attr], keepattr[attr][3])
        dest = 'xyv'
        header.add_history('Tranposing from %s format to %s format' % (origin, dest))
    elif origin == 'vyx':
        for attr in ('NAXIS', 'CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CROTA'):
            #1->3
            copycard(cards['%s3' % attr], keepattr[attr][1])
            #3->1
            copycard(cards['%s1' % attr], keepattr[attr][3])
            dest = 'xyv'
            header.add_history('Tranposing from %s format to %s format' % (origin, dest))
    elif origin == 'xyv':
        for attr in ('NAXIS', 'CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CROTA'):
            #1->2
            copycard(cards['%s2' % attr], keepattr[attr][1])
            #2->3
            copycard(cards['%s3' % attr], keepattr[attr][2])
            #3->1
            copycard(cards['%s1' % attr], keepattr[attr][3])
        dest = 'vxy'
        header.add_history('Tranposing from %s format to %s format' % (origin, dest))
    elif origin == 'yxv':
        for attr in ('NAXIS', 'CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CROTA'):
            #3->1
            copycard(cards['%s1' % attr], keepattr[attr][3])
            #1->3
            copycard(cards['%s3' % attr], keepattr[attr][1])
        dest = 'vxy'
        header.add_history('Tranposing from %s format to %s format' % (origin, dest))
    return header

def smooth_spec(f, y, window_len=5, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett',
        'blackman'
        flat window will produce a moving average smoothing.
        
    output:
        the smoothed signal

        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman,
        numpy.convolve
        scipy.signal.lfilter
        """
    xmin, xmax  = f.min(), f.max()
    data = copy.deepcopy(y)
    if window_len<3:
        return data
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is not one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    if data.size < window_len:
        raise ValueError, "data vector needs to be bigger than window size."
    s=numpy.r_[2*data[0]-data[window_len:1:-1],data,2*data[-1]-data[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w = numpy.ones(window_len,'d')
    else:
        w = eval('numpy.'+window+'(window_len)')
    y = numpy.convolve(w/w.sum(),s,mode='same')
    data = y[window_len-1:-window_len+1]
    return data

def smooth(data, header, smooth_factor, origin='vxy'):
    if origin[0] == 'v':
        shape = copy.copy(data.shape)
        vdim = int(round(shape[0]/float(smooth_factor)))
        new_shape = (vdim, shape[1], shape[2])
        data = congrid(data, new_shape)
        cards = header.cards
        card = cards['NAXIS3']
        card.value = vdim
        card = cards['CDELT3']
        card.value = card.value*smooth_factor
        header.add_history('Smoothing velocity axis by %d' % smooth_factor)
        return data, header
    else:
        shape = copy.copy(data.shape)
        vdim = int(round(shape[0]/float(smooth_factor)))
        new_shape = (shape[1], shape[2], vdim)
        data = congrid(data, new_shape)
        cards = header.cards
        card = cards['NAXIS1']
        card.value = vdim
        card = cards['CDELT1']
        card.value = card.value*smooth_factor
        header.add_history('Smoothing velocity axis by %d' % smooth_factor)
        return data, header        
    
def transpose_cube(hdu, origin='xyv',
                   smooth=None):
    """
    Given a HDU this function transposes from origin (xyv, yxv, vxy or vyx)
    and makes a new hdu into  vxy or xyv cube
    If smooth is an integer the velocity axis is also smoothed 
    by the integral value
    """
    header = hdu.header
    data = hdu.data
    if origin not in ('xyv', 'yxv', 'vxy', 'vyx'):
        raise SculptArgumentError('transpose_cube', 'Origin should be one of xyv, yxv, vxy or vxy')
    if origin[2] == 'v':
        if header.get('CTYPE3') not in ('VELO-LSR', 'VELOCITY'):
            raise SculptArgumentError('header', 'Input Cube does not have Velocity in 3rd axis')
    else:
        if header.get('CTYPE1') not in ('VELO-LSR', 'VELOCITY'):
            raise SculptArgumentError('header', 'Input Cube does not have Velocity in 1st axis')
    data1 = transpose_data(data, origin=origin)
    header1 = transpose_header(header, origin=origin)
    if smooth is not None:
        data1, header1 = smooth(data1, header1, smooth,
                                origin=origin)
    return pyfits.PrimaryHDU(data=data1, header=header1)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read input from FILE", metavar="FILE")
    parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="don't print status messages to stdout")
    parser.add_option("-o", "--output", dest="outfilename",
                      help="output to OUTFILE", metavar="OUTFILE")
    parser.add_option("-v", "--vxy", dest="xyv",
                      action="store_false", default=True,
                      help="convert output to vxy format")
    parser.add_option("-x", "--xyv", dest="xyv",
                      action="store_true", default=True,
                      help="convert output to xyv format")
    parser.add_option("-s", "--smooth", dest="smooth",
                      type="int", default=8)
    (options, args) = parser.parse_args()
    if options.xyv:
        out_ext = '_xyv'
    else:
        out_ext = '_vxy'
    if options.filename is None:
        if not args:
            parser.error("Need to give filename")
        else:
            filename = args[0]
    else:
        filename = options.filename
    if options.outfilename is None:
        base, ext = os.path.splitext(filename)
        outfilename = base+out_ext+ext
    print filename, outfilename

    if os.path.exists(filename):
        fp = pyfits.open(filename)
        header = fp[0].header
        if options.xyv and header.get('CTYPE1') != 'VELO-LSR':
            print "Input header is not of vxy type"
            sys.exit(-1)
        if options.xyv:
            data = fp[0].data
            #input data is in vxy
            data1 = transpose_data(data, origin='vxy')
            header1 = transpose_header(header, origin='vxy')
            print header1
            print "Smoothing..."
            if options.smooth:
                data, header = smooth(data1, header1, options.smooth,
                                      origin='vxy')
            fp[0].header = header
            fp[0].data = data
            print header
            fp.writeto(outfilename, clobber=True)

            
            

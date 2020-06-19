"""
Various FITS file utilities. All these utilities require
the L{pyfits} module from STScI.

Available functions are:

  - L{sxpar}: extract item(s) from FITS header
  - L{sxaddpar}: add items (including history and comments) to FITS header
  - L{getaxes}: From a FITS header get the vector for a given axis
"""

#__docformat__ = 'epytext en'
from .header_utils import sxpar, sxdelpar, sxaddpar, sxaddhist
from .ad_cdec_xy import ad_cdec_xy
from .xyad import xyad
#from sxaddpar import *
#from sxaddhist import *
#from sxdelpar import sxdelpar
from .getaxes import getaxes
from .fitsio import fitsdir

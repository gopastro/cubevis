"""
Various FITS IO utilities. These all use the nice
L{pyfits} library
"""

from astropy.io import fits as pyfits
#import pyfits
import glob
import types
from sculpt.idealpy.fits import sxpar

def _get_default_keywords(keywords, default):
    if keywords is None:
        return default
    elif type(keywords) in (types.ListType, types.TupleType):
        return keywords
    elif type(keywords) == types.StringType:
        return keywords.split(',')
    else:
        #keywords is not a string type nor a list type
        return default

def fitsdir(directory=None, keywords=None,
            nosize=False,
            alt1_keywords=None, alt2_keywords=None,
            alt3_keywords=None, ext=0):
    """
    The values of either user-specified or default FITS keywords are 
    displayed in either the primary header and/or the first extension header.
    Unless nosize option is set to True, the data size is also displayed.
    The default keywords are as follows (with keywords in 2nd row used if
    those in the first row not found, and the 3rd row if neither the keywords
    in the first or second rows found:)

     DATE-OBS     TELESCOP   OBJECT    EXPTIME       
     TDATEOBS     TELNAME    TARGNAME  INTEG        #First Alternative
     DATE         OBSERVAT             EXPOSURE     #Second Alternative
                  INSTRUME             EXPTIM       #Third Alternative

     fitsdir will also recognize gzip compressed files (must have a .gz 
     or FTZ extension).

     Example Usage:

        >>> fitsdir

     @param directory: Scalar string giving file name, disk or directory to be 
        searched. Wildcard file names are allowed. If left blank, the current
        working directory is searched for FITS files. Examples of 
        valid names include 'iraf/*.fits' (Unix).
     @type directory: string
     @param keywords: FITS keywords to display, as either a list or tuple
        of strings or as a comma delimited scalar string,
        e.g.'testname,dewar,filter'. If not supplied, then the default
        keywords are 'DATE-OBS', 'TELESCOP','OBJECT','EXPTIME'
     @type keywords: comma-delimited string or list or tuple of strings
     @param alt1_keywords: A list (either a list or tuple of strings)
        or a comma delimited strings of alternative keywords to use if the
        default keywords cannot be found.   By default, 'TDATEOBS', is the 
        alternative to DATE-OBS, 'TELNAME' for 'TELESCOP','TARGNAME'
        for 'OBJECT', and 'INTEG' for EXPTIME
     @type alt1_keywords: comma-delimited string or list or tuple of strings        
     @param alt2_keywords: A list (either a list or tuple of strings)
        or a comma delimited strings of alternative keywords to use if
        neither keywords nor alt1_keywords can be found.    
     @type alt2_keywords: comma-delimited string or list or tuple of strings
     @param alt3_keywords: A list (either a list or tuple of strings) or
        a comma delimited strings of alternative keywords to use if
        neither keywords nor alt1_keywords nor alt2_keywords can be found.
     @type alt3_keywords: comma-delimited string or list or tuple of strings        
     @param nosize: if set to True, then information about the image
        size is not displayed
     @type nosize: Boolean
     @param ext: extension number to get. Default 0
     @type ext: Integer
     """
    if directory is None:
        directory = '*.fits*'
    dirlist = glob.glob(directory)
    keywords = _get_default_keywords(keywords,
                                     ['date-obs','telescop','object','exptime'])
    alt1_keywords = _get_default_keywords(alt1_keywords,
                                          ['tdateobs','telname','targname','integ'])
    alt2_keywords = _get_default_keywords(alt2_keywords,
                                          ['date','observat','','exposure'])
    alt3_keywords = _get_default_keywords(alt3_keywords,
                                          ['','instrume','','exptim' ])
    if nosize:
        print "Filename\t " + '\t '.join(keywords)
    else:
        print "Filename\t Size\t " + '\t '.join(keywords)
    for fname in dirlist:
        header = pyfits.getheader(fname, ext=ext)
        sz = sxpar(header, 'naxis*')
        sizeint = []
        for i in range(1, sz['NAXIS']+1):
            sizeint.append("%d" % sz['NAXIS%d' % i])
        size = ' x '.join(sizeint)
        keyvalues = []
        for i, kw in enumerate(keywords):
            if header.has_key(kw):
                keyvalues.append("%s" % sxpar(header, kw))
            elif header.has_key(alt1_keywords[i]):
                keyvalues.append("%s" % sxpar(header, alt1_keywords[i]))
            elif header.has_key(alt2_keywords[i]):
                keyvalues.append("%s" % sxpar(header, alt2_keywords[i]))
            elif header.has_key(alt3_keywords[i]):
                keyvalues.append("%s" % sxpar(header, alt3_keywords[i]))
            else:
                keyvalues.append('None')
        if nosize:
            print "%s\t " % fname + '\t '.join(keyvalues)
        else:
            print "%s\t %s\t " % (fname, size) + '\t '.join(keyvalues) 
                                
    

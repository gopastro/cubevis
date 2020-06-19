from astropy.io import fits as pyfits
import types
import astropy.io
from cubevis.idealpy.utils.outvar import OutVar
from cubevis.idealpy.utils import IdealpyError


def sxaddhist(header, history, before=None, after=None):
    """
    Add history item(s) to a pyfits header object.

    @param header: FITS header object returned by pyfits. Returns exception
        if the hdr object is not a pyfits header.
    @param history: history string or list of strings.
    @param before: Optional input. The new name or key is inserted before
        the 'before' element. before can be a string with the key name
        or integer (index). If key in 'before' is not found keyError
        is raised.
    @param after:  Optional input. The new name or key is inserted after
        the 'after' element. after can be a string with the key name or
        integer (index). If key in 'after' is not  found keyError is raised.
    @return: Updated pyfits header
    """
    if not isinstance(header, astropy.io.fits.header.Header):
        raise IdealpyError("sxaddhist", "Input header is not a pyfits header type")

    if isinstance(history, tuple) or isinstance(history, list):
        for hist in history:
            header.add_history(hist, before=before, after=after)
    else:
        header.add_history(history, before=before, after=after)
    return header


def sxdelpar(header, name, debug=False):
    """
    Deletes a card from the header with a given name

    @param header: the pyfits header instance
    @type header: pyfits header
    @param name: the key name of the header item
    @type name: string
    @return: None. Prints error message if it does not
        find the name in header
        """
    if not isinstance(header, astropy.io.fits.header.Header):
        raise IdealpyError('sxdelpar', "Input header is not a pyfits header type")
    #if header.has_key(name):
    if name in header:
        #lst = header.ascardlist()
        #lst.pop(lst.index_of(name))
        del header[name]
    else:
        if debug:
            print("Header does not contain the key %s" % name)

def sxaddpar(header, name, value, comment=None,
             before=None, after=None, savecomment=True):
    """
    Add or modify a parameter in a pyfits header object.

    @param header: FITS header object returned by pyfits. Returns
        exception if the hdr object is not a pyfits header.
    @param name:  Name of parameter. If name is already in the header
        the value and possibly comment fields are modified.  Otherwise a new 
        record is added to the header.  If name is equal to 'COMMENT'
        or 'HISTORY' or a blank string then the value will be added to 
        the record without replacement.  For these cases, the comment 
        parameter is ignored.
    @param value: Value for parameter.  The value expression must be
        of the correct type, e.g. integer, floating or string.
        Python values of True or False are converted to string values
        of 'T' or 'F'
    @param comment: Optional input. String field. The '/' is added by
        this routine. Added starting in position 31. If not supplied, or
        set equal to '', or when savecomment is set to True (default)
        then the previous comment field is retained (when found) 
    @param before: Optional input. The new name or key is inserted before
        the 'before' element. before can be a string with the key name
        or integer (index). If key in 'before' is not  found keyError is raised.
    @param after:  Optional input. The new name or key is inserted after
        the 'after' element. after can be a string with the key name or
        integer (index). If key in 'after' is not  found keyError is raised.
    @param savecomment: Optional input. Default True. If True, uses existing
        comment for name (when name is being updated), when comment field is
        not given in input. 
    @return: Updated pyfits header
       """
    if not isinstance(header, astropy.io.fits.header.Header):
        raise IdealpyError("sxaddpar", "Input header is not a pyfits header type")

    if name.lower() == 'history':
        header.add_history(value, before=before, after=after)
        return header
    if name.lower() == 'comment':
        header.add_comment(value, before=before, after=after)
        return header
    oldcard = None
    #if header.has_key(name):
    if name in header:
        #header has item already
        #oldcard = header.ascardlist()[name]
        oldcard = header.cards[name]
    if not savecomment and oldcard:
        #make sure to reset existing comment
        if not comment:
            comment = ''
    header.set(name, value, comment=comment, before=before,
               after=after)
    return header

def sxpar(header, name, count=None, comment=None,
          NoContinue=False, silent=False):
    """
    Obtain the value of a parameter in a FITS header.
    The input hdr is a pyfits header and name is a string
    for the FITS keyword.

    Example Usage:
    
        >>> header = pyfits.getheader('hardmap_12co_32.fits')
        >>> sxpar(header, 'naxis1')
        768
        >>> sxpar(header, 'naxis*')
        {'NAXIS': 3, 'NAXIS1': 768, 'NAXIS2': 225, 'NAXIS3': 180}

    @param header: FITS header object returned by pyfits. Returns exception if the
        hdr object is not a pyfits header.
    @param name: String name of the parameter to return. If name is of the form
        'keyword*' then a python dictionary is returned containing keys
        as keywordN and values of the dictionary containing the corresponding
        value of the keyword.
    @param count: Optional output. count is an L{idealpy} OutVar object
        instance and will return the number of parameters found by sxpar.
        The outvar attribute of the variable instance will contain the count.
    @type count: L{idealpy} L{OutVar} object   
    @param comment: Optional Output. Dictionary of comments associated with the
        returned values.
    @type comment: L{idealpy} L{OutVar} object   
    @param silent:  default False. If set to True, produces more verbose
        output
    @type silent: Boolean
    @return: the value of parameter in header. If parameter is double precision,
        floating, long or string, the result is of that type. If the parameter
        is logical, True or False is returned. If name was of form 'keyword*'
        then a dictionary of values are returned.

     """
    if not isinstance(header, astropy.io.fits.header.Header):
        raise IdealpyError("sxpar", "Input header is not a pyfits header type")
    #Check for right output data types
    if count is not None:
        if type(count) != OutVar:
            print("count variable should object of type OutVar")
            return None
        else:
            count.clear()            
    if comment is not None:
        if type(comment) != OutVar:
            print("comment variable should object of type OutVar")
            return None
        else:
            comment.clear()
            comment.outvar = {}
            
    name = name.strip().upper()  #Copy name, make upper case

    keys = [k[0] for k in header.items()]
    # Determine if name is of form 'keyword*'.  If so, then strip off
    # the '*', and set the vector flag.  One must consider the
    # possibility that name is an empty string.
    vector = False
    if name.count('*') == 1:
        #wildcard option
        vector = True
        name = name[:-1]
        retval = {}
        ct = 0
        for key in keys:
            if key.startswith(name):
                ct += 1
                #card = header.ascardlist()[key]
                card = header.cards[key]
                retval[key] = card.value
                if comment:
                    comment.outvar[key] = card.comment
    else:
        if name not in keys:
            print("Keyword %s not found in header" % name)
            return None
        else:
            ct = 1
            #card = header.ascardlist()[name]
            card = header.cards[name]
            retval = card.value
            if comment:
                comment.outvar[name] = card.comment
    if count:
        count.outvar = ct
    return retval
        

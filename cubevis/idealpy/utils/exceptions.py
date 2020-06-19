"""
Set of idealpy exceptions that can be raised
"""


class Error(Exception):
    """
    Base class for all exceptions of sculpt
    """
    pass

class IdealpyError(Error):
    """
    Exceptions raised when the argument to an idealpy
    function call is of the wrong type
    """
    def __init__(self, argname, reason):
        """
        @param argname: The argument name that triggered the exception.
        @type argname: String
        @param reason: An explanatory text that details the error
        @type reason: String
        """
        self.argname = argname
        self.reason = reason
        self.message = self.reason
        self.args = (self.reason,)
        
    def __str__(self):
        return "%s : %s" % (self.argname, self.reason)

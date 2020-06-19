class OutVar(object):
    """
    A simple object that provides a mutable type for
    idealpy to act on.

    L{OutVar} provides a simple way to emulate IDL's way of passing outputs
    through keyword variables.

    Instantiate an L{OutVar} object by something like:

         >>> ovar = OutVar()

    Then you can pass ovar to idealpy's functions that take keywords.
    And upon return ovar.outvar contains the required response.
    """
    def __init__(self):
        self.clear()

    def __repr__(self):
        return "%s" % self.outvar

    def clear(self):
        self.outvar = None

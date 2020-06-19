from aplpy import FITSFigure

class CubeVisFITSFigure(FITSFigure):
    """
    For now a very thin wrapper around APLPy's FITSFigure
    We may want to extend this for our future use in 
    various ways
    """
    def __init__(self, data, hdu=0, figure=None,
                 subplot=(1,1,1), downsample=False,
                 north=False, convention=None,
                 dimensions=[0, 1], slices=[],
                 auto_refresh=None, **kwargs):
        FITSFigure.__init__(self, data, hdu, figure=figure,
                            subplot=subplot, downsample=downsample,
                            north=north, convention=convention,
                            dimensions=dimensions, slices=slices,
                            auto_refresh=auto_refresh, **kwargs)
        

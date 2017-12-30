class PSFPiece :
    """Declare a minimum interface for pieces of PiecewisePSF."""

    def __init__(self) :
        """
        ASserts a minimum interface to self.

        Args: None

        Returns: None
        """

        assert(callable(self))

        assert(hasattr(self, 'integrate'))
        assert(callable(self, 'integrate'))

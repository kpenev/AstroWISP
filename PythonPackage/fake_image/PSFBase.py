class PSFBase :
    """The base class for all supported PSFs."""

    def __init__(self) :
        """
        Asserts a minimum interface to self.

        Args: None

        Returns: None
        """

        assert(callable(self))

        assert(hasattr(self, 'get_left_range'))
        assert(callable(self.span_left))

        assert(hasattr(self, 'get_right_range'))
        assert(callable(self.span_right))

        assert(hasattr(self, 'get_up_range'))
        assert(callable(self.span_up))

        assert(hasattr(self, 'get_down_range'))
        assert(callable(self.span_down))

        assert(hasattr(self, 'integrate'))
        assert(callable(self.integrate))

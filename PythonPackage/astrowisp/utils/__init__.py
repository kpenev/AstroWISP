"""Place to collect miscellaneous useful functions."""

def flux_from_magnitude(magnitude, magnitude_1adu):
    """Return the flux corresponding to the given magnitude."""

    return 10.0**((magnitude_1adu - magnitude) / 2.5)

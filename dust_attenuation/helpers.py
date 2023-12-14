import warnings
import numpy as np
from astropy.utils.exceptions import AstropyUserWarning


def _test_valid_x_range(x, x_range, outname):
    """
    Test if any of the x values are outside of the valid range

    Parameters
    ----------
    x : float array
       wavelength in microns

    x_range: 2 floats
       allowed min/max of x

    outname: str
       name of curve for error message
    """
    if np.logical_or(np.any(x < x_range[0]), np.any(x > x_range[1])):
        raise ValueError(
            "Input x outside of range defined for "
            + outname
            + " ["
            + str(x_range[0])
            + " <= x <= "
            + str(x_range[1])
            + ", x has units micron]"
        )


def _positive_klambda(klam):
    """
    Check that k-lambda arrays have no negative values

    Parameters
    ----------
    klam: float array
        k-lambda array, which should only have positive values

    Returns
    -------
    kout: float array
        array clipped to zero if necessary
    """
    if not np.all(klam >= 0):
        warnings.warn(
            "k-lambda has negative values, setting them to zero.", AstropyUserWarning
        )
        return np.maximum(klam, 0.0)
    else:
        return klam

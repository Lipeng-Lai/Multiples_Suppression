import numpy as np

def ricker(f, dt):
    """
    Ricker wavelet generator.

    Parameters
    ----------
    f : float
        Dominant frequency of the wavelet (Hz).
    dt : float
        Time interval (sampling rate).

    Returns
    -------
    w : ndarray
        Ricker wavelet.
    """
    # Define the length of the wavelet in seconds
    tmax = 1.5 / f
    t = np.arange(-tmax, tmax + dt, dt)  # Time vector centered around zero

    # Ricker wavelet formula
    pi2 = np.pi ** 2
    w = (1 - 2 * pi2 * (f ** 2) * (t ** 2)) * np.exp(-pi2 * (f ** 2) * (t ** 2))

    # Normalize the wavelet to have a peak value of 1
    w = w / np.max(np.abs(w))

    return w

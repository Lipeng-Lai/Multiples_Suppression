import numpy as np

def forward_radon_freq(m, dt, h, q, N, flow, fhigh):
    """
    Frequency domain algorithm for forward linear and parabolic Radon transform.

    Parameters:
    m: Radon panel matrix, size (nt, np), where nt is the number of time samples, and np is the number of Radon parameters.
    dt: Time sampling interval (in seconds).
    h: Offset or position of each trace (in meters), with length nh.
    q: Radon parameter:
       - When N=1, it represents the ray parameter;
       - When N=2, it represents the parabolic curvature.
    N: Type of Radon transform:
       - N=1 represents linear τ-p transform;
       - N=2 represents parabolic τ-p transform.
    flow: Starting frequency (Hz).
    fhigh: Ending frequency (Hz).

    Returns:
    d: Seismic data matrix, size (nt, nh).
    """
    nt, _ = m.shape  # Number of time samples and number of parameters in the Radon panel
    nh = len(h)       # Number of offset traces

    # If performing parabolic Radon transform, normalize the squared offset
    if N == 2:
        h = (h / np.max(np.abs(h)))**2

    # Compute the number of points for Fourier transform
    nfft = 2 * (2**np.ceil(np.log2(nt)).astype(int))

    # Perform frequency domain transformation of the Radon panel
    M = np.fft.fft(m, n=nfft, axis=0)  # Fourier transform of Radon panel
    D = np.zeros((nfft, nh), dtype=complex)  # Initialize frequency domain seismic data matrix

    # Compute indices corresponding to the starting and ending frequencies
    ilow = int(np.floor(flow * dt * nfft)) + 1
    ilow = max(ilow, 0)
    ihigh = int(np.floor(fhigh * dt * nfft)) + 1
    ihigh = min(ihigh, nfft // 2)

    # Process each frequency
    for ifreq in range(ilow, ihigh + 1):
        # Compute the current frequency (radians)
        f = 2 * np.pi * (ifreq - 1) / nfft / dt

        # Construct Radon transform matrix
        L = np.exp(1j * f * np.outer(np.conj(h).T, q))  # Size (nh, nq)

        # Extract current frequency component from Radon panel
        x = np.conj(M[ifreq, :]).T  # Radon panel values corresponding to the current frequency
        
        # Perform Radon transform
        _, n2 = L.shape # To be modified
        y = L @ x[:n2]  # Compute seismic data frequency component

        # Store positive frequency component and conjugate symmetric component
        D[ifreq, :] = np.conj(y).T  # Positive frequency part
        D[nfft + 1 - ifreq, :] = y.T  # Negative frequency part (conjugate symmetry)

    # Handle DC component (mid-frequency component)
    D[nfft // 2 + 1, :] = 0  # Set DC frequency component to zero

    # Perform inverse Fourier transform on seismic data
    d = np.fft.ifft(D, axis=0).real  # Inverse transform back to time domain
    d = d[:nt, :]  # Trim to original time length

    return d

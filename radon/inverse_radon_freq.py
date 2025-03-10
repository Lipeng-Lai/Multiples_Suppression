import numpy as np

def inverse_radon_freq(d, dt, h, q, N, flow, fhigh, mu, sol):
    """
    INVERSE_RADON_FREQ: Inverse linear or parabolic Radon transform.
                        Frequency domain algorithm.

    Parameters
    ----------
    d : ndarray
        Seismic data matrix (nt, nh), where nt is the number of time samples
        and nh is the number of traces. 
    dt : float
        Time sampling interval in seconds. 
    h : ndarray
        Offset or position of traces (length nh).
    q : ndarray
        Radon parameters (length nq). 
        - If N=1, linear ray parameters. 
        - If N=2, residual moveout at far offset. 
        # 变换参数q
        example:
        [nt, nx] = size(record);
        dq = (qmax - qmin) / (nq - 1); 
        q = qmin + dq * [0:1:nq - 1];
    N : int
        Radon transform type. Radon变换类型:
        - 1: Linear Radon transform (linear τ-p domain).
        - 2: Parabolic Radon transform (parabolic τ-p domain).
    flow : float
        Starting frequency in Hz (> 0).
    fhigh : float
        Ending frequency in Hz (< Nyquist frequency). 
    mu : float
        Regularization parameter for stability.
    sol : str
        Solution method: 解法选择。
        - 'ls': Least-squares solution.
        - 'adj': Adjoint solution.

    Returns
    -------
    m : ndarray
        Radon panel (linear or parabolic τ-p domain).
    """
    # seismic shape
    nt, nh = d.shape
    nq = len(q)

    # if parabolic, offset to normalization
    if N == 2:
        h = h / np.max(np.abs(h))

    # calculate fft points
    nfft = 2 * (2 ** int(np.ceil(np.log2(nt))))

    # sure the begin and end frequency
    ilow = int(np.floor(flow * dt * nfft)) + 1
    ilow = max(ilow, 2)

    ihigh = int(np.floor(fhigh * dt * nfft)) + 1
    ihigh = min(ihigh, nfft // 2 + 1)

    # seismic data: d(t,x) transform to D(x, f)
    D = np.fft.fft(d, n=nfft, axis=0)
    
    # initialize the radon panel
    M = np.zeros((nfft, nq), dtype=complex)
    
    # initialize the regularization
    Q = np.eye(nq) * nh

    # loop for calculate every frequency dot
    for ifreq in range(ilow, ihigh+1):
        # radian
        f = 2 * np.pi * (ifreq - 1) / nfft / dt

        # structure D = L * M
        # (h ** N)[:, None] =  Matlab (h .^ N)'
        # q[None, :] = Matlab (, q) -> (1, q) 
        L = np.exp(1j * f * (h ** N)[:, None] @ q[None, :])

        y = np.conj(D[ifreq, :]).T
        
        xa = L.conj().T @ y # xa = L^H * y
        
        A = L.conj().T @ L + mu * Q # L^H * L + regular

        if sol == 'ls':
            x = np.linalg.solve(A, xa)
            # x = A \ xa -> Ax = xa
            # x = (L^H * L + regular)^{-1} * (L^H * y)
        elif sol == 'adj':
            x = xa
        else:
            raise ValueError("Unknown solution method. Use 'ls' or 'adj'.")

        M[ifreq, :] = np.conj(x).T # postive
        M[nfft - ifreq + 2, :] = x.T # negative

    # set the intermediate frequency component to zero
    M[nfft // 2 + 1, :] = 0

    # f-x transform to t-x
    m = np.fft.ifft(M, axis=0).real

    # extract the original time sampling length
    m = m[:nt, :]

    return m


def inverse_radon_freq_lambda(d, dt, h, qmin, qmax, N, flow, fhigh, mu, sol):
    
    nt, nh = d.shape
    
    if N == 2:
        h = h / np.max(np.abs(h))
    
    nfft = 2 * (2 ** int(np.ceil(np.log2(nt))))

    ilow = int(np.floor(flow * dt * nfft)) + 1
    ilow = max(ilow, 2)

    ihigh = int(np.floor(fhigh * dt * nfft)) + 1
    ihigh = min(ihigh, nfft // 2 + 1)

    freq = np.arange(ilow, ihigh + 1)
    
    nq = len(freq)
    dq = (qmax - qmin) / (nq - 1)
    q = qmin + dq * np.arange(nq)
    
    f = (freq - 1) * 2 * np.pi / nfft / dt
    
    lamb = q * f
    
    D = np.fft.fft(d, n=nfft, axis=0)
    Q = np.eye(nq) * nh

    L = np.zeros((nh, nq), dtype=complex)
    
    for x_t in range(nh):
        for q_t in range(nq):
            L[x_t, q_t] = np.exp(1j * lamb[q_t] * (h[x_t] ** N))

    
    y = np.conj(D[ilow: ihigh+1, :]).T
    
    xa = L.conj().T @ y
    
    A = L.conj().T @ L + mu * Q
    
    x = np.linalg.solve(A, xa)
    
    M = np.zeros((nfft, nq), dtype=complex)
    
    M[freq, :] = np.conj(x).T
    M[nfft - freq + 2, :] = x.T

    M[nfft // 2 + 1, :] = np.zeros((1, nq))
    
    m = np.fft.ifft(M, axis=0).real
    m = m[:nt, :]

    return m, M
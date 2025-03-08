import numpy as np

def nmo(d, dt, h, tnmo, vnmo, max_stretch):
    """
    NMO: A program for NMO correction.

    Parameters:
    d (ndarray): Input data of shape (nt, nh) (time, offsets).
    dt (float): Sampling interval in seconds.
    h (ndarray): Offsets in meters.
    tnmo (list): Intercept times in seconds.
    vnmo (list): NMO velocities in m/s.
    max_stretch (float): Maximum stretch allowed in percentage.

    Returns:
    dout (ndarray): Data after NMO correction.
    M (ndarray): Number of x-samples that survived muting at each time position.
    ti (ndarray): Interpolated time vector.
    vi (ndarray): Interpolated velocity vector.
    """
    nt, nh = d.shape

    # Interpolate t0, v pairs
    N = len(vnmo)
    if N > 1:
        t1 = np.array([0] + tnmo + [(nt - 1) * dt])
        v1 = np.array([vnmo[0]] + vnmo + [vnmo[-1]])

        ti = np.arange(0, nt) * dt
        vi = np.interp(ti, t1, v1)
    else:
        ti = np.arange(0, nt) * dt
        vi = np.full(nt, vnmo[0])

    dout = np.zeros_like(d)
    M = np.zeros(nt)

    for it in range(nt):
        for ih in range(nh):
            arg = ti[it]**2 + (h[ih] / vi[it])**2
            time = np.sqrt(arg)
            stretch = (time - ti[it]) / (ti[it] + 1e-10)

            if stretch < max_stretch / 100:
                M[it] += 1

                its = time / dt
                it1 = int(np.floor(its))
                it2 = it1 + 1
                a = its - it1

                if it2 < nt:
                    dout[it, ih] = (1 - a) * d[it1, ih] + a * d[it2, ih]

    return dout, M, ti, vi

import pylops
import numpy as np

def SRME(data_mutter, p, dt, dx):
    # data_mutter (ns, nr, nt)
    nx = data_mutter.shape[1]
    nt = data_mutter.shape[2]
    
    data_p = np.concatenate((data_mutter, np.zeros((data_mutter.shape[0],
                                                data_mutter.shape[1],
                                                data_mutter.shape[2] // p))), axis=-1)
    
    data_fft = np.fft.rfft(data_p, axis=-1)
    data_fft = data_fft.transpose(2, 0, 1).astype(np.complex64)
    
    MDCop = pylops.waveeqprocessing.MDC(
        data_fft,               
        nt=data_p.shape[2],     
        nv=data_p.shape[0],
        dt=dt,
        dr=dx,
        twosided=False,
    )
    
    multiples = MDCop @ data_p.transpose(2, 1, 0).ravel()
    
    multiples = multiples.reshape(data_p.shape[2], nx, nx)
    multiples = multiples.transpose(1, 2, 0)
    multiples = multiples[:,:, :nt]
    
    # normalization
    multiples *= -0.8 * np.max(np.abs(data_mutter)) / np.max(abs(multiples))
    
    # shift ...
    
    return multiples


def shift_up(arr, n):
    _, _, rows = arr.shape
    result = np.zeros_like(arr)
    result[:, :, :rows-n] = arr[:, :, n:]
    return result
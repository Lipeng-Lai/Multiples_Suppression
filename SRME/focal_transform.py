import numpy as np
from tqdm import tqdm

def focal_transform_forward(data, primaries, epsilon=1e-1, iter=2):
    # data(ns, nr, nt)
    # primaries(ns, nr, nt)    
    assert data.shape == primaries.shape
    
    nt = data.shape[2]
    nr = data.shape[1]
    
    nfft = iter * (2 ** int(np.ceil(np.log2(nt))))
    
    P = np.fft.rfft(data, n=nfft, axis=2)
    
    Delta_P = np.fft.rfft(primaries, n=nfft, axis=2)
    
    Q = np.zeros_like(P)
    I = np.eye(nr)
    
    for f_idx in tqdm(range(P.shape[2]), desc="forward", ncols=100):
        Delta_P_idx_f = Delta_P[:, :, f_idx] # Delta P
        P_idx_f = P[:, :, f_idx] # P
        
        term1 = np.conj(Delta_P_idx_f).T
        
        # B = \Delta P x \Delta P^H + \epsilon^2 I
        term2 = np.dot(Delta_P_idx_f, np.conj(Delta_P_idx_f).T) + epsilon**2 * I
        
        # B = term3
        term3 = np.linalg.inv(term2)
        
        # F = np.dot(term1, term3)
        Q[:, :, f_idx] = np.dot(np.dot(term1, term3), P_idx_f)
        
    q = np.fft.irfft(Q, axis=2)
    q = q[:, :, :nt]
    
    return q


def focal_transform_inverse(focal_domain_data, primaries, epsilon=1e-1, iter=2):
    assert focal_domain_data.shape == primaries.shape
    
    nt = focal_domain_data.shape[2]
    nr = focal_domain_data.shape[1]
    
    nfft = iter * (2 ** int(np.ceil(np.log2(nt))))
    
    P = np.fft.rfft(focal_domain_data, n=nfft, axis=2)
    
    Delta_P = np.fft.rfft(primaries, n=nfft, axis=2)
    
    Q = np.zeros_like(P)
    I = np.eye(nr)
    
    for f_idx in tqdm(range(P.shape[2]), desc="inverse", ncols=100):
        Delta_P_idx_f = Delta_P[:, :, f_idx] # Delta P
        P_idx_f = P[:, :, f_idx] # P
        
        term1 = np.conj(Delta_P_idx_f.T)
        
        term2 = np.dot(Delta_P_idx_f, np.conj(Delta_P_idx_f.T)) + epsilon**2 * I
        
        term3 = np.linalg.inv(term2)
        # F = np.dot(term1, term3)
        F = np.dot(term1, term3)
        
        G = np.linalg.inv(F)
        Q[:, :, f_idx] = np.dot(G, P_idx_f)
        
    q = np.fft.irfft(Q, axis=2)
    q = q[:, :, :nt]
    
    return q


def focal_transform_forward_ls(data, primaries, epsilon=1e-1, iter=2):
    
    assert data.shape == primaries.shape
    
    nt = data.shape[2]
    nr = data.shape[1]
    
    nfft = iter * (2 ** int(np.ceil(np.log2(nt))))
    
    P = np.fft.rfft(data, n=nfft, axis=2)
    
    Delta_P = np.fft.rfft(primaries, n=nfft, axis=2)
    
    Q = np.zeros_like(P)
    Delta_Q = np.zeros_like(P)
    I = np.eye(nr)
    
    for f_idx in tqdm(range(P.shape[2]), desc="forward_ls", ncols=100):
        Delta_P_idx_f = Delta_P[:, :, f_idx] # Delta P
        P_idx_f = P[:, :, f_idx] # P
        
        term1 = np.conj(Delta_P_idx_f.T)
        
        term2 = np.dot(Delta_P_idx_f, np.conj(Delta_P_idx_f.T)) + epsilon**2 * I
        
        term3 = np.linalg.inv(term2)
        # F = np.dot(term1, term3)
        Q[:, :, f_idx] = np.dot(np.dot(term1, term3), P_idx_f)
        
        term1_ls = np.conj(P_idx_f.T)
        term2_ls = np.linalg.inv(np.dot(term1_ls, P_idx_f) + epsilon**2 * I) 
        term3_ls = np.dot(term2_ls, term1_ls)
        F_ls = np.dot(term3_ls, Q[:, :, f_idx])
        
        # Delta Q = Q - np.dot(F_ls, P)
        Delta_Q[:, :, f_idx] = Q[:, :, f_idx] - np.dot(F_ls, P_idx_f)
    
    q = np.fft.irfft(Q, axis=2)
    q = q[:, :, :nt]
    
    delta_q = np.fft.irfft(Delta_Q, axis=2)
    delta_q = delta_q[:, :, :nt]
    
    return q, delta_q
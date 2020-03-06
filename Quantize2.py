import numpy as np


def quantize(signal, nbit, level):
    
    tt= (2**(nbit-1))-1
    
    nn = np.round(signal/level)
    
    nn[nn>tt] = tt
    nn[nn<-tt] = -tt
    sigQ = nn * level
    eRROR = sigQ - signal
    
    return sigQ, eRROR

"""
def _quantize1(signal, nbit, level):
    
    #if set_level:
    #   level = 0.5
    #else:
    #level = set_level

    tt = (2*(nbit -1))-1
    nn = np.round(signal/level)

    nn[nn>tt] = tt
    nn[nn<tt] =-tt
    
    SigQ =  nn * level
    eRROR = SigQ - signal
        
    return SigQ, eRROR
"""

def quantize2(signal, nbit):
    
    xi = np.round(nbit * signal)
    XQ = (1/nbit) * xi
    error = XQ - signal
        
    return XQ, error


def SNR(signal, noise):
    
    sig_vec = np.mean(signal**2, axis = 0)
    noi_vec = np.mean (noise**2, axis = 0)
    
    nsr = np.mean(sig_vec-noi_vec)/np.std(sig_vec)
    
    return nsr

---
title: "fourier transforms in different languages"
draft: true
---

```
def FFT_1D(gx, x):
    """
    Discrete fourier transform.
    
    Parameters
    ----------
    gx: `np.ndarray` 
        grid-function of position x (1d arr)
    x: `np.ndarray`
        positions
        
    Returns
    -------
    gx: `np.ndarray`
        fft of gx -> gk
    """
    dx = x[1] - x[0] # length of a cell
    L = x[-1] + dx # total length of system, including the last point
    N = len(x) # number of cells
    
    k = kn(N, L)
    gk = np.zeros((N), dtype=np.complex_)
    # summing over the wavenumber
    for i in range(N):
        gk[i] = dx * np.trapz(gx * np.exp(- 1j * k[i] * x)) 
    
    return (gk, k)

def IFFT_1D(gk, k, x):
    '''
    inverse fourier transform
    
    Parameters
    ----------
    gk: `np.ndarray`
        function in k space
    k: `np.ndarray`
        wavenumbers
    x: `np.ndarray`
        positions
    
    Returns
    -------
    gx: `np.ndarray`
        function in position space
    '''
    dx = x[1]-x[0] # length of a cell
    L = x[-1]+dx # total length of system, including the last point
    N = len(x) # number of cells
    
    gx = np.zeros((N), dtype=np.complex_)
    
    # summing over x
    for i in range(N):
        # using np.nansum because `gk` will have a nan value at k = 0
        gx[i] = (1 / L) * np.nansum(gk * np.exp(1j * k * x[i])) 
    
    return gx
```

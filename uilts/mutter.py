import numpy as np

# input = (t, x)
def mutter(input, x0, t0, k1, k2, delta):
    input = input.T
    n1, n2 = input.shape
    
    for i in range(n1):
        for j in range(n2):
            # Apply conditions for j < x0 (left side)
            if j < x0:
                if i - delta < np.floor(-k1 * j + t0 + k1 * x0):
                    input[i, j] = 0
            # Apply conditions for j >= x0 (right side)
            if j >= x0:
                if i - delta < np.floor(k2 * j + t0 - k2 * x0):
                    input[i, j] = 0
    
    return input.T
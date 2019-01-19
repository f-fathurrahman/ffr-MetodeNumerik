import numpy as np

def search_root_interval(f, a, b, dx, verbose=False):
    """
    Searches the interval (a, b) in increments dx for the 
    bounds (x1,x2) of the smallest root of f(x)
    Returns (None,None) if no roots were detected.
    """

    x1 = a
    f1 = f(x1)

    # start incrementing a by dx
    x2 = a + dx
    f2 = f(x2)

    while np.sign(f1) == np.sign(f2):

        # This should only happen when the root is not bracketed within
        # a and b
        if x1 >= b:
            return None, None

        # This is the case where the root is not found within this trial
        # interval

        # Change x1 to x2
        x1 = x2
        f1 = f2

        # new x2 by incrementing new x1
        x2 = x1 + dx
        f2 = f(x2)

    else:
        # We found the interval which contains a root
        if verbose:
            print("Found interval: (%18.10f,%18.10f)" % (x1,x2))
        return x1, x2

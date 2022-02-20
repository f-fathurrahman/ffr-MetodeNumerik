import numpy as np

def root_bairstow( a_, rr=1.0, ss=1.0, NiterMax=100, SMALL=1e-10 ):

    a = np.copy(a_) # do not modify the input

    Ndeg = len(a) - 1

    re = np.zeros(Ndeg)
    im = np.zeros(Ndeg)

    b = np.zeros(Ndeg+1)
    c = np.zeros(Ndeg+1)

    r = rr
    s = ss
    n = Ndeg
    ier = 0
    ea1 = 1.0
    ea2 = 1.0

    iteration = 0

    while True:

        #print("\nBefore factorization:")
        #print("n = ", n)
        #print("a = ", a)
        #print("b = ", b)
        #print("c = ", c)

        if (n < 3) or (iteration >= NiterMax) :
            #print("Exit main loop")
            break

        iteration = 0
        while True:
            #
            iteration = iteration + 1
            #
            if iteration >= NiterMax:
                #print("Break inner loop")
                break
            #
            #print("inner loop iter = ", iteration)
            #
            b[n] = a[n]
            b[n-1] = a[n-1] + r*b[n]
            #
            c[n] = b[n]
            c[n-1] = b[n-1] + r*c[n]
            #
            for i in range(n-2,-1,-1):
                b[i] = a[i] + r*b[i+1] + s*b[i+2]
                c[i] = b[i] + r*c[i+1] + s*c[i+2]
            #
            #print("a = ", a)
            #print("b = ", b)
            #print("c = ", c)
            det = c[2]*c[2] - c[3]*c[1]
            #print("det = ", det)
            #exit()
            #
            if abs(det) >= SMALL:
                dr = ( -b[1]*c[2] + b[0] * c[3] ) / det
                ds = ( -b[0]*c[2] + b[1] * c[1] ) / det
                r = r + dr
                s = s + ds
                #print("r = %18.10f s = %18.10f" % (r, s))
                #print("dr = %18.10f ds = %18.10f" % (dr, ds))
            else:
                # Update with different initial r and s
                print("INFO: Resetting iter")
                r = r + 1.0
                s = s + 1.0
                iteration = 0
            if (abs(dr) <= SMALL) and (abs(ds) <= SMALL):
                #print("INFO: Found roots")
                break
    
        r1, i1, r2, i2  = quadroot(r, s)
        re[n-1] = r1
        im[n-1] = i1
        re[n-2] = r2
        im[n-2] = i2
    
        n = n - 2
        for i in range(0,n+1):
            a[i] = b[i+2]
  
    if iteration < NiterMax:
        if n == 2:
            r = -a[1]/a[2]
            s = -a[0]/a[2]
            r1, i1, r2, i2 = quadroot(r, s)
            re[n-1] = r1
            im[n-1] = i1
            re[n-2] = r2
            im[n-2] = i2
        else:
            re[n-1] = -a[0]/a[1]
            im[n-1] = 0.0
    else:
        print("WARNING: there are some errors")
        ier = 1

    # Convert to complex-valued array
    zroots = np.zeros(Ndeg, dtype=np.complex128)
    # `numpy.complex128`: Complex number type composed
    # of 2 64-bit-precision floating-point numbers.
    for i in range(Ndeg):
        zroots[i] = complex(re[i], im[i])
    return zroots


def quadroot(r, s):
    disc = r**2 + 4.0*s
    if disc > 0:
        r1 = (r + np.sqrt(disc))/2.0
        r2 = (r - np.sqrt(disc))/2.0
        i1 = 0.0
        i2 = 0.0
    else:
        r1 = r/2.0
        r2 = r1
        i1 = np.sqrt(abs(disc))/2.0
        i2 = -i1
    return r1, i1, r2, i2

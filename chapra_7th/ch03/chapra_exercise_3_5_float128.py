# Using single precision

import numpy as np

N = 10000

fn_true = np.float128(np.pi)**np.float128(4)/np.float128(90)

# From 1 to 10000
fn_v1 = np.float128(0.0)
for i in range(1,N+1):
    fn_v1 = fn_v1 + np.float128(1.0)/np.float128(i)**np.float128(4)

# From 10000 to 1
fn_v2 = np.float128(0.0)
for i in range(N+1,0,-1):
    fn_v2 = fn_v2 + np.float128(1.0)/np.float128(i)**np.float128(4)

print("Using single-precision (np.float128)")

print("fn_true = ", fn_true)
print("fn_v1 = ", fn_v1)
print("fn_v2 = ", fn_v2)

print("Summing from 1 to 10000: error = ", abs(fn_v1 - fn_true))
print("Summing from 10000 to 1: error = ", abs(fn_v2 - fn_true))

print(type(fn_true))
print(type(fn_v1))
print(type(fn_v2))


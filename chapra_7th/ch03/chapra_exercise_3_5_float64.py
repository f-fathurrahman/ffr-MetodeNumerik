# Using single precision

import numpy as np

N = 10000

fn_true = np.pi**4/90

# From 1 to 10000
fn_v1 = 0.0
for i in range(1,N+1):
    fn_v1 = fn_v1 + 1.0/i**4

# From 10000 to 1
fn_v2 = 0.0
for i in range(N+1,0,-1):
    fn_v2 = fn_v2 + 1.0/i**4

print("Using double-precision (np.float64)")

print("fn_true = ", fn_true)
print("fn_v1 = ", fn_v1)
print("fn_v2 = ", fn_v2)

print("Summing from 1 to 10000: error = ", abs(fn_v1 - fn_true))
print("Summing from 10000 to 1: error = ", abs(fn_v2 - fn_true))

print(type(fn_true))
print(type(fn_v1))
print(type(fn_v2))


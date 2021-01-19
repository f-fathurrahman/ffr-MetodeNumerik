def root_newton_raphson(f, df, x0, TOL=1e-10, NiterMax=100, verbose=False):
    x = x0
    for i in range(1,NiterMax+1):
        xnew = x  - f(x)/df(x) # Newton-Raphson formula
        fxnew = f(xnew) # should be zero when xnew is a root
        if verbose:
            print("%3d %18.10f %18.10e" % (i, xnew, fxnew))
        if abs(fxnew) <= 1e-10:
            if verbose:
                print("Convergence is achieved")
            break
        x = xnew
    return x

R = 0.082054 # L atm /(mol K)

# T in K, p in atm
# V is returned in L
def calc_V_ideal_gas(T, p):
    return R*T/p

# vdw eq, RHS is zero
def vdw_eq(a, b, T, p, v):
    return (p + a/v**2)*(v - b) - R*T

def d_vdw_eq(a, b, T, p, v):
    return p - a/v**2 + 2*a*b/v**3


# using Newton-Raphson
def solve_for_v(a, b, T, p, v0=1.0, verbose=False):
    f = lambda v: vdw_eq(a, b, T, p, v)
    df = lambda v: d_vdw_eq(a, b, T, p, v)
    v = root_newton_raphson(f, df, v0, verbose=verbose)
    return v

# Constants for vdW equation

a_CO2 = 3.592
b_CO2 = 0.04267

a_O2 = 1.360
b_O2 = 0.03183

p_list = [1.0, 10.0, 100.0]

# Test for some values using verbose=True
#T = 300.0
#v = solve_for_v(a_O2, b_O2, T, p_list[1], verbose=True) 

print("")
print("p in atm, v in L/mol")
print("")

T = 300.0
print("T = %4.0f K" % T)
print("p     v_ideal_gas    v_CO2       v_O2")
print("----  -----------    -----       ----")
for p in p_list:
    v_ideal_gas = calc_V_ideal_gas(T, p)
    v_CO2 = solve_for_v(a_CO2, b_CO2, T, p)
    v_O2 = solve_for_v(a_O2, b_O2, T, p)
    print("%4.0f %10.4f  %10.4f  %10.4f" % (p, v_ideal_gas, v_CO2, v_O2))


print()
T = 500.0
print("T = %4.0f K" % T)
print("p     v_ideal_gas    v_CO2       v_O2")
print("----  -----------    -----       ----")
for p in p_list:
    v_ideal_gas = calc_V_ideal_gas(T, p)
    v_CO2 = solve_for_v(a_CO2, b_CO2, T, p)
    v_O2 = solve_for_v(a_O2, b_O2, T, p)
    print("%4.0f %10.4f  %10.4f  %10.4f" % (p, v_ideal_gas, v_CO2, v_O2))


print()
T = 700.0
print("T = %4.0f K" % T)
print("p     v_ideal_gas    v_CO2       v_O2")
print("----  -----------    -----       ----")
for p in p_list:
    v_ideal_gas = calc_V_ideal_gas(T, p)
    v_CO2 = solve_for_v(a_CO2, b_CO2, T, p)
    v_O2 = solve_for_v(a_O2, b_O2, T, p)
    print("%4.0f %10.4f  %10.4f  %10.4f" % (p, v_ideal_gas, v_CO2, v_O2))

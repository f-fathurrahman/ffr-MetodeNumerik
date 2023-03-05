import math, cmath

print("Kelompok XXX")
print("NIM1 Nama1")

a = float(input("Masukkan koefisien a: "))
b = float(input("Masukkan koefisien b: "))
c = float(input("Masukkan koefisien c: "))

D = b**2 - 4*a*c
akar_kompleks = False
if D >= 0:
    akar_kompleks = False
    x1 = (-b + math.sqrt(D))/(2*a)
    x2 = (-b - math.sqrt(D))/(2*a)
else:
    akar_kompleks = True
    x1 = (-b + cmath.sqrt(D))/(2*a)
    x2 = (-b - cmath.sqrt(D))/(2*a)

if akar_kompleks:
    print("Akar-akar kompleks")
else:
    print("Akar-akar real")

print("x1 = ", x1)
print("x2 = ", x2)

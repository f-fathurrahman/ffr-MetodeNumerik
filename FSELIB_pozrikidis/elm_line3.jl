"""
 FSELIB

 Symmetric discretization of a line segment
 into a graded mesh of n elements
 subtended between the left end-point x1
 and the right end-point x2

 n is odd: 1,3,5,...

 The variable "ratio" is the ratio
 of the size of the mid-elements
 to the size of the first element
"""
function elm_line3( x1,x2, N, ratio )

  if N % 2 != 1
    @printf("ERROR: N must be odd. This N = %d", N)
    exit()
  end

  if N==1
    xe = zeros(2)
    xe[1] = x1
    xe[2] = x2
    return xe
  end

  if ratio==1
    α  = 1.0
    factor = 1.0/N
  else
    texp = 2.0/(N-1.0)
    α = ratio^texp
    tmp1 = (N + 1.0)/2.0
    tmp2 = (N - 1.0)/2.0
    factor = (1.0 - α) / (2.0 - α^tmp1 - α^tmp2)
  end

  xe = zeros(N+1)

  Δx = (x2-x1) * factor   # length of first element

  xe[1] = x1    # first point

  for i = 2:round(Int64,(N+3)/2)
    xe[i] = xe[i-1] + Δx
    Δx = Δx*α
  end

  Δx = Δx/α^2

  for i = round(Int64,(N+5)/2):N+1
    xe[i] = xe[i-1] + Δx
    Δx = Δx/α
  end

  return xe

end

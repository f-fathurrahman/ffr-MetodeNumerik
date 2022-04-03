"""
 FSELIB

 Symmetric discretization of a line segment
 into a graded mesh of n elements
 subtended between the left end-point x1
 and the right end-point x2

 n is 1 or even: 2,4,6,...

 The variable "ratio" is the ratio
 of the size of the mid
 elements to the size of the first element
"""

function elm_line2( x1, x2, N::Int64, ratio )

  if N == 1
    xe = zeros(2)
    xe[1] = x1
    xe[2] = x2
    return xe
  end


  if N==2
    xe = zeros(3)
    xe[1] = x1
    xe[2] = 0.5*(x1 + x2)
    xe[3] = x2
    return xe
  end

  xe = zeros(N+1)
  Nh = Int64(N/2)
  xh = 0.5*( x1 + x2 )   # mid-point

  if ratio==1
    α = 1.0;
    factor = 1.0/Nh
  else
    texp   = 1.0 / (Nh - 1.0)
    α = ratio^texp
    factor = (1.0-α) / (1.0 - α^Nh)
  end

  # length of first element
  Δx = (xh-x1) * factor

  # first point
  xe[1] = x1
  for i = 2:Nh+1
    xe[i] = xe[i-1] + Δx
    Δx = Δx*α
  end

  Δx = Δx/α

  for i = Nh+2:N+1
    xe[i]  = xe[i-1] + Δx
    Δx = Δx/α
  end

  return xe

end

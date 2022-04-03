"""
 FSELIB

 Unsymmetrical discretization of a line segment
 into a graded mesh of n elements
 subtended between the left end-point x1
 and the right end-point x2

 The variable "ratio" is the ratio
 of the length of the last
 element to the length of the first element.

 xe: element end-nodes
"""

function elm_line1( x1, x2, N, ratio )

  if N == 1
    xe = zeros(2)
    xe[1] = x1
    xe[2] = x2
    return xe
  end

  if ratio == 1
    α = 1.0
    factor = 1.0/N
  else
    texp = 1/(N-1)
    α = ratio^texp
    factor = (1.0 - α)/(1.0 - α^N)
  end

  # length of the first element
  Δx = (x2 - x1)*factor

  xe = zeros(N+1)
  xe[1] = x1
  for i = 2:N+1
    xe[i] = xe[i-1] + Δx
    Δx = α*Δx
  end

  return xe

end

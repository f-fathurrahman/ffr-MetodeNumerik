PROGRAM example_18_7
  IMPLICIT NONE 
  REAL(8) :: t(5), v(5)
  REAL(8) :: tt, vv

  t = (/ 1.d0, 3.d0, 5.d0, 7.d0, 13.d0 /)
  v = (/ 800.d0, 2310.d0, 3090.d0, 3940.d0, 4755.d0 /)

  tt = 10.d0
  CALL lagrange_interp( 5, t, v, tt, vv )
  WRITE(*,*) 'vv = ', vv
END PROGRAM 


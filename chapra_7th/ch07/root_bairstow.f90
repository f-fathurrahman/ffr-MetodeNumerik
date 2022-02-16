subroutine root_bairstow(a,nn,es,rr,ss,maxit,re,im,ier)
  implicit none
  integer :: nn
  real(8) :: a(0:nn)
  real(8) :: rr, ss
  real(8) :: es
  integer :: maxit
  real(8) :: re(nn)
  real(8) :: im(nn)
  integer :: ier
  ! Local
  real(8), allocatable :: b(:), c(:)
  real(8), parameter :: SMALL=1.d-12
  integer :: n, i, iter
  real(8) :: ea1, ea2
  real(8) :: det, dr, ds
  real(8) :: i1, i2, s, r, r1, r2

  allocate( b(0:nn) )
  allocate( c(0:nn) )

  r = rr
  s = ss
  n = nn
  ier = 0
  ea1 = 1.d0
  ea2 = 1.d0

  iter = 0
  write(*,*) 'iter = ', iter
  write(*,*) 'maxit = ', maxit
  write(*,*) 'n = ', n

  do

    if( (n < 3) .or. (iter >= maxit) ) then
      write(*,*) 'iter = ', iter
      write(*,*) 'Exit main loop'
      exit
    endif
    
    iter = 0
    do
      iter = iter + 1
      write(*,*) 'inner loop iter = ', iter
      b(n) = a(n)
      b(n-1) = a(n-1) + r*b(n)
    
      c(n) = b(n)
      c(n - 1) = b(n-1) + r*c(n)
      do i = n-2,0,-1
        b(i) = a(i) + r*b(i+1) + s*b(i+2)
        c(i) = b(i) + r*c(i+1) + s*c(i+2)
      enddo
    
      det = c(2)*c(2) - c(3)*c(1)
      if(abs(det) >= SMALL) then
        dr = ( -b(1)*c(2) + b(0) * c(3) ) / det
        ds = ( -b(0)*c(2) + b(1) * c(1) ) / det
        r = r + dr
        s = s + ds
        !
        if(abs(r) >= SMALL) then
          ea1 = ABS(dr/r) * 100
        endif
        !
        if(abs(s) >= SMALL) then
          ea2 = ABS(ds/s) * 100
        endif
      else
        ! what's this?
        write(*,*) 'Resetting iter'
        r = r + 1.d0
        s = s + 1.d0
        iter = 0
      endif
    
      IF( (ea1 <= es) .AND. (ea2 <= es) ) then
        exit
      endif
    enddo
    
    call Quadroot(r,s,r1,i1,r2,i2)
    re(n) = r1
    im(n) = i1
    re(n-1) = r2
    im(n-1) = i2
    
    n = n - 2
    do i = 0,n
      a(i) = b(i+2)
    enddo
  
  enddo


  if( iter < maxit ) then
    if(n == 2) then
      r = -a(1)/a(2)
      s = -a(0)/a(2)
      call quadroot(r,s,r1,i1,r2,i2)
      re(n) = r1
      im(n) = i1
      re(n-1) = r2
      im(n-1) = i2
    else
      re(n) = -a(0)/a(1)
      im(n) = 0.d0
    endif
  else
    ier = 1
  endif

  deallocate(b)
  deallocate(c)
end subroutine 

subroutine quadroot(r,s,r1,i1,r2,i2)
  implicit none
  real(8) :: r,s,r1,i1,r2,i2
  ! local
  real(8) :: disc
  disc = r**2 + 4.d0*s
  
  if(disc > 0) then
    r1 = (r + SQRT(disc))/2.d0
    r2 = (r - SQRT(disc))/2.d0
    i1 = 0.d0
    i2 = 0.d0
  else
    r1 = r/2.d0
    r2 = r1
    i1 = SQRT(ABS(disc))/2.d0
    i2 = -i1
  endif
  return
end subroutine


program test
  implicit none
  integer, parameter :: nn = 5
  real(8), allocatable :: a(:), re(:), im(:)
  real(8), parameter :: es=1e-10
  real(8) :: rr, ss
  integer, parameter :: maxit = 100
  integer :: ier
  integer :: i

  allocate(a(0:nn))
  allocate(re(nn))
  allocate(im(nn))

  ! coefficients are a0, a1, ..., aN
  ! N is the degree of the polinomial
  a(:) = (/ 1.25d0, -3.875d0, 2.125d0, 2.75d0, -3.5d0, 1.d0 /)
  re(:) = 0.d0
  im(:) = 0.d0

  rr = 0.d0
  ss = 0.d0

  call root_bairstow(a,nn,es,rr,ss,maxit,re,im,ier)

  do i = 1,nn
    write(*,'(1x,I4,2F18.10)') i, re(i), im(i)
  enddo

  deallocate(a)
  deallocate(re)
  deallocate(im)

end program

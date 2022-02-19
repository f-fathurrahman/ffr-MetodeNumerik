subroutine polydiv(a, n, d, m, q, r)
  implicit none
  integer :: n
  real(8) :: a(0:n)
  integer :: m
  real(8) :: d(0:m)
  real(8) :: q(0:n)
  real(8) :: r(0:n)
  ! local
  integer :: i, j, k
  
  do j = 0, n
    !write(*,*) 'j = ', j
    r(j) = a(j)
    q(j) = 0.d0
  enddo

  do k = n-m,0,-1
    !write(*,*) 'k = ', k
    q(k+1) = r(m+k)/d(m)
    do j = m+k-1, k, -1
      write(*,*) 'j = ', j
      !r(j) = r(j) - q(k+1) * b(j-k)
      r(j) = r(j) - q(k+1) * d(j-k)
    enddo
  enddo
  
  do j = m, n
    write(*,*) '2nd part j = ', j
    r(j) = 0.d0
  enddo

  !do i = 0, n-m
  !  a(i) = q(i+1)
  !enddo
end subroutine



!-----------
program test
!-----------
  implicit none

  integer, parameter :: n = 2
  real(8) :: a(0:n)
  integer, parameter :: m = 1
  real(8) :: d(0:m)
  real(8) :: q(0:n)
  real(8) :: r(0:n)

  a = (/ -24.d0, 2.d0, 1.d0 /)
  d = (/ -4.d0, 1.d0 /)

  call polydiv(a, n, d, m, q, r)

  write(*,*) 'after polydiv:'
  write(*,*) 'a = ', a
  write(*,*) 'd = ', d
  write(*,*) 'q = ', q
  write(*,*) 'r = ', r

end program
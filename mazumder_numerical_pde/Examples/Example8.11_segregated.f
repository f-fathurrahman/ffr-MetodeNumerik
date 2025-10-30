! Segregated Solution of 3 PDEs 
! Example 8.11

  program main

    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 3
    INTEGER(4), PARAMETER :: nx = 101, ny = 101
    REAL(8), PARAMETER :: l_x = 1.0d0, l_y = 1.0d0
    INTEGER(4), PARAMETER :: max_out=100000, max_in = 1000
    INTEGER(4) :: i,j,it_in,it_out,k,kl
    REAL(8), PARAMETER :: one=1.0d0,half=0.5d0,zero=0.0d0,two=2.0d0
    REAL(8), PARAMETER :: tol_in = 1.0D-3, tol_out = 1.0D-6
    REAL(8), DIMENSION(nx) :: x,ax,bx,cx,dx,sx
    REAL(8), DIMENSION(ny) :: y,ay,by,cy,dy,sy
    REAL(8), DIMENSION(n,n) :: gij
    REAL(8), DIMENSION(n) :: res
    REAL(8) :: phi(n,nx,ny)
    REAL(8) :: sums,sumk,resk,d_x,d_y,dx2,dy2,ao,aw,ae,as,an
    REAL(4) :: ta(2),time_start,time_end
    LOGICAL :: converge

! Open Output file
    OPEN(unit=20,file="Example8_11_segregated.out",status="unknown")
    OPEN(unit=30,file="Example8_11_segregated.rsl",status="unknown")

! Geometry
    d_x = l_x/(nx-1)
    d_y = l_y/(ny-1)
    dx2 = d_x*d_x
    dy2 = d_y*d_y
    DO i = 1,nx
      x(i) = (i-1)*d_x
    ENDDO
    DO j = 1,ny
      y(j) = (j-1)*d_y
    ENDDO

! Initial Guess for all phi
    phi(:,:,:) = half

! Boundary Conditions
! Left
    phi(1,1,:) = one; phi(2,1,:) = zero; phi(3,1,:) = zero
! Right
    phi(1,nx,:) = zero; phi(2,nx,:) = zero; phi(3,nx,:) = one
! Bottom
    phi(1,:,1) = zero; phi(2,:,1) = one; phi(3,:,1) = zero
! Top
    phi(1,:,ny) = zero; phi(2,:,ny) = zero; phi(3,:,ny) = one

! Set Gamma
    gij(:,:) = one
    gij(1,2) = 0.4d0; gij(1,3) = -0.2d0
    gij(2,1) = 0.4d0; gij(2,3) = -0.3d0
    gij(3,1) = -0.2d0; gij(3,2) = -0.3d0

! Link coefficients
    ao = -(two/dx2 + two/dy2)
    ae = one/dx2
    aw = one/dx2
    an = one/dy2
    as = one/dy2

    time_start = etime(ta)

! Start of Global Iteration Loop For Phi
    DO it_out = 1,max_out

! Loop over all three dependent variables
      DO k = 1,n

! Inner Iteration
        DO it_in = 1,max_in

! Row-wise sweep

! Interior Rows
          DO j = 2,ny-1
            
            ax(:)=zero;bx(:)=zero;cx(:)=zero;dx(:)=zero;sx(:)=zero
            dx(1) = one; bx(1) = phi(k,1,j)
            DO i = 2,nx-1
               dx(i) = ao*gij(k,k)
               cx(i) = ae*gij(k,k)
               ax(i-1) = aw*gij(k,k)
               bx(i) = (-an*phi(k,i,j+1)-as*phi(k,i,j-1))*gij(k,k)
               sums = zero
               DO kl = 1,n
                 IF(k == kl)CYCLE
                 sums = sums - gij(k,kl)*((phi(kl,i+1,j)-two*phi(kl,i,j)+phi(kl,i-1,j))/dx2 &
                             + (phi(kl,i,j+1)-two*phi(kl,i,j)+phi(kl,i,j-1))/dy2)
               ENDDO
               bx(i) = bx(i) + sums
            ENDDO
            dx(nx) = one; bx(nx) = phi(k,nx,j)

            CALL TRI(nx,ax,dx,cx,bx,sx)

            phi(k,:,j) = sx(:)

          ENDDO  ! End Row-wise sweep

! Column-wise sweep

! Interior Columns
          DO i = 2,nx-1
            
            ay(:)=zero;by(:)=zero;cy(:)=zero;dy(:)=zero;sy(:)=zero
            dy(1) = one; by(1) = phi(k,i,1)
            DO j = 2,ny-1
               dy(j) = ao*gij(k,k)
               cy(j) = an*gij(k,k)
               ay(j-1) = as*gij(k,k)
               by(j) = (-ae*phi(k,i+1,j)-aw*phi(k,i-1,j))*gij(k,k)
               sums = zero
               DO kl = 1,n
                 IF(k == kl)CYCLE
                 sums = sums - gij(k,kl)*((phi(kl,i+1,j)-two*phi(kl,i,j)+phi(kl,i-1,j))/dx2 &
                             + (phi(kl,i,j+1)-two*phi(kl,i,j)+phi(kl,i,j-1))/dy2)
               ENDDO
               by(j) = by(j) + sums
            ENDDO
            dy(ny) = one; by(ny) = phi(k,i,ny)

            CALL TRI(ny,ay,dy,cy,by,sy)

            phi(k,i,:) = sy(:)

          ENDDO  ! End Column-wise sweep
          
! Compute residual
          sumk = zero
          DO i = 2,nx-1
            DO j = 2,ny-1
              sums = zero
              DO kl = 1,n
                sums = sums + gij(k,kl)*((phi(kl,i+1,j)-two*phi(kl,i,j)+phi(kl,i-1,j))/dx2 &
                             + (phi(kl,i,j+1)-two*phi(kl,i,j)+phi(kl,i,j-1))/dy2)
              ENDDO
              sumk = sumk + sums**2
            ENDDO
          ENDDO
          resk = SQRT(MAX(zero,sumk))

          IF(resk < tol_in)EXIT

        ENDDO  ! Inner iteration loop

      ENDDO  ! Loop over dependent variables

! Compute residual
        DO k = 1,n
          sumk = zero
          DO i = 2,nx-1
            DO j = 2,ny-1
              sums = zero
              DO kl = 1,n
                sums = sums + gij(k,kl)*((phi(kl,i+1,j)-two*phi(kl,i,j)+phi(kl,i-1,j))/dx2 &
                             + (phi(kl,i,j+1)-two*phi(kl,i,j)+phi(kl,i,j-1))/dy2)
              ENDDO
              sumk = sumk + sums**2
            ENDDO
          ENDDO
          res(k) = SQRT(MAX(zero,sumk))
        ENDDO

! Monitor Outer Iteration residual
      converge = .true.
      DO k = 1,n
        IF(res(k) > tol_out) converge = .false.
      ENDDO
      WRITE(30,12) it_out, (res(k), k = 1,n)
      WRITE(*,12) it_out, (res(k), k = 1,n)
      IF(converge) EXIT
        
    ENDDO   ! Outer Iteration Loop

    time_end = etime(ta)

    print *, time_end-time_start

! Output Printout
    WRITE(*,*) ""
    WRITE(20,*) 'VARIABLES = "X", "Y", "Phi1", "Phi2", "Phi3"'
    WRITE(20,*) 'ZONE I=', nx, ', J=', ny, ', DATAPACKING=POINT'
    
    DO j = 1,ny
      DO i = 1,nx
        WRITE(20,11) x(i),y(j),(phi(k,i,j),k=1,n)
      ENDDO
    ENDDO
    

 11 FORMAT(5(1x,E14.6))
 12 FORMAT(I6,3(1x,E14.6))

    CLOSE(unit=20)
    CLOSE(unit=30)

  END program main
    
!------------------------------------------------------------------

! Tridiagonal solver (uses Thomas Algorithm)

    SUBROUTINE TRI(n,a,d,c,b,x)

    INTEGER(4), INTENT(in) :: n
    INTEGER(4) :: i
    REAL(8) :: xmult
    REAL(8), DIMENSION(n) :: a,b,c,d,x

    DO i = 2,n
      xmult = a(i-1)/d(i-1)
      d(i) = d(i) - xmult*c(i-1)
      b(i) = b(i) - xmult*b(i-1)
    ENDDO

    x(n) = b(n)/d(n)
    DO i = n-1,1,-1
      x(i) = (b(i)-c(i)*x(i+1))/d(i)
    ENDDO

    END SUBROUTINE TRI

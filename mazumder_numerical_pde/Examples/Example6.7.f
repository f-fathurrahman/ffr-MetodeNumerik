! Start of Main Program
! Example 6.7

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: nc = 40, mc = 40, max_iter = 100000
    INTEGER(4) :: i,j,iter
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0, &
                               tol = 1.0D-6
    REAL(8), PARAMETER :: phi_right = 0.0D0, &
                             & phi_left = 0.0d0, &
                               phi_bot = 1.0d0, &
                               phi_top = 0.0d0, &
                               length = 1.0D0, height = 1.0d0, &
                               costh = 0.5d0
    REAL(8), DIMENSION(nc) :: xi1,x
    REAL(8), DIMENSION(nc) :: xi2,y
    REAL(8), DIMENSION(nc,mc) :: phi
    REAL(8) :: dxi1,d12,dxi2,d21,sinth,sumres,res,c2
    REAL(8), DIMENSION(nc,mc) :: ao,ae,aw,an,as,ane,ase,anw,asw,rhs
    
! Open output file
    OPEN(unit=10, file="Example6_7.out", status = "unknown")
    OPEN(unit=20, file="Example6_7.rsl", status = "unknown")

! Geometry
    dxi1 = length/nc; dxi2 = height/mc
    d12 = dxi1/dxi2; d21 = dxi2/dxi1
    sinth = sqrt(one-costh*costh)
    c2 = costh/2.0d0
    xi1(1) = 0.5d0*dxi1
    DO i = 2,nc
      xi1(i) = xi1(i-1)+dxi1
    ENDDO
    xi2(1) = 0.5d0*dxi2
    DO j = 2,mc
      xi2(j) = xi2(j-1)+dxi2
    ENDDO

! Fill in Link Coefficients (they do not change with iteration)
    ao(:,:) = zero; rhs(:,:) = zero
    ae(:,:) = zero; aw(:,:) = zero; an(:,:) = zero; as(:,:) = zero
    ane(:,:) = zero; anw(:,:) = zero; ase(:,:) = zero; asw(:,:) = zero
! Left Bottom
    i = 1; j = 1
    ao(i,j) = 4.0d0*(d21+d12)+c2
    an(i,j) = -(4.0d0*d12/3.0d0-c2)
    ae(i,j) = -(4.0d0*d21/3.0d0-c2)
    ane(i,j) = c2
    rhs(i,j) = 8.0d0*(d21*phi_left+d12*phi_bot)/3.0d0+costh*(phi_left+phi_bot)
! Right Bottom
    i = nc; j = 1
    ao(i,j) = 4.0d0*(d21+d12)+c2
    an(i,j) = -(4.0d0*d12/3.0d0-c2)
    aw(i,j) = -(4.0d0*d21/3.0d0-c2)
    anw(i,j) = c2
    rhs(i,j) = 8.0d0*(d21*phi_right+d12*phi_bot)/3.0d0+costh*(phi_right+phi_bot)
! Left Top
    i = 1; j = mc
    ao(i,j) = 4.0d0*(d21+d12)+c2
    as(i,j) = -(4.0d0*d12/3.0d0-c2)
    ae(i,j) = -(4.0d0*d21/3.0d0-c2)
    ase(i,j) = c2
    rhs(i,j) = 8.0d0*(d21*phi_left+d12*phi_top)/3.0d0+costh*(phi_left+phi_top)
! Right Top
    i = nc; j = mc
    ao(i,j) = 4.0d0*(d21+d12)+c2
    as(i,j) = -(4.0d0*d12/3.0d0-c2)
    aw(i,j) = -(4.0d0*d21/3.0d0-c2)
    asw(i,j) = c2
    rhs(i,j) = 8.0d0*(d21*phi_right+d12*phi_top)/3.0d0+costh*(phi_right+phi_top)
! Left Wall
    i = 1
    DO j = 2,mc-1
      ao(i,j) = 4.0d0*d21+2.0d0*d12
      an(i,j) = -d12
      as(i,j) = -d12
      ae(i,j) = -(4.0d0*d21/3.0d0)
      ane(i,j) = c2
      ase(i,j) = -c2
      rhs(i,j) = 8.0d0*d21*phi_left/3.0d0
    ENDDO
! Right Wall
    i = nc
    DO j = 2,mc-1
      ao(i,j) = 4.0d0*d21+2.0d0*d12
      an(i,j) = -d12
      as(i,j) = -d12
      aw(i,j) = -(4.0d0*d21/3.0d0)
      anw(i,j) = c2
      asw(i,j) = -c2
      rhs(i,j) = 8.0d0*d21*phi_right/3.0d0
    ENDDO
! Bottom Wall
    j = 1
    DO i = 2,nc-1
      ao(i,j) = 4.0d0*d12+2.0d0*d21
      ae(i,j) = -d21
      aw(i,j) = -d21
      an(i,j) = -(4.0d0*d12/3.0d0)
      ane(i,j) = c2
      anw(i,j) = -c2
      rhs(i,j) = 8.0d0*d12*phi_bot/3.0d0
    ENDDO
! Top Wall
    j = mc
    DO i = 2,nc-1
      ao(i,j) = 4.0d0*d12+2.0d0*d21
      ae(i,j) = -d21
      aw(i,j) = -d21
      as(i,j) = -(4.0d0*d12/3.0d0)
      ase(i,j) = c2
      asw(i,j) = -c2
      rhs(i,j) = 8.0d0*d12*phi_top/3.0d0
    ENDDO
! Interior Cells
    DO j = 2,mc-1
      DO i = 2,nc-1
        ao(i,j) = 2.0d0*(d12+d21)
        ae(i,j) = -d21
        aw(i,j) = -d21
        an(i,j) = -d12
        as(i,j) = -d12
        ane(i,j) = c2
        anw(i,j) = -c2
        ase(i,j) = -c2
        asw(i,j) = c2
        rhs(i,j) = zero
      ENDDO
    ENDDO

! Start of Gauss-Seidel Iteration Loop

    phi(:,:) = zero  ! Initial guess

    DO iter = 1,max_iter

! Update solution
      DO i = 1,nc
        DO j = 1,mc
          phi(i,j) = (rhs(i,j) &
                   - ae(i,j)*phi(i+1,j) - aw(i,j)*phi(i-1,j) &
                   - an(i,j)*phi(i,j+1) - as(i,j)*phi(i,j-1) &
                   - ane(i,j)*phi(i+1,j+1) - ase(i,j)*phi(i+1,j-1) &
                   - anw(i,j)*phi(i-1,j+1) - asw(i,j)*phi(i-1,j-1))/ao(i,j)
        ENDDO
      ENDDO

! Residual calculation
      sumres = zero
      DO i = 1,nc
        DO j = 1,mc
          sumres = sumres + (rhs(i,j) - ao(i,j)*phi(i,j) & 
                   - ae(i,j)*phi(i+1,j) - aw(i,j)*phi(i-1,j) &
                   - an(i,j)*phi(i,j+1) - as(i,j)*phi(i,j-1) &
                   - ane(i,j)*phi(i+1,j+1) - ase(i,j)*phi(i+1,j-1) &
                   - anw(i,j)*phi(i-1,j+1) - asw(i,j)*phi(i-1,j-1))**2
        ENDDO
      ENDDO

      res = SQRT(MAX(zero,sumres))

      WRITE(*,*) iter, res
      WRITE(20,*) iter, res

      IF(res < tol)EXIT

    ENDDO  ! Iteration loop

! Print results

    WRITE(10,*) 'VARIABLES = "X", "Y", "PHI"'
    WRITE(10,*) 'ZONE I=', nc+2, ', J=', mc+2, ', DATAPACKING=POINT'

!Bottom Row
    WRITE(10,10) zero,zero,(phi_bot+phi_left)/2.0d0
    DO i = 1,nc
      x(i) = xi1(i)
      WRITE(10,10) x(i),zero,phi_bot
    ENDDO
    WRITE(10,10) length,zero,(phi_bot+phi_right)/2.0d0

!Interior Rows
    DO j = 1,mc
      y(j) = xi2(j)*sinth
      WRITE(10,10) xi2(j)*costh, y(j), phi_left
      DO i = 1,nc
         x(i) = xi1(i) + xi2(j)*costh
         WRITE(10,10) x(i),y(j),phi(i,j)
      ENDDO
      WRITE(10,10) length+xi2(j)*costh, y(j), phi_right
    ENDDO

!Top Row
    WRITE(10,10) height*costh,height*sinth,(phi_left+phi_top)/2.0d0
    DO i = 1,nc
      x(i) = xi1(i) + height*costh
      WRITE(10,10) x(i),height*sinth,phi_top
    ENDDO
    WRITE(10,10) height*costh+length,height*sinth,(phi_right+phi_top)/2.0d0


 10 format(3(1x,F13.6))

    CLOSE(unit=10)

    END Program main

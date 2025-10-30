! Coupled Solution of 3 PDEs 
! Example 8.11

  program main

    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 3
    INTEGER(4), PARAMETER :: nx = 101, ny = 101
    REAL(8), PARAMETER :: l_x = 1.0d0, l_y = 1.0d0
    INTEGER(4), PARAMETER :: max_out=100000
    INTEGER(4) :: i,j,it_out,k,kl
    REAL(8), PARAMETER :: one=1.0d0,half=0.5d0,zero=0.0d0,two=2.0d0
    REAL(8), PARAMETER :: tol_out = 1.0D-6
    REAL(8), DIMENSION(nx) :: x
    REAL(8), DIMENSION(nx,n,n) :: ax,bx,cx
    REAL(8), DIMENSION(nx,n) :: rx,sx
    REAL(8), DIMENSION(ny) :: y
    REAL(8), DIMENSION(ny,n,n) :: ay,by,cy
    REAL(8), DIMENSION(ny,n) :: ry,sy
    REAL(8), DIMENSION(n,n) :: gij
    REAL(8), DIMENSION(n) :: res
    REAL(8) :: phi(n,nx,ny)
    REAL(8) :: sums,sumk,resk,d_x,d_y,dx2,dy2,ao,aw,ae,as,an
    REAL(4) :: ta(2),time_start,time_end
    LOGICAL :: converge

! Open Output file
    OPEN(unit=20,file="Example8_11_coupled.out",status="unknown")
    OPEN(unit=30,file="Example8_11_coupled.rsl",status="unknown")

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

! Start of Iteration Loop For Coupled Solution of Phi
    DO it_out = 1,max_out

! Row-wise sweep

! Interior Rows
          DO j = 2,ny-1
            
            ax(:,:,:)=zero;bx(:,:,:)=zero;cx(:,:,:)=zero;rx(:,:)=zero;sx(:,:)=zero
            DO k = 1,n
              bx(1,k,k) = one
              rx(1,k) = phi(k,1,j)
            ENDDO
            DO i = 2,nx-1

               DO k = 1,n
                 sums = zero
                 DO kl = 1,n
                   bx(i,k,kl) = ao*gij(k,kl)
                   cx(i,k,kl) = ae*gij(k,kl)
                   ax(i,k,kl) = aw*gij(k,kl)
                   sums = sums -gij(k,kl)*(phi(kl,i,j+1)+phi(kl,i,j-1))/dy2
                 ENDDO
                 rx(i,k) = sums
               ENDDO

            ENDDO
            DO k = 1,n
              bx(nx,k,k) = one
              rx(nx,k) = phi(k,nx,j)
            ENDDO

            CALL BLKTRI(nx,n,n,nx,ax,bx,cx,sx,rx)

            DO i = 1,nx
              DO k = 1,n
                phi(k,i,j) = sx(i,k)
              ENDDO
            ENDDO

          ENDDO  ! End Row-wise sweep

! Column-wise sweep

! Interior Column
          DO i = 2,nx-1
            
            ay(:,:,:)=zero;by(:,:,:)=zero;cy(:,:,:)=zero;ry(:,:)=zero;sy(:,:)=zero
            DO k = 1,n
              by(1,k,k) = one
              ry(1,k) = phi(k,i,1)
            ENDDO
            DO j = 2,ny-1

               DO k = 1,n
                 sums = zero
                 DO kl = 1,n
                   by(j,k,kl) = ao*gij(k,kl)
                   cy(j,k,kl) = an*gij(k,kl)
                   ay(j,k,kl) = as*gij(k,kl)
                   sums = sums -gij(k,kl)*(phi(kl,i+1,j)+phi(kl,i-1,j))/dx2
                 ENDDO
                 ry(j,k) = sums
               ENDDO

            ENDDO
            DO k = 1,n
              by(ny,k,k) = one
              ry(ny,k) = phi(k,i,ny)
            ENDDO

            CALL BLKTRI(ny,n,n,ny,ay,by,cy,sy,ry)

            DO j = 1,ny
              DO k = 1,n
                phi(k,i,j) = sy(j,k)
              ENDDO
            ENDDO

          ENDDO  ! End Column-wise sweep

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
    

!-----------------------------------------------------------------------
      SUBROUTINE BLKTRI(NDIM,MDIM,M,N,A,B,C,X,R)

!     .......THIS PROGRAM SOLVES A TRI-BLOCK-DIAGONAL MATRIX.
!      .......NDIM = DIMENSION OF THE A,B AND C IN THE MAIN PROGRAM
!      .......MDIM = DIMENSION OF BLOCKS IN MAIN PROGRAM
!      .......M = SIZE OF BLOCKS
!      .......N = NUMBER OF BLOCKS
!      .......A,B,C = M*M MATRICES; A GOES FROM 2 TO N (sub-diagonal), B GOES
!      .......        FROM 1 TO N (diagonal) AND C GOES FROM 1 TO N-1 (super-diagonal)
!      .......X = RESULT VECTOR
!      .......R = RIGHT HAND SIDE VECTOR

      IMPLICIT NONE

      INTEGER(4) :: NDIM,MDIM,M,N,I,J,K,KK,MP,MPM1,JP1,MM1,JB, &
                        NM1,JR,MR,KR
      REAL(8) :: T,S
      REAL(8):: A(NDIM,MDIM,MDIM),B(NDIM,MDIM,MDIM),C(NDIM,MDIM,MDIM), &
      X(NDIM,MDIM),R(NDIM,MDIM)

!.....FORWARD ELIMINATION-BLOCKS
      DO I=1,N   ! loop 1
!.....GAUSS JORDAN FOR B(I)
        MM1 = M-1
        DO J=1,MM1 !loop 2
!.....NORMALIZE PIVOT ROW
          JP1=J+1
          DO K=JP1,M
            T=1./B(I,J,J)
            B(I,J,K) = T*B(I,J,K)
          ENDDO
          DO K=1,M
            IF (I == N) CYCLE
            C(I,J,K) = T*C(I,J,K)
          ENDDO
          R(I,J) = T*R(I,J)
!.....LOWER ELIMINATION IN B
          DO K=JP1,M
            DO KK=JP1,M
              B(I,K,KK) = B(I,K,KK)-B(I,K,J)*B(I,J,KK)
            ENDDO
            DO KK=1,M
              IF (I == N) CYCLE
              C(I,K,KK) = C(I,K,KK)-B(I,K,J)*C(I,J,KK)
            ENDDO
            R(I,K) = R(I,K)-B(I,K,J)*R(I,J)
          ENDDO
        ENDDO !loop 2
!.....UPPER (GAUSS JORDAN) ELIMINATION IN B
        T= 1./B(I,M,M)
        DO KK = 1,M
          C(I,M,KK) = T*C(I,M,KK)
        ENDDO
        R(I,M) = T*R(I,M)
        MM1=M-1
        DO J=1,MM1
          MP = M-J+1
          MPM1=MP-1
          DO K=1,MPM1
            MR = MP-K
            DO KK=1,M
              IF (I == N)CYCLE
              C(I,MR,KK) = C(I,MR,KK)-B(I,MR,MP)*C(I,MP,KK)
            ENDDO
            R(I,MR) = R(I,MR)-B(I,MR,MP)*R(I,MP)
          ENDDO
        ENDDO
!.....B(I) IS NOW THE UNIT MATRIX
!.....ELIMINATE A(I+1)
        IF (I == N) CYCLE
        DO J=1,M
          DO K=1,M
            DO KK=1,M
              B(I+1,K,KK) = B(I+1,K,KK)-A(I+1,K,J)*C(I,J,KK)
            ENDDO
            R(I+1,K) = R(I+1,K)-A(I+1,K,J)*R(I,J)
          ENDDO
        ENDDO
      ENDDO !loop 1

!.....BACK SUBSTITUTION
      DO K=1,M
        X(N,K) = R(N,K)
      ENDDO
      NM1 = N-1
      DO J=1,NM1
        JB=N-J
        DO K=1,M
          KR = M-K+1
          S=0.
          DO KK=1,M
            S = C(JB,KR,KK)*X(JB+1,KK)+S
          ENDDO
          X(JB,KR) = R(JB,KR)-S
        ENDDO
      ENDDO

      END SUBROUTINE BLKTRI

! This is a program to solve the 2D Diffusion Equation in 
! a rectangular domain using the Unstructured FVM.
! This particular problem has a material property discontinuity
! in the center of the domain.
! Exercise 7.7
!***********************************************************************
! Definitions of parameters used in the geometry

 MODULE geometry

   IMPLICIT NONE
   SAVE

   INTEGER(4) :: n_c,n_v,n_f,n_f_b
   INTEGER(4), PARAMETER :: n = 4, m = 2   ! number of cells along each dir.
   REAL(8), PARAMETER :: length=2.0d0, height=1.0d0
   INTEGER(4), DIMENSION(:,:), POINTER :: lcf,lcv,lfv,lfc,snsign
   INTEGER(4), DIMENSION(:), POINTER :: bf_to_f,bface,bnode,f_to_bf
   REAL(8), DIMENSION(:), POINTER :: xc,yc,xv,yv,xf,yf,area,wfun,delta,vol,sumwt
   REAL(8), DIMENSION(:,:), POINTER :: sn,st,cwfun

! lcf = cell to face connectivity
! lcv = cell to vertex connectivity
! lfv = face to vertex connectivity
! lfc = face to cell connectivity
! area = area of face
! vol = volume of cell
! wfun = cell-to-face weighting function
! cwfun = cell-to-vertex weighting function
! xc = x coordinate of cell centroid
! yc = y coordinate of cell centroid
! xv = x coordinate of vextex
! yv = y coordinate of vertex
! xf = x coordinate of face
! yf = y coordinate of face
! sn = surface normal of face, has 3 components
! st = surface tangent of face, has 3 components
! bf_to_f = boundary face number to global face number 
! f_to_bf = global face number to boundary face number
! bface = flag stating if boundary face (0 for interior faces)
! bnode = flag stating if boundary node (0 for interior nodes)

 END MODULE geometry

!***********************************************************************
       program main

         IMPLICIT NONE

         CALL grid
         CALL solve

       END program main
!-------------------------------------------------------------------------
!      Calculation of grid and all geometry related data

       SUBROUTINE grid
       
       USE geometry

       IMPLICIT NONE
       INTEGER(4), PARAMETER :: info = 1
       INTEGER(4) :: i,j,ic,iv,ifc,ifc1,ifc2,ifc3,ifc4,ifc_b,jc,k, &
                         v1,v2,v3,v4,c1,c2
       REAL(8), PARAMETER :: one=1.0d0,two=2.0d0,four=4.0d0, &
                                  three=3.0d0,half=0.5d0,zero=0.0d0
       REAL(8) :: xx,yy,dx,dy,dx1,dy1,dx2,dy2,dlen,dlen1,dlen2, &
                       d1,d2,d3,d4,ndotl

       dx = length/n
       dy = height/m 
       
       n_c = 8  ! Total number of cells
       n_v = 15  ! Total number of vertices
       n_f = 22 ! Total number of faces

       ALLOCATE(xc(n_c),yc(n_c))
       ALLOCATE(xv(n_v),yv(n_v))

! Cell center coordinates. Do this manually by following the figure.
       xc(1) = dx/two; yc(1) = height-dy/two
       xc(2) = xc(1)+dx; yc(2) = yc(1)
       xc(3) = dx/two; yc(3) = dy/two
       xc(4) = xc(3)+dx; yc(4) = yc(3)
       xc(5) = xc(2)+dx; yc(5) = yc(2)
       xc(6) = xc(5)+dx; yc(6) = yc(5)
       xc(7) = xc(4)+dx; yc(7) = yc(4)
       xc(8) = xc(7)+dx; yc(8) = yc(7)

! Vertex coordinates. Do this manually by following the figure.
      xv(11) = zero; yv(11) = zero
      xv(10) = dx; yv(10) = zero
      xv(9) = xv(10)+dx; yv(9) = zero
      xv(8) = xv(9)+dx; yv(16) = zero
      xv(7) = xv(8)+dx; yv(7) = zero
      xv(12) = zero; yv(12) = dy
      xv(13) = xv(12)+dx; yv(13) = yv(12)
      xv(14) = xv(13)+dx; yv(14) = yv(13)
      xv(15) = xv(14)+dx; yv(15) = yv(14)
      xv(6) = xv(15)+dx; yv(6) = yv(15)
      xv(1) = zero; yv(1) = height
      xv(2) = xv(1)+dx; yv(2) = yv(1)
      xv(3) = xv(2)+dx; yv(3) = yv(2)
      xv(4) = xv(3)+dx; yv(4) = yv(3)
      xv(5) = xv(4)+dx; yv(5) = yv(4)
      

! Grid Connectivity
! Cell to face connectivity
       ALLOCATE(lcf(n_c,4))
       lcf(1,1)=13; lcf(1,2)=1; lcf(1,3)=12; lcf(1,4)=19
       lcf(2,1)=14; lcf(2,2)=2; lcf(2,3)=13; lcf(2,4)=20
       lcf(3,1)=16; lcf(3,2)=19; lcf(3,3)=11; lcf(3,4)=10
       lcf(4,1)=17; lcf(4,2)=20; lcf(4,3)=16; lcf(4,4)=9
       lcf(5,1)=15; lcf(5,2)=3; lcf(5,3)=14; lcf(5,4)=21
       lcf(6,1)=5; lcf(6,2)=4; lcf(6,3)=15; lcf(6,4)=22
       lcf(7,1)=18; lcf(7,2)=21; lcf(7,3)=17; lcf(7,4)=8
       lcf(8,1)=6; lcf(8,2)=22; lcf(8,3)=18; lcf(8,4)=7
       
! Cell to vertex connectivity
       ALLOCATE(lcv(n_c,4))
       lcv(1,1)=2; lcv(1,2)=1; lcv(1,3)=12; lcv(1,4)=13
       lcv(2,1)=3; lcv(2,2)=2; lcv(2,3)=13; lcv(2,4)=14
       lcv(3,1)=13; lcv(3,2)=12; lcv(3,3)=11; lcv(3,4)=10
       lcv(4,1)=14; lcv(4,2)=13; lcv(4,3)=10; lcv(4,4)=9
       lcv(5,1)=4; lcv(5,2)=3; lcv(5,3)=14; lcv(5,4)=15
       lcv(6,1)=5; lcv(6,2)=4; lcv(6,3)=15; lcv(6,4)=6
       lcv(7,1)=15; lcv(7,2)=14; lcv(7,3)=9; lcv(7,4)=8
       lcv(8,1)=6; lcv(8,2)=15; lcv(8,3)=8; lcv(8,4)=7

! Face to Vertex Connectivity
       ALLOCATE(lfv(n_f,2))
       lfv(1,1)=2; lfv(1,2)=1
       lfv(2,1)=3; lfv(2,2)=2
       lfv(3,1)=4; lfv(3,2)=3
       lfv(4,1)=5; lfv(4,2)=4
       lfv(5,1)=6; lfv(5,2)=5
       lfv(6,1)=7; lfv(6,2)=6
       lfv(7,1)=8; lfv(7,2)=7
       lfv(8,1)=9; lfv(8,2)=8
       lfv(9,1)=10; lfv(9,2)=9
       lfv(10,1)=11; lfv(10,2)=10
       lfv(11,1)=12; lfv(11,2)=11
       lfv(12,1)=1; lfv(12,2)=12
       lfv(13,1)=2; lfv(13,2)=13
       lfv(14,1)=3; lfv(14,2)=14
       lfv(15,1)=4; lfv(15,2)=15
       lfv(16,1)=13; lfv(16,2)=10
       lfv(17,1)=14; lfv(17,2)=9
       lfv(18,1)=15; lfv(18,2)=8
       lfv(19,1)=13; lfv(19,2)=12
       lfv(20,1)=14; lfv(20,2)=13
       lfv(21,1)=15; lfv(21,2)=14
       lfv(22,1)=6; lfv(22,2)=15

! Compute Face Center Coordinates
       ALLOCATE(xf(n_f),yf(n_f))
       DO ifc = 1,n_f
         v1 = lfv(ifc,1)
         v2 = lfv(ifc,2)
         xf(ifc) = half*(xv(v1)+xv(v2))
         yf(ifc) = half*(yv(v1)+yv(v2))
       ENDDO

! Boundary Face Indexing
       n_f_b = 12   ! Total number of boundary faces
       ALLOCATE(bf_to_f(n_f_b),f_to_bf(n_f))
       ALLOCATE(bface(n_f))  ! Binary indicating if boundary face
       bface(:) = 0     ! interior face
       f_to_bf(:) = -1
       bf_to_f(:) = -1
       bface(1) = 1; bface(2) = 1; bface(3) = 1; bface(4) = 1
       bface(5) = 1; bface(6) = 1; bface(7) = 1; bface(8) = 1
       bface(9) = 1; bface(10) = 1; bface(11) = 1; bface(12) = 1
       bf_to_f(1) = 1; f_to_bf(1) = 1
       bf_to_f(2) = 2; f_to_bf(2) = 2
       bf_to_f(3) = 3; f_to_bf(3) = 3
       bf_to_f(4) = 4; f_to_bf(4) = 4
       bf_to_f(5) = 5; f_to_bf(5) = 5
       bf_to_f(6) = 6; f_to_bf(6) = 6
       bf_to_f(7) = 7; f_to_bf(7) = 7
       bf_to_f(8) = 8; f_to_bf(8) = 8
       bf_to_f(9) = 9; f_to_bf(9) = 9
       bf_to_f(10) = 10; f_to_bf(10) = 10
       bf_to_f(11) = 11; f_to_bf(11) = 11
       bf_to_f(12) = 12; f_to_bf(12) = 12

! Face to Cell Connectivity
      ALLOCATE(lfc(n_f,2))
      lfc(:,:) = -999
      lfc(1,1) = 1; lfc(1,2) = 1
      lfc(2,1) = 2; lfc(2,2) = 2
      lfc(3,1) = 5; lfc(3,2) = 5
      lfc(4,1) = 6; lfc(4,2) = 6
      lfc(5,1) = 6; lfc(5,2) = 6
      lfc(6,1) = 8; lfc(6,2) = 8
      lfc(7,1) = 8; lfc(7,2) = 8
      lfc(8,1) = 7; lfc(8,2) = 7
      lfc(9,1) = 4; lfc(9,2) = 4
      lfc(10,1) = 3; lfc(10,2) = 3
      lfc(11,1) = 3; lfc(11,2) = 3
      lfc(12,1) = 1; lfc(12,2) = 1
      lfc(13,1) = 1; lfc(13,2) = 2
      lfc(14,1) = 2; lfc(14,2) = 5
      lfc(15,1) = 5; lfc(15,2) = 6
      lfc(16,1) = 3; lfc(16,2) = 4
      lfc(17,1) = 4; lfc(17,2) = 7
      lfc(18,1) = 7; lfc(18,2) = 8
      lfc(19,1) = 1; lfc(19,2) = 3
      lfc(20,1) = 2; lfc(20,2) = 4
      lfc(21,1) = 5; lfc(21,2) = 7
      lfc(22,1) = 6; lfc(22,2) = 8

! Area, unit surface normal, and unit surface tangent
! Surface Normal always pointed from cell1 to cell2
! At Boundary, surface normal is pointed outward
! Surface tangent is pointed from vertex1 to vertex2

      ALLOCATE(area(n_f))
      ALLOCATE(sn(n_f,2))
      ALLOCATE(st(n_f,2))

      DO ifc = 1,n_f
        v1 = lfv(ifc,1)
        v2 = lfv(ifc,2)
        dx1 = xv(v2) - xv(v1)
        dy1 = yv(v2) - yv(v1)
        area(ifc) = SQRT(dx1*dx1+dy1*dy1)
        sn(ifc,1) = dy1/area(ifc)
        sn(ifc,2) = -dx1/area(ifc)
        st(ifc,1) = dx1/area(ifc)
        st(ifc,2) = dy1/area(ifc)
      ENDDO

! Make sure surface normal points from cell 1 to cell 2
      DO ifc = 1,n_f
        c1 = lfc(ifc,1)
        c2 = lfc(ifc,2)
        IF(c1 == c2) THEN  ! Make external normals point out
          dx1 = xf(ifc)-xc(c1)
          dy1 = yf(ifc)-yc(c1)
          ndotl = dx1*sn(ifc,1)+dy1*sn(ifc,2)
          IF(ndotl < 0.0d0)THEN
            sn(ifc,:) = -sn(ifc,:)
            st(ifc,:) = -st(ifc,:)
          ENDIF
        ELSE  ! Point from cell 1 to cell 2
          dx1 = xc(c2)-xc(c1)
          dy1 = yc(c2)-yc(c1)
          ndotl = dx1*sn(ifc,1)+dy1*sn(ifc,2)
          IF(ndotl < 0.0d0)THEN
            sn(ifc,:) = -sn(ifc,:)
            st(ifc,:) = -st(ifc,:)
          ENDIF
        ENDIF
      ENDDO

! Calculate delta = n.l
      ALLOCATE(delta(n_f))
      DO ifc = 1,n_f
        c1 = lfc(ifc,1)
        c2 = lfc(ifc,2)
        IF(c1 == c2)THEN  ! Boundary faces
          dx1 = xf(ifc)-xc(c1)
          dy1 = yf(ifc)-yc(c1)
          delta(ifc) = ABS(dx1*sn(ifc,1)+dy1*sn(ifc,2))
        ELSE              ! Interior faces
          dx1 = xc(c2)-xc(c1)
          dy1 = yc(c2)-yc(c1)
          delta(ifc) = ABS(dx1*sn(ifc,1)+dy1*sn(ifc,2))
        ENDIF
      ENDDO

! Make surface normal point outward from cell's perspective
      ALLOCATE(snsign(n_c,4))
      DO ic = 1,n_c
        DO j = 1,4  ! Loop over faces of cell (assuming quads)
          ifc = lcf(ic,j)
          c1 = lfc(ifc,1)
          IF(ic == c1)THEN  ! Surface normal already out
            snsign(ic,j) = 1
          ELSE               ! Flip direction
            snsign(ic,j) = -1
          ENDIF
        ENDDO
      ENDDO

! Calculate cell-to-face interpolation weights 
! Note that the interpolation function has been computed such that
! face_value = wfun*cell_value(1) + (1-wfun)*cell_value(2)
          
      ALLOCATE(wfun(n_f))
      wfun(:) = zero
      DO ifc = 1,n_f
        c1 = lfc(ifc,1)
        c2 = lfc(ifc,2)
        dx1 = xf(ifc)-xc(c1)
        dy1 = yf(ifc)-yc(c1)
        dlen1 = SQRT(dx1*dx1+dy1*dy1)
        dx2 = xf(ifc)-xc(c2)
        dy2 = yf(ifc)-yc(c2)
        dlen2 = SQRT(dx2*dx2+dy2*dy2)
        wfun(ifc) = dlen2/(dlen1+dlen2)
      ENDDO
        
! Calculate cell-to-vertex interpolation weights 
! The interpolation weight is 1/distance
! Vextex value = SUM(cwfun*cell_value)

      ALLOCATE(cwfun(n_c,4),sumwt(n_v))
      cwfun(:,:) = zero
      sumwt(:) = zero
      DO ic = 1,n_c
        DO iv = 1,4 ! assuming 4 vertices
          v1 = lcv(ic,iv)
          d1 = SQRT((xv(v1)-xc(ic))**2 + (yv(v1)-yc(ic))**2)
          cwfun(ic,iv) = one/d1
          sumwt(v1) = sumwt(v1) + one/d1
        ENDDO        
      ENDDO

      DO ic=1,n_c
        DO iv=1,4
          v1 = lcv(ic,iv)
          cwfun(ic,iv) = cwfun(ic,iv)/sumwt(v1)
        END DO
      END DO

! Calculation of Cell Volume

      ALLOCATE(vol(n_c))
      vol(:) = zero
      DO ic = 1,n_c
        v1 = lcv(ic,1)
        v2 = lcv(ic,2)
        v3 = lcv(ic,3)
        v4 = lcv(ic,4)
        dx1 = xv(v1)-xv(v3)
        dy1 = yv(v1)-yv(v3)
        dx2 = xv(v2)-xv(v3)
        dy2 = yv(v2)-yv(v3)
        vol(ic) = half*(dx1*dy2 - dy1*dx2)
        dx1 = xv(v1)-xv(v4)
        dy1 = yv(v1)-yv(v4)
        dx2 = xv(v2)-xv(v4)
        dy2 = yv(v2)-yv(v4)
        vol(ic) = vol(ic) + half*(dx1*dy2 - dy1*dx2)
      ENDDO

! Tag Boundary nodes
      ALLOCATE(bnode(n_v))
      bnode(:) = 0
      DO ifc=1,n_f
        IF (bface(ifc) == 1) THEN
          DO j=1,2
            v1 = lfv(ifc,j)
            bnode(v1) = 1
          END DO
        END IF
      END DO

! Write out Grid Data

      IF(info > 0)THEN

        OPEN(unit = 10, file = "Exercise7_7_grid.dat", status = "unknown")

        WRITE(10,*) "Number of Cells =", n_c
        WRITE(10,*) "Number of Faces =", n_f
        WRITE(10,*) "Number of Vertices =", n_v
        WRITE(10,*) "Number of Boundary Faces =", n_f_b

        WRITE(10,*) "Cell Center Data"
        DO ic = 1,n_c
          WRITE(10,11) ic, xc(ic),yc(ic),vol(ic)
        ENDDO

        WRITE(10,*) "Face Center Data"
        DO ifc = 1,n_f
          WRITE(10,11) ifc, xf(ifc),yf(ifc),area(ifc)
        ENDDO
  
        WRITE(10,*) "Surface Normal Data"
        DO ifc = 1,n_f
          WRITE(10,13) ifc, xf(ifc),yf(ifc),sn(ifc,1),sn(ifc,2)
        ENDDO

        WRITE(10,*) "Delta"
        DO ifc = 1,n_f
          WRITE(10,11) ifc, xf(ifc),yf(ifc),delta(ifc)
        ENDDO

        WRITE(10,*) "Vertex Data"
        DO iv = 1,n_v
          WRITE(10,12) iv, xv(iv),yv(iv)
        ENDDO

        WRITE(10,*) "Cell to Face Connectivity"
        DO ic = 1,n_c
          WRITE(10,*) ic, (lcf(ic,ifc), ifc = 1,4)
        ENDDO

        WRITE(10,*) "Face to Cell Connectivity"
        DO ifc = 1,n_f
          WRITE(10,*) ifc, (lfc(ifc,ic), ic = 1,2)
        ENDDO

        CLOSE(unit = 10)

      ENDIF

 11   format(i5,3(F13.6))
 12   format(i5,2(F13.6))
 13   format(i5,4(F13.6))
      
     END SUBROUTINE grid

!...........................................................................
     SUBROUTINE solve

     USE geometry

     IMPLICIT NONE
     INTEGER(4) :: ic,ifc,iv,ifb,j,icn,iter,c1,c2,v1,v2
     INTEGER(4), PARAMETER :: max_iter = 100
     REAL(8), DIMENSION(:), POINTER :: ap,sc,phib,phi,phinode,weight,scskew,gamf,gam
     REAL(8), DIMENSION(:,:), POINTER :: anb
     REAL(8) :: sumf,sumr,res,dx1,dy1,tdotl,flux
     REAL(8) :: gam1,gam2,phian
     REAL(8), PARAMETER :: tol = 1.0D-8, phi_left = 1.0d0

     ALLOCATE(ap(n_c),sc(n_c),anb(n_c,4),phib(n_f_b),phi(n_c))
     ALLOCATE(phinode(n_v),weight(n_v))
     ALLOCATE(scskew(n_c),gamf(n_f),gam(n_c))
     phi(:) = 0.0d0


! Open Residual File
     OPEN(unit=30,file="Exercise7_7.rsl", status="unknown")

! Hardwire BC at Left Boundary
     phib(:) = 0.0d0
     DO ifb = 1,n_f_b
       ifc = bf_to_f(ifb)
       IF(xf(ifc) < 1.0D-3)THEN  ! Left Boundary
         phib(ifb) = phi_left
       ENDIF
     ENDDO

! Hardwire gamma values
     gam1 = 1.0d0; gam2 = 10.0d0
     gam(1:4) = gam1; gam(5:8) = gam2

! Compute gamma at faces
     DO ifc = 1,n_f
       c1 = lfc(ifc,1)
       c2 = lfc(ifc,2)
       gamf(ifc) = wfun(ifc)*gam(c1) + (1.0d0-wfun(ifc))*gam(c2) !Arithmatic mean
       IF(ifc == 14 .or. ifc == 17) gamf(ifc) = gam(c1)*gam(c2)/gamf(ifc) !Harmonic mean
     ENDDO

! Calculate Link Coefficients
     DO ic = 1, n_c
        ap(ic) = 0.0d0
        sc(ic) = 0.0d0
        DO j = 1,4
          ifc = lcf(ic,j)
          IF(bface(ifc) == 0)THEN ! Interior
            ap(ic) = ap(ic) + area(ifc)*gamf(ifc)/delta(ifc)
            anb(ic,j) = -area(ifc)*gamf(ifc)/delta(ifc)
          ELSE
            ifb = f_to_bf(ifc)
            anb(ic,j) = 0.0d0
            IF(ifb == 5 .or. ifb == 6 .or. ifb == 11 .or. ifb == 12)THEN ! Hardwire
              ap(ic) = ap(ic) + area(ifc)*gamf(ifc)/delta(ifc)
              sc(ic) = sc(ic) + phib(ifb)*area(ifc)*gamf(ifc)/delta(ifc)
            ENDIF
          ENDIF
        ENDDO
      ENDDO

! Gauss-Seidel Iteration
      DO iter = 1,max_iter

! Compute Vertex Values
        phinode(:) = 0.0d0; weight(:) = 0.0d0
        DO ic = 1,n_c
          DO j=1,4
             iv = lcv(ic,j)
             IF(bnode(iv) == 1)CYCLE
             weight(iv) = weight(iv) + cwfun(ic,j)*gam(ic)
             phinode(iv) = phinode(iv) + phi(ic)*cwfun(ic,j)*gam(ic)
          END DO
        END DO
        DO iv = 1,n_v
          IF(bnode(iv) == 1)CYCLE
          phinode(iv) = phinode(iv)/weight(iv)
        END DO
! Hardwire left boundary nodes
       DO iv = 1,n_v
         IF(bnode(iv) == 1)THEN
           IF(xv(iv) < 1.0d-3)phinode(iv) = phi_left
         ENDIF
       ENDDO
     
! Compute tangential flux (skew) source

      DO ic = 1,n_c
        scskew(ic) = 0.0d0
        sumf = 0.0d0
        DO j = 1,4
          ifc = lcf(ic,j)
          IF(bface(ifc) == 1)CYCLE
          c1 = lfc(ifc,1)
          c2 = lfc(ifc,2)
          v1 = lfv(ifc,1)
          v2 = lfv(ifc,2)
          dx1 = xc(c2)-xc(c1)
          dy1 = yc(c2)-yc(c1)
          tdotl = st(ifc,1)*dx1+st(ifc,2)*dy1
          sumf = sumf + tdotl*(phinode(v2)-phinode(v1))*snsign(ic,j)*gamf(ifc)/delta(ifc)
        ENDDO
        scskew(ic) = sumf
      ENDDO
          
! Update solution
        DO ic = 1,n_c
          sumf = 0.0d0
          DO j = 1,4
            ifc = lcf(ic,j)
            IF(snsign(ic,j) == 1)THEN
              icn = lfc(ifc,2)
            ELSE
              icn = lfc(ifc,1)
            ENDIF
            sumf = sumf + anb(ic,j)*phi(icn)
          ENDDO
          phi(ic) = (sc(ic)+scskew(ic)-sumf)/ap(ic)
        ENDDO

! Update boundary values for Top and Bottom boundaries
        DO ifb = 1,n_f_b
          IF(ifb == 5 .or. ifb == 6 .or. ifb == 11 .or. ifb == 12)CYCLE ! Hardwire
          ifc = bf_to_f(ifb)
          ic = lfc(ifc,1)
          phib(ifb) = phi(ic)  ! Zero flux (first order treatment)
        ENDDO

! Compute Residual
        sumr = 0.0d0
        DO ic = 1,n_c
          sumf = 0.0d0
          DO j = 1,4
            ifc = lcf(ic,j)
            IF(snsign(ic,j) == 1)THEN
              icn = lfc(ifc,2)
            ELSE
              icn = lfc(ifc,1)
            ENDIF
            sumf = sumf + anb(ic,j)*phi(icn)
          ENDDO
          sumr = sumr + (sc(ic)+scskew(ic)-ap(ic)*phi(ic)-sumf)**2
        ENDDO
        res = SQRT(MAX(0.0d0,sumr))

        print *, iter,res
        WRITE(30,*) iter, res

        IF(res < tol) EXIT
      
      ENDDO ! Iteration loop

! Compute Vertex Values
      phinode(:) = 0.0d0; weight(:) = 0.0d0
      DO ic = 1,n_c
        DO j=1,4
           iv = lcv(ic,j)
           IF(bnode(iv) == 1)CYCLE
           weight(iv) = weight(iv) + cwfun(ic,j)*gam(ic)
           phinode(iv) = phinode(iv) + phi(ic)*cwfun(ic,j)*gam(ic)
        END DO
      END DO
      DO iv = 1,n_v
        IF(bnode(iv) == 1)CYCLE
        phinode(iv) = phinode(iv)/weight(iv)
      ENDDO
! Hardwire Left boundary nodes
     DO iv = 1,n_v
       IF(bnode(iv) == 1)THEN
         IF(xv(iv) < 1.0d-3)phinode(iv) = phi_left
       ENDIF
     ENDDO

! Write Solution for TECPLOT in FE format
      OPEN(unit=20,file="Exercise7_7_contour.dat", status="unknown")
      WRITE(20,*) 'VARIABLES = "X", "Y", "PHI"'
      WRITE(20,*) 'ZONE N=', n_v,', E=',n_c,', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
      DO iv = 1,n_v
        WRITE(20,*) xv(iv),yv(iv),phinode(iv)
      END DO
      DO ic = 1,n_c
        WRITE(20,*) (lcv(ic,j), j=1,4)
      END DO

! Print out Cell Center Values
      WRITE(20,*) "Cell Center Values"
      DO ic = 1,n_c
        IF(xc(ic) <= 1.0d0)THEN
          phian = -gam2*xc(ic)/(gam1+gam2)+1.0d0
        ELSE
          phian = -gam1*xc(ic)/(gam1+gam2) + 2.0d0*gam1/(gam1+gam2)
        ENDIF
        WRITE(20,*) ic, phi(ic), phian
      ENDDO
! Print out Fluxes
      WRITE(20,*) "Fluxes"
      sumf = 0.0d0
      DO ifb = 1,n_f_b
       IF(ifb == 5 .or. ifb == 6 .or. ifb == 11 .or. ifb == 12)THEN !Hardwire
         ifc = bf_to_f(ifb)
         ic = lfc(ifc,1)
         flux = gamf(ifc)*(phib(ifb)-phi(ic))/delta(ifc)
         WRITE(20,*) ifb,flux
         sumf = sumf + flux*area(ifc)
       ENDIF
     ENDDO
     WRITE(20,*) "Imbalance=", sumf


 11   format(i5,3(f13.6))

     END SUBROUTINE solve

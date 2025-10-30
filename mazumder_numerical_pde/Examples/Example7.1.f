! Example 7.1
! Unstructured FV solution of the Poisson Equation in a 2D rhombus.

!***********************************************************************
! Definitions of parameters used in the geometry

 MODULE geometry

   IMPLICIT NONE
   SAVE

   INTEGER(4) :: n_c,n_v,n_f,n_f_b
   INTEGER(4), PARAMETER :: n = 40   ! number of cells along each dir.
   REAL(8), PARAMETER :: len=1.0d0

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
       REAL(8) :: xx,yy,d,pi,theta,dx1,dy1,dx2,dy2,dlen,dlen1,dlen2, &
                       d1,d2,d3,d4,ndotl

       pi = four*atan(one)
       theta = pi/three

       d = len/n    ! size of each cell
       ic = 0
       DO j = 1,n
         DO i = 1,n
           ic = ic + 1
         ENDDO
       ENDDO
       iv = 0
       DO j = 1,n+1
         DO i = 1,n+1
           iv = iv + 1
         ENDDO
       ENDDO
       n_c = ic  ! Total number of cells
       n_v = iv  ! Total number of vertices

       ALLOCATE(xc(n_c),yc(n_c))
       ALLOCATE(xv(n_v),yv(n_v))

! Cell center coordinates
       ic = 0
       DO j = 1,n
         DO i = 1,n
           ic = ic + 1
           xc(ic) = (i-0.5d0)*d + (j-0.5d0)*d*cos(theta)
           yc(ic) = (j-0.5d0)*d*sin(theta)
         ENDDO
       ENDDO

! Vertex coordinates
       iv = 0
       DO j = 1,n+1
         DO i = 1,n+1
           iv = iv + 1
           xv(iv) = (i-1)*d + (j-1)*d*cos(theta)
           yv(iv) = (j-1)*d*sin(theta)
         ENDDO
       ENDDO

! Grid Connectivity
! Cell to face connectivity
       ALLOCATE(lcf(n_c,4))
       ic = 0
       DO j = 1,n
         DO i = 1,n
           ic = ic + 1
           lcf(ic,1) = ic
           lcf(ic,2) = n_c + n + ic + j 
           lcf(ic,3) = ic+n
           lcf(ic,4) = n_c + n + ic + j - 1
         ENDDO
       ENDDO
       n_f = lcf(n_c,2) ! Total number of faces
         
! Cell to vertex connectivity
       ALLOCATE(lcv(n_c,4))
       ic = 0
       DO j = 1,n
         DO i = 1,n
           ic = ic + 1
           lcv(ic,1) = ic + j - 1
           lcv(ic,2) = ic + j
           lcv(ic,3) = ic + n + j + 1
           lcv(ic,4) = ic + n + j
         ENDDO
       ENDDO

! Face to Vertex Connectivity: Vertices counter-clockwise
       ALLOCATE(lfv(n_f,2))
       DO ic = 1,n_c
         ifc1 = lcf(ic,1)
         lfv(ifc1,1) = lcv(ic,1)   ! southern face vertices
         lfv(ifc1,2) = lcv(ic,2)   ! southern face vertices
         ifc3 = lcf(ic,3)
         lfv(ifc3,1) = lcv(ic,3)   ! northern face vertices
         lfv(ifc3,2) = lcv(ic,4)   ! northern face vertices
         ifc2 = lcf(ic,2)
         lfv(ifc2,1) = lcv(ic,2)   ! eastern face vertices
         lfv(ifc2,2) = lcv(ic,3)   ! eastern face vertices
         ifc4 = lcf(ic,4)
         lfv(ifc4,1) = lcv(ic,4)   ! western face vertices
         lfv(ifc4,2) = lcv(ic,1)   ! western face vertices
       ENDDO

! Compute Face Center Coordinates
       ALLOCATE(xf(n_f),yf(n_f))
       DO ifc = 1,n_f
         v1 = lfv(ifc,1)
         v2 = lfv(ifc,2)
         xf(ifc) = half*(xv(v1)+xv(v2))
         yf(ifc) = half*(yv(v1)+yv(v2))
       ENDDO
       
! Boundary Face Indexing
       n_f_b = 2*(n+n)   ! Total number of boundary faces
       ALLOCATE(bf_to_f(n_f_b),f_to_bf(n_f))
       ALLOCATE(bface(n_f))  ! Binary indicating if boundary face
       bface(:) = 0     ! interior face
       f_to_bf(:) = 0
       ifc_b = 0
       DO i = 1,n
         ic = i 
         ifc_b = ifc_b + 1
         bf_to_f(ifc_b) = lcf(ic,1)
         bface(lcf(ic,1)) = 1 ! boundary face
         f_to_bf(lcf(ic,1)) = ifc_b
       ENDDO
       DO i = 1,n
         ic = i + (n-1)*n
         ifc_b = ifc_b + 1
         bf_to_f(ifc_b) = lcf(ic,3)
         bface(lcf(ic,3)) = 1
         f_to_bf(lcf(ic,3)) = ifc_b
       ENDDO
       DO j = 1,n
         ic = 1 + (j-1)*n
         ifc_b = ifc_b + 1
         bf_to_f(ifc_b) = lcf(ic,4)
         bface(lcf(ic,4)) = 1
         f_to_bf(lcf(ic,4)) = ifc_b
       ENDDO 
       DO j = 1,n
         ic = n + (j-1)*n
         ifc_b = ifc_b + 1
         bf_to_f(ifc_b) = lcf(ic,2)
         bface(lcf(ic,2)) = 1
         f_to_bf(lcf(ic,2)) = ifc_b
       ENDDO 

! Face to Cell Connectivity
      ALLOCATE(lfc(n_f,2))
      lfc(:,:) = -999
      DO ic = 1,n_c
        DO k = 1,4
          ifc1 = lcf(ic,k)
          IF(lfc(ifc1,1) /= -999)CYCLE
          lfc(ifc1,1) = ic
          IF(bface(ifc1) == 1)THEN
            lfc(ifc1,2) = ic
          ELSE
            IF(k == 1)THEN
              jc = ic - n
            ELSEIF(k == 3)THEN
              jc = ic + n
            ELSEIF(k == 2)THEN
              jc = ic + 1
            ELSEIF(k == 4)THEN
              jc = ic - 1
            ENDIF
            lfc(ifc1,2) = jc
          ENDIF
        ENDDO
      ENDDO

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

        OPEN(unit = 10, file = "grid.dat", status = "unknown")

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
     INTEGER(4), PARAMETER :: max_iter = 100000
     REAL(8), DIMENSION(:), POINTER :: ap,sc,phib,phi,phinode,weight,scskew
     REAL(8), DIMENSION(:,:), POINTER :: anb
     REAL(8) :: sumf,sumr,res,dx1,dy1,tdotl
     REAL(8), PARAMETER :: tol = 1.0D-8, phi_bot = 1.0d0

     ALLOCATE(ap(n_c),sc(n_c),anb(n_c,4),phib(n_f_b),phi(n_c))
     ALLOCATE(phinode(n_v),weight(n_v))
     ALLOCATE(scskew(n_c))
     phi(:) = 0.0d0


! Open Residual File
     OPEN(unit=30,file="Example7_1.rsl", status="unknown")

! Hardwire BC at bottom Wall
     phib(:) = 0.0d0
     DO ifb = 1,n_f_b
       ifc = bf_to_f(ifb)
       IF(yf(ifc) < 1.0D-3)THEN  ! Bottom Wall
         phib(ifb) = phi_bot
       ENDIF
     ENDDO

! Calculate Link Coefficients
     DO ic = 1, n_c
        ap(ic) = 0.0d0
        sc(ic) = 0.0d0
        DO j = 1,4
          ifc = lcf(ic,j)
          IF(bface(ifc) == 0)THEN ! Interior
            ap(ic) = ap(ic) + area(ifc)/delta(ifc)
            anb(ic,j) = -area(ifc)/delta(ifc)
          ELSE
            ifb = f_to_bf(ifc)
            ap(ic) = ap(ic) + area(ifc)/delta(ifc)
            anb(ic,j) = 0.0d0
            sc(ic) = sc(ic) + phib(ifb)*area(ifc)/delta(ifc)
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
             weight(iv) = weight(iv) + cwfun(ic,j)
             phinode(iv) = phinode(iv) + phi(ic)*cwfun(ic,j)
          END DO
        END DO
! Hardwire lower boundary nodes
       DO iv = 1,n_v
         IF(bnode(iv) == 1)THEN
           IF(yv(iv) < 1.0d-3)phinode(iv) = phi_bot
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
          sumf = sumf + tdotl*(phinode(v2)-phinode(v1))*snsign(ic,j)/delta(ifc)
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
           weight(iv) = weight(iv) + cwfun(ic,j)
           phinode(iv) = phinode(iv) + phi(ic)*cwfun(ic,j)
        END DO
      END DO
! Hardwire lower boundary nodes
     DO iv = 1,n_v
       IF(bnode(iv) == 1)THEN
         IF(yv(iv) < 1.0d-3)phinode(iv) = phi_bot
       ENDIF
     ENDDO

! Write Solution for TECPLOT in FE format
      OPEN(unit=20,file="Example7_1.dat", status="unknown")
      WRITE(20,*) 'VARIABLES = "X", "Y", "PHI"'
      WRITE(20,*) 'ZONE N=', n_v,', E=',n_c,', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
      DO iv = 1,n_v
        WRITE(20,*) xv(iv),yv(iv),phinode(iv)
      END DO
      DO ic = 1,n_c
        WRITE(20,*) (lcv(ic,j), j=1,4)
      END DO

 11   format(i5,3(f13.6))

     END SUBROUTINE solve

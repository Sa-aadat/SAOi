! SAOi: 
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch 
! SAOi: Quazi-Cauch routines coded by Pierre Duysinx, slightly modified 
! SAOi: for SAOi 
! SAOi: 


      SUBROUTINE QCAUCHY (ITER,IQC,X,S,D,G,G0,Y,N,
     &                    iuser,luser,cuser,ruser)
C     ******************
C
C Historique de la routine (creation,modification,correction)
C +--------------------------------------------------------------------+
C !Programmeur! Commentaires                          ! Date   !Version!
C +-----------!---------------------------------------!--------!-------+
C ! DUYSINX   ! Creation                              !02-04-01!       !
C +--------------------------------------------------------------------+
C 
C  QUASI CAUCHY UPDATING
C  *********************
C  
C  IQC = Type of update                                         (I)
C        1 updating diagonal matrix D 
C        2 updating root of diagonalmatrix D^1/2
C
C  D(N) Diagonal Hessian matrix at stage  [k]                   (I)
C                               at stage [k+1] (or +)           (O)
C  X(N) Design Point [k+1]                                      (I)
C  S(N) = X(N)[k+1] - X(N)[k] Step in design space              (I)
C  Y(N) = GRAD F(X[k])                                          (I)
C         GRAD F(X[k+1]) - GRAD F(X[k])                         (O)
C         Difference of the function gradients.
C  G(N) = GRAD F(X[k+1]) function gradient at stage k+1         (I)
C  G0(N) = GRAD F(X[k])  function gradient at stage k           (I)
C
C----------------------------------------------------------------------+
C
      IMPLICIT         NONE
      include          'ctrl.h'
      INTEGER          ITER,IQC,N
      DOUBLE PRECISION X(N),S(N),D(N),G(N),G0(N),Y(N)
C
      INTEGER I
      DOUBLE PRECISION ZERO,ONE,EPS,TOL,SDS,SY,DII,ZX,ZG,ALPHA
C
      include          'ctrl_get.inc'
      ZERO  = 0.D 00
      ONE   = 1.D 00
      EPS   = 1.D-06
      TOL   = 1.D-10
      ALPHA = 3.3D+00
C
C ********************************************************
C * ITERATION 1 : INITIALISATION OF THE DIAGONALE MATRIX *
C ********************************************************
C   Initial value of Dii = 1
C   or for structures with a CONLIN curvature
C
      IF(ITER.EQ.1) THEN
C
        IF (IQC.EQ.0) THEN
          ZX = ZERO
          ZG = ZERO
          DO I = 1,N
            ZX = ZX + X(I)*X(I)
            ZG = ZG + G(I)*G(I)
          ENDDO
          ZX = DSQRT(ZX)
          ZG = DSQRT(ZG)
          DII = ALPHA *ZX/ZG
c         D(I) = ALPHA*DABS(G(I)/X(I))
c         DII = one
          DO I = 1,N
            D(I) = DII
          ENDDO
        ELSE
          DO 1 I=1,N
!           D(I) = ALPHA*DABS(G(I)/X(I))      ! ag2009 (original)
!           D(I) = one                        ! ag2009 (new 1)
            D(I) = 2.d0*DABS(G(I)/X(I))       ! ag2009 (new 2 = t2:r)
 1        CONTINUE
        ENDIF
C
      ELSE ! end ITER.EQ.1
C
C *************************************************
C * ITERATION k>1 : UPDATING THE DIAGONAL HESSIAN *
C *************************************************
C
C S  : Design variable step
C Y  : Step of function gradients 
C D  : Current diagonal estimation of second order terms
C
        DO 2 I=1,N
 2      Y(I)  = G(I) - G0(I)
c       write(6,*) 'y= ',(y(i),i=1,n)
c       write(2,*) 'y= ',(y(i),i=1,n)
C
        SDS = ZERO
        SY = ZERO
        DO 3 I=1,N
          SDS   = SDS + D(I)*S(I)*S(I)
          SY    = SY  + Y(I)*S(I) 
 3      CONTINUE
c       write(6,*) 'sds= ',sds,' sy= ',sy
c       write(2,*) 'sds= ',sds,' sy= ',sy
C
C SY < 0 : non convex function !
C Updating is impossible
C Reduce trust region: D:= D*2.
C
        IF (SY.LT.ZERO) THEN
          !write(6,*) 'Warning: qcauchy: s*y <0 non convex function !',sy
          write(2,*) 'Warning: qcauchy: s*y <0 non convex function !',
     &                sy
          IF (IQC.EQ.0) THEN
            ZX = ZERO
            ZG = ZERO
            DO I = 1,N
              ZX = ZX + X(I)*X(I)
              ZG = ZG + G(I)*G(I)
            ENDDO
            ZX = DSQRT(ZX)
            ZG = DSQRT(ZG)
            DII = ALPHA *ZX/ZG
c           DII = one
            DO I = 1,N
              D(I) = DII
              D(I) = 2.d0*DABS(G(I)/X(I))       ! ag2009 (new 2 = t2:r)
            ENDDO
          ELSE
            DO 4 I=1,N
C             D(I) = D(I) * 1.4
C             D(I) = EPS
              D(I) = ALPHA*DABS(G(I)/X(I))
c             D(I) = one
              D(I) = 2.d0*DABS(G(I)/X(I))       ! ag2009 (new 2 = t2:r)
 4          CONTINUE
          ENDIF
          RETURN
        ENDIF
C
C Quasi Cauchy Update
C
        IF (IQC.EQ.0) THEN
          DII = SY/DMAX1(SDS,TOL)
          DO 5 I = 1,N
            D(I) = DII
 5        CONTINUE
        ELSE IF (IQC.EQ.1) THEN
          CALL CAUCHY1 (S,D,N,SDS,SY)
C         ============
        ELSE
          CALL CAUCHY2 (S,D,N,TOL,SDS,SY)
C         ============
        ENDIF
C
C ... MAX CURVATURE
C
        DO 11 I = 1,N
          D(I) = DMIN1(D(I),10.*DABS(G(I)/X(I)))
          D(I) = DMAX1(D(I),0.1*DABS(G(I)/X(I)))
          D(I) = 2.d0*DABS(G(I)/X(I))       ! ag2009 (new 2 = t2:r)
 11     CONTINUE
C 
      ENDIF  ! end ITER.GT.1
      RETURN
C
      END
C----------------------------------------------------------------------+
      SUBROUTINE CAUCHY0 (ITER,X,S,D,G,G0,Y,N,
     &                    iuser,luser,cuser,ruser)
C     ******************
C
C Historique de la routine (creation,modification,correction)
C +--------------------------------------------------------------------+
C !Programmeur! Commentaires                          ! Date   !Version!
C +-----------!---------------------------------------!--------!-------+
C ! DUYSINX   ! Creation                              !02-04-01!       !
C +--------------------------------------------------------------------+
C 
C  QUASI CAUCHY UPDATING : OREN LUENBERGER UPDATE
C  **********************************************
C  
C
C  D(N) Diagonal Hessian matrix at stage  [k]                   (I)
C                               at stage [k+1] (or +)           (O)
C  X(N) Design Point [k+1]                                      (I)
C  S(N) = X(N)[k+1] - X(N)[k] Step in design space              (I)
C  Y(N) = GRAD F(X[k])                                          (I)
C         GRAD F(X[k+1]) - GRAD F(X[k])                         (O)
C         Difference of the function gradients.
C  G(N) = GRAD F(X[k+1]) function gradient at stage k+1         (I)
C  G0(N) = GRAD F(X[k])  function gradient at stage k           (I)
C
C----------------------------------------------------------------------+
C
      IMPLICIT         NONE
      include          'ctrl.h'
      INTEGER          ITER,N
      DOUBLE PRECISION X(N),S(N),D(N),G(N),G0(N),Y(N)
C
      INTEGER I
      DOUBLE PRECISION ZERO,ONE,EPS,TOL,SDS,SY,DII,ZX,ZG
      include          'ctrl_get.inc'
C
      ZERO  = 0.D 00
      ONE    = 1.D 00
      EPS   = 1.D-06
      TOL   = 1.D-10
C
      ZX = ZERO
      ZG = ZERO
      DO 1 I = 1,N
          ZX = ZX + X(I)*X(I)
          ZG = ZG + G(I)*G(I)
 1    CONTINUE
        ZX = DSQRT(ZX)
        ZG = DSQRT(ZG)
C
C ********************************************************
C * ITERATION 1 : INITIALISATION OF THE DIAGONALE MATRIX *
C ********************************************************
C   Initial value of Dii = 1
C   or for structures with a CONLIN curvature
C
      IF(ITER.EQ.1) THEN
C
          DII = 3.3D+00*ZX/ZG
          DO 2 I = 1,N
            D(I) = DII
 2      CONTINUE
      ENDIF
C
C *************************************************
C * ITERATION k>1 : UPDATING THE DIAGONAL HESSIAN *
C *************************************************
C
C S  : Design variable step
C Y  : Step of function gradients 
C D  : Current diagonal estimation of second order terms
C
      DO 11 I=1,N
 11   Y(I)  = G(I) - G0(I)
C
      SDS = ZERO
      SY = ZERO
      DO 12 I=1,N
        SDS   = SDS + D(I)*S(I)*S(I)
        SY    = SY  + Y(I)*S(I) 
 12   CONTINUE
c      write(6,*) 'sds= ',sds,' sy= ',sy
c      write(2,*) 'sds= ',sds,' sy= ',sy
C
C SY < 0 : non convex function !
C Updating is impossible
C Reduce trust region: D:= D*2.
C
      IF (SY.LT.ZERO) THEN
c        write(6,*) 'Warning: qcauchy: s*y <0 non convex function !',sy
c        write(2,*) 'Warning: qcauchy: s*y <0 non convex function !',sy
          DII = 3.3D+00*ZX/ZG
        DO 13 I=1,N
          D(I) = DII
 13     CONTINUE
        RETURN    
        ENDIF
C
C Quasi Cauchy Update Oren Luenberger
C
        DII = SY/DMAX1(SDS,TOL)
        DO 14 I = 1,N
          D(I) = DII
 14   CONTINUE
C
      RETURN
C
      END
C----------------------------------------------------------------------+
      subroutine cauchy1 (s,d,n,a,b)
c     ==========
C
C Historique de la routine (creation,modification,correction)
C +--------------------------------------------------------------------+
C !Programmeur! Commentaires                          ! Date   !Version!
C +-----------!---------------------------------------!--------!-------+
C ! DUYSINX   ! Creation                              !02-04-01!       !
C +--------------------------------------------------------------------+
C 
C  D(N) Diagonal Hessian matrix at stage  [k]                   (I)
C                               at new stage [k+1] (or +)       (O)
C  X(N) New Design Point [k+1]                                  (I)
C  S(N) = X(N)[k+1] - X(N)[k] Direction de progression          (I)
C  Y(N) = GRAD F(X[k+1]) - GRAD F(X[k])                         (I)
C  a   = S.D.S                                                  (I)
C  b   = S.Y                                                                                                      (I)
c
c ******************************************************************
c * Quasi Cauchy update of diagonal matrix D = diag{d_i}           *
c ******************************************************************
c
      implicit         none
      integer          n
      double precision s(*),d(*),a,b
c       
      integer i
      double precision zero,c
      data zero/0.0d+00/
c
      !write(6,*) 'Cauchy 1'
c
c ... Update D
c
c     D+ = D + Lambda
c        Lambda = {(b-a) over tr(E^2)} E
c       with E = diag{s_i^2}
c
c
      c = zero
        do 2 i=1,n
          c = c + s(i)**4
 2    continue
c
      c = (b-a)/c  
        do 3 i = 1,n
          d(i) = d(i) + c*s(i)*s(i)
3     continue    
c
c Is D definite positive ? May be not !
c
      !write(6,*) 'Cauchy 1 : di updated ',(d(i),i=1,n)
      write(2,*) 'Cauchy 1 : di updated ',(d(i),i=1,n)
c
      return
      end
C----------------------------------------------------------------------+
      subroutine cauchy2 (s,d,n,tol,a,b)
c     ==========
C
C Historique de la routine (creation,modification,correction)
C +--------------------------------------------------------------------+
C !Programmeur! Commentaires                          ! Date   !Version!
C +-----------!---------------------------------------!--------!-------+
C ! DUYSINX   ! Creation                              !02-04-01!       !
C +--------------------------------------------------------------------+
C 
C  D(N) Diagonal Hessian matrix at stage  [k]                   (I)
C                               at new stage [k+1] (or +)       (O)
C  X(N) New Design Point [k+1]                                  (I)
C  S(N) = X(N)[k+1] - X(N)[k] Direction de progression          (I)
C  Y(N) = GRAD F(X[k+1]) - GRAD F(X[k])                         (I)
C  A   = S.D.S                                                  (I)
C  B   = S.Y                                                                                                      (I)
c
c ******************************************************************
c * Quasi Cauchy update of diagonal matrix D = diag{d_i}           *
c ******************************************************************
c
      implicit         none
      double precision s(*),d(*),tol
      integer          n,i
c       
      integer imax
      double precision zmu,f,df,zero,one,smax,pole,a,b,zmin,zmax,aux,
     &                 eps
      data zero/0.0d+00/,one/1.0d+00/
c
      write(6,*) 'Cauchy 2'
c
c ... Case b = sT y = 0
c     mu* tends to infinity
c     d(i) goes to zero
c
c      write(6,*) 'Cauchy 2 : sy=0'
c
        if (dabs(b).lt.tol) then
          eps = dsqrt(tol)
          do 1 i = 1,n
            d(i) = eps
 1        continue
        return
        endif
c
c ... Largest pole
c
c
c      write(6,*) 'Cauchy 2 : Largest pole'
c
      smax = -one/tol
        imax = 0
        do 2 i = 1,n
          pole = -one/(s(i)*s(i))
          if (pole.ge.smax) then
            smax = pole
            imax = i
          endif
 2    continue
c
c ... Solve non-linear equation F(zmu) = b
c
c
c      write(6,*) 'Cauchy 2 : solve linear equation'
c
      zmu = zero
        call cauchy3 (s,d,n,tol,f,zmu)
c       ==============
        call cauchy4 (s,d,n,tol,df,zmu)
c     ============
      if (dabs(f).lt.tol) goto 5
c 
      if (f.lt.b) then 
          zmax = zero
          zmin = smax + tol
        else
          zmin = zero
          zmax = 10.*(b-f)/df
        endif
c
        call cauchy5 (s,d,n,tol,b,zmin,zmax,zmu,f,df)
c     ============
c
 5    continue
c
c ... Update D
c
c     D+ = D                 if b=a
c        = (I + mu*E)^-2 D   if b =/ a
c       with E = diag{s_i^2}
c
c
c      write(6,*) 'Cauchy 2 : update D'
c
        if (dabs(b-a).gt.tol) then
          do 10 i = 1,n
            aux =       one + zmu*s(i)*s(i)
            aux = one / (aux*aux)
            d(i) = aux * d(i)
10      continue          
        endif
c
C      write(6,*) 'Cauchy 2 : di updated ',(d(i),i=1,n)
c
      return
      end
C----------------------------------------------------------------------+
      subroutine cauchy3 (s,d,n,tol,f,zmu)
c     ==========
c
c ******************************************************************
c * Computation of                                                 *
c * F(zmu) = Sum_{i=1}^{n} d_i s_i^2 over (1+zmu s_i^2)^2          *
c ******************************************************************
c
      double precision s(*),d(*),tol,zmu,f
      double precision deno,zero,one,prec
      integer n,i
      data zero/0.0d+00/,one/1.0d+00/
c
c      write(6,*) 'Cauchy 3'
c
      prec = tol*tol
      f = zero
      do 1 i = 1,n
          if (dabs(s(i)).gt.tol) then
            deno = one + zmu*s(i)*s(i)
            deno = deno*deno
            if (dabs(deno).lt.prec) deno = dsign(prec,deno) 
            f = f + d(i)*s(i)*s(i)/deno
          endif
 1      continue
c
      return
      end
C----------------------------------------------------------------------+
      subroutine cauchy4 (s,d,n,tol,df,zmu)
c     ==========
c
c ******************************************************************
c * Computation of the derivative of                               *
c * F(zmu) = Sum_{i=1}^{n} d_i s_i^2 over (1+zmu s_i^2)^2          *
c * dF over dzmu = -2 Sum_{i=1}^{n} d_i s_i^4 over (1+zmu s_i^2)^3 *
c ******************************************************************
c
      double precision s(*),d(*),tol,zmu,df
      double precision deno,zero,one,two,aux,prec
      integer n,i
      data zero/0.0d+00/,one/1.0d+00/,two/2.0d+00/
c
c      write(6,*) 'Cauchy 4'
c
      prec = tol*tol
      df = zero
      do 1 i = 1,n
          if (dabs(s(i)).gt.tol) then
            aux = s(i)*s(i)
            deno = one + zmu*aux
            deno = deno*deno*deno
            if (dabs(deno).lt.prec) deno = dsign(prec,deno) 
            df = df + d(i)*aux*aux/deno
          endif
 1      continue
      df = -two*df
c
      return
      end
C----------------------------------------------------------------------+
      subroutine cauchy5 (s,d,n,tol,b,zmin,zmax,zmu,f,df)
c     ==========
c
c ******************************************************************
c * Solve non-linear equation F(zmu) = b  in [zmin,zmax]           *
c * Method: bissection                                             *
c ******************************************************************
c
      implicit         none
      double precision s(*),d(*),b,tol,zmin,zmax,zmu,f,df
      integer          n,i
c       
      double precision zero,one,two,z1,z2,z3,f1,df1,f2,df2,f3,df3
      data zero/0.0d+00/,one/1.0d+00/,two/2.0d+00/
c
c      write(6,*) 'Cauchy 5'
c
c ... Initialization
c ... Check admissibility of interval [bmin,bmax]
c
c      write(6,*) 'Cauchy 5 : find interval'
c
      z1 = zmin
        z2 = zmax
        call cauchy3 (s,d,n,tol,f1,z1)
c       =============
      f1 = f1-b
        call cauchy4 (s,d,n,tol,df1,z1)
c       =============
      if (dabs(f1).lt.tol) then
          zmu = z1
        f = f1
         df = df1
         return
        endif
c
      i = 0
 10   i = i + 1
        call cauchy3 (s,d,n,tol,f2,z2)
c       =============
      f2 = f2-b
        call cauchy4 (s,d,n,tol,df2,z2)
c       =============
        if (f2.gt.zero) then
           z1 = z2
           z2 = z2 + 2.0d+00*dabs(z2-z1)
           if (i.lt.10000) then
             goto 10
           else
             zmu = z2
           f = f2
             df = df2
             return
           endif
        endif
c
c ... Bissection
c
c      write(6,*) 'Cauchy 5 : bissection'
c
      i = 0
 20   i = i+1
        z3 = (z1+z2)/two
      call cauchy3 (s,d,n,tol,f3,z3)
c       =============
      f3 = f3-b
        call cauchy4 (s,d,n,tol,df3,z3)
c       =============
        if (dabs(f3).gt.tol) then
          if (f3.gt.zero) then
            z1 = z3
            f1 = f3
            df1 = df3
          else
            z2 = z3
            f2 = f3
            df2 = df3
          endif
          if (i.lt.10000) goto 20
        endif
c
        zmu = z3 
      f = f3
        df = df3
c       write(6,*) 'zmu',zmu
c
      return
      end
C----------------------------------------------------------------------+


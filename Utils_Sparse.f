! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: utilities for the sparse SAOi routines
! SAOi:

 
      subroutine store_iters (n,ni,f_h,f,v_h,viol,x_h,x_h2,x,c_h,c,
     &                        cplus_norm,cplus_norm_h,gf_h,gf,gc_h,gc,
     &                        xlam,xlam_h,ic,ig,z,z_h,
     &                        iuser,luser,cuser,ruser,ictrl,lctrl,
     &                        rctrl,cctrl)
      implicit         none
      include          'ctrl.h'
      integer          i,j,n,ni,ic,ig
      double precision f_h,f,v_h,viol,x_h(n),x_h2(n),x(n),c_h(ni),c(ni)
      double precision gf_h(*),gf(n),gc_h(*),gc(*),xlam(ni),xlam_h(ni)
      double precision cplus_norm,cplus_norm_h,z(n),z_h(n)
      include          'ctrl_get.inc'
!
      if (ic.eq.1) then
        cplus_norm_h=cplus_norm
        f_h=f
        v_h=viol
        call dcopy (n,x_h,1,x_h2,1)
        call dcopy (n,x,1,x_h,1)
        call dcopy (ni,c,1,c_h,1)
        call dcopy (ni,xlam,1,xlam_h,1)
        call dcopy (n,z,1,z_h,1)
!
        if (ig.eq.0) return
!
!        call dcopy (n,gf,1,gf_h,1)
!        call dcopy (n*ni,gc,1,gc_h,1)
!
      elseif (ic.eq.-1) then
        cplus_norm=cplus_norm_h
        f=f_h
        viol=v_h
        call dcopy (n,x_h,1,x,1)
        call dcopy (ni,c_h,1,c,1)
        call dcopy (ni,xlam_h,1,xlam,1)
        call dcopy (n,z_h,1,z,1)
!
        if (ig.eq.0) return
!
!        call dcopy (n,gf_h,1,gf,1)
!        call dcopy (n*ni,gc_h,1,gc,1)
!
      else
        stop ' undefined in store_it'
      endif
!
      return
      end subroutine store_iters
c     --------------------------------------------------------------
c                SIMPLE DRIVER FOR L-BFGS-B (version 2.4)
c     --------------------------------------------------------------
       subroutine drive_lbfgsb24f_l (nprimal,xprimal,ni,nq,f_a,c_a,
     &                          h_a,iact,nact,x,x_h,
     &                          acurv,bcurv,ccurv,acurvn,
     &                          bcurvn,ccurvn,x_l,x_u,
     &                          f,c,h,gf,gc,gh,ksubiter,
     &                          xkkt,flam,cstage0,ostring,xlam_h,
     &                          ns,ictrl,lctrl,rctrl,cctrl,message,
     &                          nnz,Acol,Aptr,yhi,outerloop,eqn,lin,
     &                          amult,bmult,iuser,luser,cuser,ruser)

!  This driver calls the l-bfgs-b code

      implicit         none
      include          'ctrl.h'
      integer          umax
      parameter        (umax = 8)
      character*60     task,csave
      character*90     ostring
      logical          eqn(*), lin(*)
      logical          lsave(4),kkt_local,cstage0
      integer          ni,nq,nprimal,ksubiter,n,m,iprints,nnz
      integer          i,j,nact,iact(nact),nbd(ni),iwa(3*ni)
      integer          isave(44),isparse,outerloop
      integer          ns,icnt,message
      integer          Acol(nnz),Aptr(ni+1),ifree,jni
      double precision gcvec(1),bvec(1)
      double precision flambda,factr,pgtol,yhi
      double precision l(ni),u(ni),g(ni),time1,time2
      double precision dsave(29),xprimal(nprimal),flam
      double precision wa(2*umax*ni + 4*ni + 11*umax*umax + 8*umax)
      double precision x_h(nprimal),x_l(nprimal),x_u(nprimal),x(ni)
      double precision acurv,bcurv(ni),ccurv(nq),xlam_h(ni)
      double precision acurvn(1),bcurvn(1,1)
      double precision ccurvn(*),xlamsml,xlambig
      double precision f,c(ni),h(nq),gf(nprimal),gc(nnz)
      double precision f_a,c_a(ni),h_a(nq),gh_a(nq,nprimal)
      double precision xkkt,gf_a(nprimal),gc_a(nnz)
      double precision gh(nq,nprimal)
      double precision timef1,timef2,timef3
      double precision amult,bmult(ni)
      data             kkt_local /.false./
      include          'ctrl_get.inc'

! construct the sparse matrices
!       call sparset (nprimal,ni,gf,gc,ns,Aptr,nrow,gcvec,bvec,bcurvn,
!      &              ictrl,lctrl,rctrl,cctrl,
!      &              iuser,luser,cuser,ruser)

!      write(*,*) ' Entering subproblem '

! no output whatsoever
      iprints = -10

! specify the tolerances in the stopping criteria
      factr  =  1.0d1
      pgtol  =  1.0d-10

! specify the dimension n and the number m of limited memory corrections stored
      n      = ni
      m      = 8
!       
      if (m.gt.umax ) stop ' m.gt.umax  in Develop/j02.f'

! set bounds on the dual variables
      do 10 i = 1,ni
        nbd(i) = 2
        if (eqn(i)) then
          l(i)   =  max(-biglam,x(i)-DualTrustRadius)
        else    
          l(i)   =  max(0.d0   ,x(i)-DualTrustRadius)
        endif
        u(i)   = min(biglam,x(i)+DualTrustRadius)
   10 continue


      !if (.not.strict_struct) then
      !  jni = 0
      !  do i=1,ni
      !    if (Aptr(i+1)-Aptr(i).eq.0) then
      !      u(i)=0.d0
      !      jni=jni+1
      !    endif
      !  enddo
      !endif

! start the iteration by initializing task and timer
      call timer(time1)
      task = 'START'
      ostring=' subproblem terminated with standard settings'
      message=0

! entry point of the subproblem optimization loop
  111 continue

! call the l-bfgs-b code
      call setulb24(n,m,x,l,u,nbd,flambda,g,factr,pgtol,wa,iwa,task,
     &              iprints,csave,lsave,isave,dsave)
!

      !if (flambda.lt.1.d5.and.outerloop.ge.57) then
      !  write(*,*) 
      !  stop
      !endif


      if (task(1:2) .eq. 'FG') then

! check the time limit
        call timer(time2)

        if (time2 - time1 .gt. tlimit) then
          task='STOP: CPU time exceeds limit.'
          ostring=' subproblem terminated on time limit'
          message=1

        else

! construct the falk dual
          call falk_l (nprimal,ni,xprimal,x,x_l,x_u,x_h,f,c,
     &               gf,gc,f_a,c_a,g,
     &               ksubiter,iact,nact,flambda,
     &               ns,nnz,Aptr,Acol,bvec,ictrl,lctrl,rctrl,
     &               cctrl,acurv,bcurv,acurvn,bcurvn,yhi,xlam_h,
     &               amult,bmult,iuser,luser,cuser,ruser)
!
          goto 111
        endif

      elseif (task(1:5) .eq. 'NEW_X') then
! the minimization routine has returned with a new iterate,and we have opted to (conditionally) continue the iterations

! terminate if the total number of f and g evaluations exceeds ksubmax
        if (isave(34) .ge. ksubmax) then
          task='STOP: allowed no. of subproblem evaluations exceeded'
          ostring=' subproblem terminated on number of samples of the du
     &al'
          message=2
        endif

! terminate if  |proj g|/(1 + |f|) < 1.0d-10.
        if (dsave(13) .le. 1.d-10*(1.0d0 + dabs(f))) then
          task='STOP: projected gradient is sufficiently small'
          ostring=' subproblem terminated with |proj g|/(1 + |f|) < 1.0d
     &-10'
          message=0
        endif
        goto 111

! terminate execution when task is neither FG nor NEW_X
      else

! end of the loop
      endif

! calculate the approximate KKT conditions on the subproblem level
      flam=-flambda
!
      return
      end
c======================= The end of driver1 ============================
       subroutine falk_l (nprimal,ni,xprimal,xdual,x_l,x_u,x_h,
     &                 f,c,gf,gc,f_a,c_a,g,
     &                 kounts,iact,nact,flambda,
     &                 ns,nnz,Aptr,Acol,bvec,ictrl,lctrl,rctrl,
     &                 cctrl,acurv,bcurv,acurvn,bcurvn,yhi,xdual_h,
     &                 amult,bmult,iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,nprimal,ni,kounts,nact,iact(nact)
      double precision acurv,bcurv(ni),acurvn(1),bcurvn(1,1)
      integer          nnz,Acol(nnz),Aptr(ni+1)
      integer          ns,icnt
      double precision gc(nnz),bvec(1),xdual_h(ni)
      double precision xdual(ni),x_l(nprimal),x_u(nprimal)
      double precision xprimal(nprimal),gf(nprimal)
!     double precision acurv,bcurv(ni),acurvn(nprimal)
!     double precision bcurvn(ni,nprimal)
      double precision c(ni),f,f_a,c_a(ni),yhi
      double precision temp1,temp2,beta,flambda,g(ni),x_h(nprimal)
      double precision amult,bmult(ni)
      include          'ctrl_get.inc'

!  calculate the Falk dual - sparse implementation
      if (ifalk.eq.0) then
        call falk_dq_sprse_l (nprimal,ni,xprimal,xdual,x_l,x_u,
     &                        x_h,f,c,gf,f_a,c_a,g,kounts,iact,nact,
     &                        flambda,ns,nnz,Aptr,Acol,gc,bvec,ictrl,
     &                        lctrl,rctrl,cctrl,acurv,bcurv,acurvn,
     &                        bcurvn,yhi,amult,bmult,
     &                        iuser,luser,cuser,ruser)
      elseif (ifalk.eq.1) then
        call falk_dq_sprse_lq (nprimal,ni,xprimal,xdual,x_l,x_u,
     &                        x_h,f,c,gf,f_a,c_a,g,kounts,iact,nact,
     &                        flambda,ns,nnz,Aptr,Acol,gc,bvec,ictrl,
     &                        lctrl,rctrl,cctrl,acurv,bcurv,acurvn,
     &                        bcurvn,yhi,xdual_h,amult,bmult,
     &                        iuser,luser,cuser,ruser)
      else
        stop ' ifalk not defined in falk_l'
      endif
!
      return
      end
!----------------------------------------------------------------------c
      subroutine falk_dq_sprse_l (nprimal,ni,xprimal,xdual,x_l,x_u,
     &                            x_h,f,c,gf,f_a,c_a,g,kounts,iact,nact,
     &                            flambda,ns,nnz,Aptr,Acol,gc,bvec,
     &                            ictrl,lctrl,rctrl,cctrl,acurv,bcurv,
     &                            acurvn,bcurvn,yhi,amult,bmult,
     &                            iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,nprimal,ni,kounts,nact,iact(nact)
      integer          ij,i1
      double precision acurv,bcurv(ni),acurvn(1),bcurvn(1,1)
      integer          nnz,Acol(nnz),Aptr(ni+1)
      integer          ns,icnt
      double precision gc(nnz),bvec(1)
      double precision xdual(ni),x_l(nprimal),x_u(nprimal),zero,one,two
      double precision xprimal(nprimal),gf(nprimal)
!     double precision acurv,bcurv(ni),acurvn(nprimal)
      double precision gcji,bcji,yhi
      double precision c(ni),f,f_a,c_a(ni)
      double precision temp1(nprimal),temp2(nprimal)
      double precision beta,flambda,g(ni)
      double precision x_h(nprimal),xdualj,y,sum1
      double precision amult,bmult(ni)
      data             zero /0.d0/,one /1.d0/,two /2.d0/
      include          'ctrl_get.inc'

!  construct the bounded dual, looping over ndual
      do i=1,nprimal
        temp1(i)=gf(i)
        if (approx_f.eq.1) then
          temp2(i)=max(atol1,acurv*amult)
        elseif (approx_f.eq.4) then
          temp2(i)=max(atol1,2.d0*amult/x_h(i)*dabs(gf(i)))
        elseif (approx_f.eq.14) then
          temp2(i)=max(atol1,-2.d0*amult/x_h(i)*gf(i))
        else
          stop ' undefined in falk_dq_sprse_l '
        endif
      enddo
!
      ij=0
      do j=1,ni
        xdualj=xdual(j)
        do i1=1,Aptr(j+1)-Aptr(j)
          ij=ij+1
          i=Acol(ij)
          gcji=gc(ij)
          if (xdualj.gt.0.d0) then
            if (approx_c.eq.1) then
              bcji=max(btol1,bcurv(j)*bmult(j))
            elseif (approx_c.eq.4) then
              bcji=max(btol1,2.d0*bmult(j)/x_h(i)*dabs(gcji))
            elseif (approx_c.eq.14) then
              bcji=max(btol1,-2.d0*bmult(j)/x_h(i)*gcji)
            else
              stop ' undefined in falk_dq_sprse_l '
            endif
            temp1(i) = temp1(i) + xdualj*gcji
            temp2(i) = temp2(i) + xdualj*bcji
          endif
        enddo
      enddo
!
      do i=1,nprimal
        beta = x_h(i)-temp1(i)/temp2(i)
        xprimal(i) = min(x_u(i),max(x_l(i),beta))
      enddo

!  evaluate the relaxation variables if needed
      if (relax) then
        sum1 = zero
        do k = 1,nact
          j = iact(k)
          sum1=sum1+xdual(j)
        enddo
        y = (sum1-pen1)/pen2       !   we desire biglam > pen1   !
        y = min(ymax,max(zero,y))
      else
        y = zero
      endif

!  have updated xprimal; now get function value to maximize the dual
      kounts = kounts + 1
      call fun_as (nprimal, xprimal, f_a, x_h, acurv, acurvn,f,
     &             gf,ictrl,lctrl,rctrl,cctrl,amult,
     &             iuser,luser,cuser,ruser)
      call conin_as(nprimal, ni, xprimal, c_a, iact, nact, x_h, bcurv,
     &              bcurvn,c,gc,nnz,Aptr,Acol,ictrl,lctrl,rctrl,cctrl,
     &              bmult,iuser,luser,cuser,ruser)
      flambda = -f_a
      do k = 1,nact
        j = iact(k)
        if (Aptr(j+1)-Aptr(j).ne.0) then
          flambda = flambda - xdual(j)*c_a(j)
        else
          ! do nothing
        endif
      enddo

!  compute gradient g
      do k = 1,nact
        j = iact(k)
        if (Aptr(j+1)-Aptr(j).ne.0) then
          g(j) = -c_a(j)
        else
          g(j) = 0.d0
        endif
      enddo

!  do some corrections for relaxation
      if (relax) then
        flambda = flambda - pen1*y - pen2*y**2/two
        do k = 1,nact
          j = iact(k)
          if (Aptr(j+1)-Aptr(j).ne.0) then
            flambda = flambda + xdual(j)*y
            g(j) = g(j) + y
          endif  
        enddo
        yhi = y
      endif
!
      return
      end
!----------------------------------------------------------------------c
      subroutine falk_dq_sprse_lq (nprimal,ni,xprimal,xdual,x_l,x_u,
     &                            x_h,f,c,gf,f_a,c_a,g,kounts,iact,nact,
     &                            flambda,ns,nnz,Aptr,Acol,gc,bvec,
     &                            ictrl,lctrl,rctrl,cctrl,acurv,bcurv,
     &                            acurvn,bcurvn,yhi,xdual_h,
     &                            amult,bmult,iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,nprimal,ni,kounts,nact,iact(nact)
      integer          ij,i1
      double precision acurv,bcurv(ni),acurvn(1),bcurvn(1,1)
      integer          nnz,Acol(nnz),Aptr(ni+1)
      integer          ns,icnt
      double precision gc(nnz),bvec(1)
      double precision xdual(ni),x_l(nprimal),x_u(nprimal),zero,one,two
      double precision xprimal(nprimal),gf(nprimal)
!     double precision acurv,bcurv(ni),acurvn(nprimal)
      double precision gcji,bcji,yhi
      double precision c(ni),f,f_a,c_a(ni),xdual_h(ni)
      double precision temp1(nprimal),temp2(nprimal)
      double precision beta,flambda,g(ni)
      double precision x_h(nprimal),xdualt,xdualb,y,sum1
      double precision amult,bmult(ni)
      data             zero /0.d0/,one /1.d0/,two /2.d0/
      include          'ctrl_get.inc'

!  construct the bounded dual, looping over ndual
      do i=1,nprimal
        temp1(i)=gf(i)
        if (approx_f.eq.1) then
          temp2(i)=max(atol1,acurv*amult)
        elseif (approx_f.eq.4) then
          temp2(i)=max(atol1,2.d0*amult/x_h(i)*dabs(gf(i)))
        elseif (approx_f.eq.14) then
          temp2(i)=max(atol1,-2.d0*amult/x_h(i)*gf(i))
        else
          stop ' undefined in falk_dq_sprse_lq '
        endif
      enddo
!
      ij=0
      do j=1,ni
        xdualt=xdual(j)
        xdualb=xdual_h(j)
        do i1=1,Aptr(j+1)-Aptr(j)
          ij=ij+1
          i=Acol(ij)
          gcji=gc(ij)
          if (xdualt.ne.0.d0) then
            if (approx_c.eq.1) then
              bcji=max(btol1,bcurv(j)*bmult(j))
            elseif (approx_c.eq.4) then
              bcji=max(btol1,2.d0*bmult(j)/x_h(i)*dabs(gcji))
            elseif (approx_c.eq.14) then
              bcji=max(btol1,-2.d0*bmult(j)/x_h(i)*gcji)
            else
              stop ' undefined in falk_dq_sprse_lq '
            endif
            temp1(i) = temp1(i) + xdualt*gcji
          endif
          if (xdualb.ne.0.d0) then
            if (approx_c.eq.1) then
              bcji=max(btol1,bcurv(j)*bmult(j))
            elseif (approx_c.eq.4) then
              bcji=max(btol1,2.d0*bmult(j)/x_h(i)*dabs(gcji))
            elseif (approx_c.eq.14) then
              bcji=max(btol1,-2.d0*bmult(j)/x_h(i)*gcji)
            else
              stop ' undefined in falk_dq_sprse_l '
            endif
            temp2(i) = temp2(i) + xdualb*bcji
          endif
        enddo
      enddo
!
      do i=1,nprimal
        beta = x_h(i)-temp1(i)/max(1.d-8,temp2(i))
        xprimal(i) = min(x_u(i),max(x_l(i),beta))
      enddo

!  evaluate the relaxation variables if needed
      if (relax) then
        sum1 = zero
        do k = 1,nact
          j = iact(k)
          sum1=sum1+xdual(j)
        enddo
        y = (sum1-pen1)/pen2       !   we desire biglam > pen1   !
        y = min(ymax,max(zero,y))
      else
        y = zero
      endif

!  have updated xprimal; now get function value to maximize the dual
      kounts = kounts + 1
      call fun_asq (nprimal,ni,xprimal,f_a,x_h,acurv,acurvn,f,xdual_h,
     &              gf,ictrl,lctrl,rctrl,cctrl,bcurv,nnz,Aptr,Acol,
     &              iact,nact,amult,iuser,luser,cuser,ruser)
      call conin_asq(nprimal,ni,xprimal,c_a,iact,nact,x_h,c,gc,nnz,Aptr,
     &               Acol,ictrl,lctrl,rctrl,cctrl,bmult,
     &               iuser,luser,cuser,ruser)
      flambda = -f_a
      do k = 1,nact
        j = iact(k)
        flambda = flambda - xdual(j)*c_a(j)
      enddo

!  compute gradient g
      do k = 1,nact
        j = iact(k)
        g(j) = -c_a(j)
      enddo

!  do some corrections for relaxation
      if (relax) then
        flambda = flambda - pen1*y - pen2*y**2/two
        do k = 1,nact
          j = iact(k)
          flambda = flambda + xdual(j)*y
          g(j) = g(j) + y
        enddo
        yhi = y
      endif

!
      return
      end
!----------------------------------------------------------------------c
      subroutine fun_as (n,x,f_a,x_h,acurv,acurvn,f,gf,
     &                   ictrl,lctrl,rctrl,cctrl,amult,
     &                   iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          n,i
      double precision x(n),delx(n),x_h(n),f,gf(n)
      double precision acurv,acurvn(1),delx2(n),temp,dp,dp1,f_a
      double precision amult
      include          'ctrl_get.inc'
!
      do i=1,n
        temp     = x(i)-x_h(i)
        delx(i)  = temp
        delx2(i) = temp**2
      end do
!
      !xax = ddot(n,delx2,1,acurvn,1)
      !dp  = ddot(n,delx,1,gf,1)
      dp=0.d0
      dp1=0.d0
      do i=1,n
        dp=dp+delx(i)*gf(i)
        if (approx_f.eq.1) then
          dp1=dp1+delx2(i)*acurv*amult
        elseif (approx_f.eq.4) then
          dp1=dp1+delx2(i)*max(atol1,2.d0*amult/x_h(i)*dabs(gf(i)))
        elseif (approx_f.eq.14) then
          dp1=dp1+delx2(i)*max(atol1,-2.d0*amult/x_h(i)*gf(i))
        else
          write(*,*) ' illegal in conin_as'
        endif
      end do
!
      f_a = f+dp+dp1/2.d0
!
      return
      end subroutine fun_as
!-----------------------------------------------------------------------
      subroutine conin_as(n,ni,x,c_a,iact,nact,x_h,bcurv,bcurvn,
     &                    c,gc,nnz,Aptr,Acol,ictrl,lctrl,rctrl,cctrl,
     &                    bmult,iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          nnz,n,ni,nact,iact(nact),ij,i,j,i1,j1
      double precision x(n),delx(n),c_a(ni),x_h(n),c(ni),gc(nnz)
      double precision bcurv(ni),bcurvn(1,1),delx2(n),dp,dp1,temp
      double precision bmult(ni)
      integer          Acol(nnz),Aptr(ni+1)
      include          'ctrl_get.inc'
!
      do i=1,n
        temp    = x(i)-x_h(i)
        delx(i) = temp
        delx2(i)= temp**2
      end do
!
      ij=0
      do j1=1,nact
        j=iact(j1)
        dp=0.d0
        dp1=0.d0
        do i1=1,Aptr(j+1)-Aptr(j)
          ij=ij+1
          i=Acol(ij)
          dp=dp+delx(i)*gc(ij)
          if (approx_c.eq.1) then
            dp1=dp1+delx2(i)*max(btol1,bcurv(j)*bmult(j))
          elseif (approx_c.eq.4) then
           dp1=dp1+delx2(i)*max(btol1,2.d0*bmult(j)/x_h(i)*dabs(gc(ij)))
           elseif (approx_c.eq.14) then
           dp1=dp1+delx2(i)*max(btol1,-2.d0*bmult(j)/x_h(i)*gc(ij))
         else
            write(*,*) ' illegal in conin_as'
          endif
        end do
!
        c_a(j)=c(j)+dp+dp1/2.d0
!
      end do
!
      return
      end subroutine conin_as
!----------------------------------------------------------------------c
      subroutine fun_asq (n,ni,x,f_a,x_h,acurv,acurvn,f,xdual_h,
     &                    gf,ictrl,lctrl,rctrl,cctrl,bcurv,nnz,Aptr,Acol
     &                    ,iact,nact,amult,iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          n,i,ni,nnz,j,ij,j1,i1,nact,iact(nact)
      double precision x(n),delx(n),x_h(n),f,gf(n),bcurv(ni)
      double precision acurv,acurvn(1),delx2(n),temp,dp,dp1,f_a
      double precision c(ni),gc(nnz),xdual_h(ni)
      double precision bcurvn(1,1)
      double precision amult
      integer          Acol(nnz),Aptr(ni+1)
      include          'ctrl_get.inc'
!
      do i=1,n
        temp     = x(i)-x_h(i)
        delx(i)  = temp
        delx2(i) = temp**2
      end do
!
      !xax = ddot(n,delx2,1,acurvn,1)
      !dp  = ddot(n,delx,1,gf,1)
      dp=0.d0
      dp1=0.d0
      do i=1,n
        dp=dp+delx(i)*gf(i)
        if (approx_f.eq.1) then
          dp1=dp1+delx2(i)*max(atol1,acurv*amult)
        elseif (approx_f.eq.4) then
          dp1=dp1+delx2(i)*max(atol1,2.d0*amult/x_h(i)*dabs(gf(i)))
        elseif (approx_f.eq.14) then
          dp1=dp1+delx2(i)*max(atol1,-2.d0*amult/x_h(i)*gf(i))
        else
          write(*,*) ' illegal in fun_asq'
        endif
      end do
!
!       ij=0
!       do j1=1,nact
!         j=iact(j1)
!         dp1=0.d0
!         do i1=1,Aptr(j+1)-Aptr(j)
!           ij=ij+1
!           i=Acol(ij)
!           if (approx_c.eq.1) then
!             dp1=dp1+delx2(i)*max(btol1,bcurv(j))*xdual_h(j)
!           elseif (approx_c.eq.4) then
!             dp1=dp1+delx2(i)*max(btol1,2.d0/x_h(i)*dabs(gc(ij)))
!      &             *xdual_h(j)
!           elseif (approx_c.eq.14) then
!             dp1=dp1+delx2(i)*max(btol1,-2.d0/x_h(i)*(gc(ij)))
!      &             *xdual_h(j)
!           else
!             write(*,*) ' illegal in fun_asq'
!           endif
!         end do
!       end do
!
      f_a = f+dp+dp1/2.d0
!
      return
      end subroutine fun_asq
!-----------------------------------------------------------------------
      subroutine conin_asq (n,ni,x,c_a,iact,nact,x_h,c,gc,nnz,Aptr,
     &                      Acol,ictrl,lctrl,rctrl,cctrl,bmult,
     &                      iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          nnz,n,ni,nact,iact(nact),ij,i,j,i1,j1
      double precision x(n),delx(n),c_a(ni),x_h(n),c(ni),gc(nnz)
      double precision bcurv(ni),bcurvn(1,1),delx2(n),dp,dp1,temp
      double precision bmult(ni)
      integer          Acol(nnz),Aptr(ni+1)
      include          'ctrl_get.inc'
!
      do i=1,n
        temp    = x(i)-x_h(i)
        delx(i) = temp
        delx2(i)= temp**2
      end do
!
      ij=0
      do j1=1,nact
        j=iact(j1)
        dp=0.d0
        dp1=0.d0
        do i1=1,Aptr(j+1)-Aptr(j)
          ij=ij+1
          i=Acol(ij)
          dp=dp+delx(i)*gc(ij)
!           if (approx_c.eq.1) then
!             dp1=dp1+delx2(i)*bcurv(j)
!           elseif (approx_c.eq.4) then
!             dp1=dp1+delx2(i)*max(btol1,2.d0/x_h(i)*dabs(gc(ij)))
!           elseif (approx_c.eq.14) then
!             dp1=dp1+delx2(i)*max(btol1,-2.d0/x_h(i)*(gc(ij)))
!           else
!             write(*,*) ' illegal in conin_as'
!           endif
        end do
!
        c_a(j)=c(j)+dp!+dp1/2.d0
!
      end do
!
      return
      end subroutine conin_asq
!----------------------------------------------------------------------c
       subroutine form_kkts (xkkt,n,ni,nq,x,iact,nact,x_l0,x_u0,
     &                       gf,gc,gh,xlam,xlamsml,xlambig,nnz,
     &                       Acol,Aptr,ifree,
     &                       ictrl,lctrl,rctrl,cctrl,z,
     &                       iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          ni,nnz,Acol(nnz),Aptr(ni+1)
      integer          i,j,ij,j1,n,nq,nact,iact(nact),ifree,i1
      double precision x(n),x_l0(n),x_u0(n),gf(n),gc(nnz),gh(nq,n),gcji
      double precision xkkt,kkt(n),xlam(*),xlamsml,xlambig,xlamj,z(n)
!
      include          'ctrl_get.inc'
!
      ifree   = 0
      xkkt    = 0.d0
      xlamsml = 1.d16
      xlambig =-1.d16
!
      do j=1,ni
        xlamj=xlam(j)
        if (xlamj.gt.xlambig) xlambig=xlamj
        if (xlamj.lt.xlamsml) xlamsml=xlamj
      enddo
!
      do i=1,n
        kkt(i)=gf(i)
      end do
!
      ij=0
      do j1=1,nact ! previously np
        j=iact(j1)
        xlamj=xlam(j)
        do i1=1,Aptr(j+1)-Aptr(j)
          ij=ij+1
          i=Acol(ij)
          kkt(i)=kkt(i)+xlamj*gc(ij)
        end do
      end do
!
      xkkt=0.d0
      if (subsolver.ge.25.and.subsolver.le.28) then
        do i=1,n
          kkt(i)=kkt(i)-z(i)
          xkkt=xkkt+kkt(i)**2
          ifree=ifree+1
        enddo
      else  
        do i=1,n
          if     (x(i).le.x_l0(i)+0.d-6) then
!
          elseif (x(i).ge.x_u0(i)-0.d-6) then
!
          else
            xkkt=xkkt+kkt(i)**2
            ifree=ifree+1
          endif
        enddo
      endif
!
      xkkt=dsqrt(xkkt)
!
      return
      end subroutine form_kkts
!----------------------------------------------------------------------c
      subroutine diaHess_sq1s (n,ni,ne,acurv,bcurv,ccurv,x,x_h,
     &                        gf,gc,gh,f,c,h,f_h,c_h,h_h,
     &                        ictrl,lctrl,rctrl,cctrl,outerloop,
     &                        Acol,Aptr,nnz,
     &                        iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,outerloop,i1,ij,nnz
      integer          Acol(nnz),Aptr(ni+1)
      double precision acurv,bcurv(ni),ccurv(ne),x(n),x_h(n),delx(n)
      double precision gf(n),gc(nnz),gh(ne,n),dp,zero,one,delxnorm2
      double precision f,c(ni),h(ne),f_h,c_h(ni),h_h(ne),gcji
      data             zero /0.d0/
      data             one  /1.d0/
      include          'ctrl_get.inc'
!
! the classical spherical quadratic approximation
!
      if (outerloop.eq.0) then    ! initialize the approximation

        if (approx_f.eq.1) then
          acurv = one
        elseif (approx_f.eq.4) then               ! T2:R
          !  do nothing; formed in drive_galQPs
        elseif (approx_f.eq.14) then              ! T2:R
          !  do nothing; formed in drive_galQPs
        elseif (approx_f.ge.100) then             ! exact non-sparse
          !  do nothing; formed in drive_galQPs
        else
          write(*,*) 'approx_f: ',approx_f,'  ',cname1
          stop ' approx_f unknown in diaHess_sq1s'
        endif
!
        if (approx_c.eq.1) then
          do j=1,ni
            bcurv(j) = one ! btol1
          end do
        elseif (approx_c.eq.4) then               ! T2:R
          !  do nothing; formed in drive_galQPs
        elseif (approx_c.eq.14) then              ! T2:R
          !  do nothing; formed in drive_galQPs
        elseif (approx_c.ge.100) then             ! exact non-sparse
          !  do nothing; formed in drive_galQPs
        else
          write(*,*) 'approx_c: ',approx_c,'  ',cname1
          stop ' approx_c unknown in diaHess_sq1s'
        endif
!
      else     ! the iteration number > 0
!
        delxnorm2=0.d0
        do i=1,n
          delx(i)=x_h(i)-x(i)
          delxnorm2=delxnorm2+delx(i)**2
        end do
        if (approx_f.eq.1) then
          dp=0.d0
          do i=1,n
            dp=dp+gf(i)*delx(i)
          end do
          acurv=2.d0*(f_h-f-dp)/max(1.d-10,delxnorm2)
          acurv=max(atol1,acurv)
        elseif (approx_f.eq.4) then               ! T2:R
          !  do nothing; formed in drive_galQPs
        elseif (approx_f.eq.14) then              ! T2:R
          !  do nothing; formed in drive_galQPs
        elseif (approx_f.ge.100) then             ! exact non-sparse
          !  do nothing; formed in drive_galQPs
        else
          write(*,*) 'approx_f: ',approx_f,'  ',cname1
          stop ' approx_f unknown in diaHess_sq1s'
        endif
!
        if (approx_c.eq.1) then
          ij=0
          do j=1,ni
            dp=0.d0
            do i1=1,Aptr(j+1)-Aptr(j)
              ij=ij+1
              i=Acol(ij)
              gcji=gc(ij)
              dp=dp+gcji*delx(i)
            enddo
            bcurv(j)=2.d0*(c_h(j)-c(j)-dp)/max(1.d-10,delxnorm2)
            bcurv(j)=max(btol1,bcurv(j))
          enddo
        elseif (approx_c.eq.4) then               ! T2:R
          !  do nothing; formed in drive_galQPs
        elseif (approx_c.eq.14) then              ! T2:R
          !  do nothing; formed in drive_galQPs
        elseif (approx_c.ge.100) then             ! exact non-sparse
          !  do nothing; formed in drive_galQPs
        else
          write(*,*) 'approx_c: ',approx_c,'  ',cname1
          stop ' approx_c unknown in diaHess_sq1s'
        endif
      endif
!
      return
      end subroutine diaHess_sq1s
!----------------------------------------------------------------------c
      subroutine filter_ises (n,ni,accept,f,f_h,f_a,v,v_h,c,c_h,c_a,
     &                       kloop,acurv,acurvn,bcurv,
     &                       bcurvn,eps_locala,eps_localb,
     &                       filter_list,infilter,init,rho,
     &                       ictrl,lctrl,rctrl,cctrl,lloop,amult,bmult,
     &                       iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          kloop,infilter,i,n,ni,il,lloop
      double precision c(ni),c_h(ni),c_a(ni),filter_list(2,nmax)
      double precision acurv,acurvn(n),bcurv(ni),bcurvn(1,1)
      double precision rho,beta,gamma,delta_f
      double precision sigma,fact_red,delta_q
      double precision test1left,test2left,v,f,v_h,f_h
      double precision test1right,test2right,f_a,eps_locala,eps_localb
      double precision amult,bmult(ni)
      logical          accept,acceptable,update_filt,init,iprntl
!
      data             gamma                /1.00d-3/     !  /1.00d-5/
      data             sigma                /1.00d-2/     !  /1.00d-1/
      data             fact_red             /2.00d00/     !  /2.00d00/
      data             iprntl               /.false./     !  /.false./
!
      include          'ctrl_get.inc'

      beta=1.d0-gamma
!
      if (kloop.eq.1.and.init) then  ! initialize the filter and check
        do i=1,outermax
          filter_list(1,i)=0.d0
          filter_list(2,i)=0.d0
       enddo
        if (sigma.lt.gamma) stop ' sigma < gamma in filter_ise'
        filter_list(1,1)=filter_hi
        filter_list(2,1)=-big
        infilter=1
        init=.false.
      endif

!  determine if acceptabile to the filter
      test1left  = v
      test2left  = f + gamma*v
      do i=1,infilter
        test1right = beta*filter_list(1,i)
        test2right =      filter_list(2,i)
        if (test1left.le.test1right.or.test2left.le.test2right) then
          ! do nothing
        else
          acceptable=.false.
          accept=.false.
          if (classic_trust) then
            rho=rho/fact_red
            if (rho.le.xtol*10.d0) goto 1000
          else
            stop 'not possible in filter '
          call conserve2 (n,ni,accept,f,f_a,c,c_a,eps_locala,
     &                   eps_localb,acurv,bcurv,
     &                   rho,ictrl,lctrl,rctrl,cctrl,lloop,amult,
     &                   bmult,iuser,luser,cuser,ruser)
          endif
          update_filt=.false.
          if (iprntl) write(*,*)' case b: reducing rho;',
     &                            ' dominated by filter',rho,f,v
          goto 1000
        endif
      enddo

!  determine acceptability of the new optimal candidate point
!   x^k+1 = x^k + d, which has not yet been accepted or rejected,
!   w.r.t. the current point x^k <--> (h^k, f^k).
!   x^k has not yet been considered for inclusion in the filter;
!   this will only happen lower down in this subroutine.
!   (alternatively, test if h <= h^k here.)
      test1right = beta*v_h
      test2right = f_h
      if (test1left.le.test1right.or.test2left.le.test2right) then
        ! do nothing
      else
        acceptable=.false.
        accept=.false.
        if (classic_trust) then
          rho=rho/fact_red
        else
          stop 'not possible in filter '
!           call conserve(n,ni,accept,f,f_a,c,c_a,eps_locala,
!      &                    eps_localb,acurv,acurvn,bcurv,bcurvn,
!      &                    rho,ictrl,lctrl,rctrl,cctrl,lloop,
!      &                    iuser,luser,cuser,ruser)
        endif
        update_filt=.false.
                    if (iprntl) write(*,*)' case c: reducing rho;',
     &                           ' dominated by (h^k, f^k)',rho,f,v
        goto 1000
      endif

!  check for sufficient improvement
      delta_q = f_h-f_a
      delta_f = f_h-f
!
      if (delta_f.lt.sigma*delta_q.and.delta_q.gt.0.d0) then
        acceptable=.true.
        accept=.false.
        if (classic_trust) then
          rho=rho/fact_red
        else
          stop 'not possible in filter '
!           call conserve(n,ni,accept,f,f_a,c,c_a,eps_locala,
!      &                    eps_localb,acurv,acurvn,bcurv,bcurvn,
!      &                    rho,ictrl,lctrl,rctrl,cctrl,lloop,
!      &                    iuser,luser,cuser,ruser)
        endif
        update_filt=.false.
        if (iprntl) write(*,*)' case d: reducing rho;',
     &                          ' inadequate improvement',rho,f,v
        goto 1000
      endif
!
      if (delta_q.lt.0.d0) then
        acceptable=.true.
        accept=.true.
        update_filt=.true.
        if (iprntl) write(*,*)' case e: updating the filter;',
     &                          ' h-type iteration',rho,f,v
        goto 1000
      else
        acceptable=.true.
        accept=.true.
        update_filt=.false.
        if (iprntl) write(*,*)' case f: f-type iteration',rho,f,v
        goto 1000
      endif

! return, after checking rho and possibly updating the filter
 1000 continue
      rho=max(rho,rho_l_min)
      rho=min(rho,rho_l_max)

! add the last starting point (v_h,f_h) to the filter
      if (update_filt) then
        call update_filter_ise (filter_list,infilter,v_h,f_h,beta,gamma,
     &                          ictrl,lctrl,rctrl,cctrl,
     &                          iuser,luser,cuser,ruser)
      endif
!
      return
      end
!----------------------------------------------------------------------!
      subroutine conserve2 (n,ni,accept,f,f_a,c,c_a,eps_locala,
     &                      eps_localb,acurv,bcurv,rho,
     &                      ictrl,lctrl,rctrl,cctrl,lloop,amult,bmult,
     &                      iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,n,ni,lloop
      double precision bcurv(ni),fact_red,fact_inc,factx
      double precision rho,c(ni),c_a(ni),acurv,acurvn(n),zero,one,two
      double precision f,f_a,eps_locala,eps_localb,cbig
      double precision delxnorm2,delx(n),dp,fbig
      double precision amult,bmult(ni)
      logical          accept
!
      data             fbig                 /1.00d03/
      data             cbig                 /1.00d03/
      data             zero                 /0.00d00/
      data             one                  /1.00d00/
      data             two                  /2.00d00/
      data             fact_inc             /2.00d00/
      data             fact_red             /2.00d00/
!
      include          'ctrl_get.inc'
!
      if (iaggressive.eq.0) then
        factx=one
      elseif (iaggressive.eq.1) then
        factx=one*dble(lloop+1)
      elseif (iaggressive.eq.2) then
        factx=one*dble(lloop+1)**2
      else
        factx=one ! illegal; reset
      endif
!
      accept=.true.
      if (f_a.lt.f-eps_locala) then
        accept=.false.
        !acurv=min(fbig,max(acurv,atol2)*fact_inc*factx)
        amult = amult*fact_inc*factx
        !write(*,*) lloop,acurvn(1)
      endif
!
      if (.not.unconstrained) then
        do j=1,ni
          if (c_a(j).lt.c(j)-eps_localb.and.c(j).gt.-1.d-4) then
            accept=.false.
            !bcurv(j)=min(cbig,max(bcurv(j),btol2)*fact_inc*factx)
            !write(*,*) lloop,bcurvn(j,1)
            bmult(j) = bmult(j)*fact_inc*factx
          endif
        enddo
      endif
!
      if (.not.classic_conserv.and..not.accept) then
        rho=rho/fact_red
        rho=max(rho,rho_l_min)
        rho=min(rho,rho_l_max)
      endif
!
      return
      end subroutine conserve2
!
      subroutine FinDiffGradHessSparse (n,x,m,ni,ne,f,c,ictrl,lctrl,
     &           rctrl,cctrl,do_diff,iuser,luser,cuser,ruser,eqn,lin,
     &           sloop,acurv,nnzh,xlam)
!----------------------------------------------------------------------!
!                                                                      !
!  Compute the Jacobian and Hessian using differences. Coded for       !
!  clarity, not efficiency.                                            !
!  Please see the users manual for type declarations and comments      !
!                                                                      !
!----------------------------------------------------------------------!
      implicit          none
      include           'ctrl.h'
      logical           eqn(*), lin(*)
      logical           do_diff
      logical           debug_local
      integer           i,j,k,l,m,n,ni,ne,nfe,sloop
      double precision  x(n),dx1(n),dx2(n),temp
      double precision  dx12(n),dx21(n)
      double precision  c_x(m),f_x
      double precision  c_d1(m),f_d1
      double precision  c_d12(m),f_d12,c_d21(m),f_d21
      double precision  c_d2(m),f_d2
      double precision  f,c(m),gf(n),gc(m,n)
      double precision  gf2(n,n),gc2(m,n,n),acurv
      double precision  H_val(nnzh),xlam(m),xlamj

      integer           H_row(nnzh),H_col(nnzh),H_ptr(n+1)
      integer           h_ne,nnzh,k1,l1,m1

!
      data              debug_local /.true./
!
      include          'ctrl_get.inc'
!
      if (sloop.eq.1) then
!
        do_diff=.true.
!
        !if (structure.ne.1) stop ' structure.ne.1 in finDiffGradHessDense'
!
        call SAOi_funcs(n,m,ni,ne,x,f_x,c_x,iuser,luser,
     &                cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
!  Construct the Jacobian using forward differences
!
        nfe = 0
        dx1 = x
!
        do i=1,n
          dx1(i)=x(i)+deltx
          call SAOi_funcs(n,m,ni,ne,dx1,f_d1,c_d1,iuser,luser,
     &                    cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
          nfe=nfe+1
          gf(i)=(f_d1-f_x)/deltx
!
          if (debug_local) write(35,1001)i,gf(i) ! same format as in CheckSAOi_grad.out
!
          do j=1,m
            k=(i-1)*ni+j                     ! for debugging only
            gc(j,i)=(c_d1(j)-c_x(j))/deltx
!
            if (debug_local) write(35,1002)i,j,k,gc(j,i) ! same format as in CheckSAOi_grad.out
!
          end do
          dx1(i)=x(i)
        end do
!
!  Construct the diagonal (lower triangle) Hessian terms using forward differences
!
        !deltxH = 1.d-4
!
        gf2 = 0.d0
        gc2 = 0.d0
!
        dx1 = x
        dx2 = x
        do i=1,n
!
          dx1(i)=x(i)+deltxH
          dx2(i)=x(i)+2.d0*deltxH
!
          call SAOi_funcs(n,m,ni,ne,dx1,f_d1,c_d1,iuser,luser,
     &                    cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
          call SAOi_funcs(n,m,ni,ne,dx2,f_d2,c_d2,iuser,luser,
     &                    cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
          nfe=nfe+2
!           gf2(i,i)=(f_d2-2.d0*f_d1+f_x)/deltxH**2
!
          do j=1,m
            gc2(j,i,i)=(c_d2(j)-2.d0*c_d1(j)+c_x(j))/deltxH**2
          end do
!
          dx1(i)=x(i)
          dx2(i)=x(i)
!
        end do
!
!    Construct the off-diagonal lower triangle Hessian terms using forward differences
!
        k = 0
        dx12 = x
        dx21 = x
        dx2 = x
        do i=1,n
          do l=i+1,n
!
            dx2(i)=x(i)+deltxH
            dx2(l)=x(l)+deltxH
!
            call SAOi_funcs(n,m,ni,ne,dx2,f_d2,c_d2,iuser,luser,
     &                      cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
            dx2(i)=x(i)
            dx2(l)=x(l)
!
            dx12(i)=x(i)+deltxH
!
            call SAOi_funcs(n,m,ni,ne,dx12,f_d12,c_d12,iuser,luser,
     &                      cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
            dx12(i)=x(i)
!
!    f{(x(i))(x(l) + h}
!
            dx21(l)=x(l)+deltxH

            call SAOi_funcs(n,m,ni,ne,dx21,f_d21,c_d21,iuser,luser,
     &                      cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
            dx21(l)=x(l)
!
            nfe=nfe+3
!
!  Now construct the lower triangle
!
           gf2(l,i)=(f_d2-f_d12-f_d21+f_x)/deltxH**2 !not needed in Sand
!
            do j=1,m
              gc2(j,l,i)=(c_d2(j)-c_d12(j)-c_d21(j)+c_x(j))/deltxH**2
            end do
!
          enddo
        end do

        do i=1,n
          do l=1,i
            k=k+1
            ruser(k) = gf2(i,l) !not needed in Sand
          enddo
        enddo
!
        m1=0
        do j=1,m !constraint
          do i=1,n !row
            do l=1,i !column
              k=k+1
              ruser(k) = gc2(j,i,l)
            enddo
          enddo
        enddo
!
        if (k.gt.rmax) then
          write(*,*) 'k.gt.rmax in FinDiffGradHessDense;',
     &               ' rmax is defined in size.h',k,rmax
          stop
        endif
!
        do_diff=.false.
!
      endif
!
      H_val = 0.d0

      k = 0
      k1 = 0
      do i=1,n
        do l1=1,i
          k=k+1
          k1=k1+1
          H_val(k1) = ruser(k)
          H_row(k1) = i
          H_col(k1) = l1
        enddo
      enddo

      do j=1,m! constraint number
        k1 = 0
        xlamj = xlam(j)
        do i=1,n ! row number
          do l1=1,i ! column number
            k=k+1
            k1=k1+1
            H_val(k1) = H_val(k1) + xlamj*ruser(k)
          enddo
        enddo
      enddo

      k1=0
      k=0
      do i=1,n ! row number
        do l1=1,i ! column number
          k=k+1
          if (abs(H_val(k)) .gt. 1.1d-3) then
            k1=k1+1
            ruser(k1) = H_val(k)
            H_row(k1) = H_row(k)
            H_col(k1) = H_col(k)
          endif
        enddo
      enddo

      if (abs(sum(xlam)) .gt. 1.1d-3) then

        write(36,*)'Finite Difference '
        j=5
        write(36,*)'xlam = '
        do i= 1,m-j,j
          write(36,"(10ES12.3)") xlam(i:i+j-1)
        enddo
        write(36,"(100ES12.3)") xlam(i:m)
        write(36,*)''
        write(36,*) 'constraints = ',m
        write(36,*) 'variables   = ',n
        write(36,*) 'nnzh        = ',k1
        write(36,*)''
        write(36,"(a12,a13,a8)") "Hval",  "hrow",  "hcol"
        do i =1,n
          do j= 1,k1
            if ((h_col(j) == i) .and. (abs(ruser(j)) .gt. 2.d-2)) then
              Write(36,"(ES16.4,i8,i8)") ruser(j), H_row(j), H_col(j)
            endif
          enddo
        enddo
      endif

      return
!
 1000 format(1000i4)
 1001 format(1i8,16x,1e22.10)
 1002 format(3i8,1e22.10)
!
      end subroutine FinDiffGradHessSparse
!--------------------------------------------------------------------------


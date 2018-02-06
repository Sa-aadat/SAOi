! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: The secondary SAO driver for the sparse implementation  
! SAOi:


      subroutine sao_sparse(n,ni,ne,x,x_l0,x_u0,timef,timeg,times,
     &                     f,viol,xkkt,kloop,llooptot,
     &                     ictrl,lctrl,rctrl,cctrl,warn,nfe,nge,
     &                     feasible,iuser,luser,cuser,ruser,nnz,nnzh,
     &                     eqn,lin,nw,nm,mxloop,ierr,time0,iseed,xlam,
     &                     iname,sloop,z)
      implicit         none
      include          'ctrl.h'
      logical          eqn(*), lin(*)
      character*90     ostring
      integer          kloopmax,lloopmax,kloop,llooptot,lloop,ig,ns
      integer          iactblo,iactbhi,iactc,ipr,n5,i,j,nact,nfe,nge,n1
      integer          n,ni,ne,iact(ni),ksubitertot,infilter,ksubiter
      integer          outerloop,message,ifree,nw,nm,mxloop,ierr,iseed
      integer          Acol(nnz),Aptr(ni+1),nnz,nnzh ! unity in dense
      integer          storind(nnz),iname,nsx,sloop
      double precision x(n),x0(n),xlam(ni),c_h(ni),h_h(ne),v_h
      double precision h_a(ne),c_a(ni),x_h(n),x_h2(n),f_h,yhi
      double precision gf_h(1),gc_h(1),gh_h(ne,n),x_l0(n),x_u0(n)
      double precision acurv,bcurv(ni),ccurv(ne),x_l(n),x_u(n)
      double precision acurvn(1),bcurvn(1,1),ccurvn(ne,n),drandu1
      double precision f,c(ni),h(ne),gf(n),filter_list(2,nmax)
      double precision gc(nnz),gh(ne,n),cplus_norm_h,timey
      double precision timef1,timeg1,times1,timef2,timeg2,times2
      double precision eps_locala,eps_localb,rfd,dml_inf0,rho
      double precision current_tol,current_kkt,f_a,cplus_norm,norm2b
      double precision SAOi_seconds,epsmch,dpmeps,viol,drange,timex
      double precision pred_f,real_f,qindicator,delxnormc,xkkt,flam
      double precision xlamsml,xlambig,phibw,timef,timeg,times,xkktl
      double precision delxnormi,norm2f!,Li(n),Ui(n),Li_h(n),Ui_h(n)
      double precision hsum,xlam_h(ni),s(n),time0,time_tmp
      double precision amult,bmult(ni),z(n),z_h(n)
      double precision delta_t, t_g1, t_g2, gelta_t
      
      integer          Hrow(nnzh),Hcol(nnzh),Hptr(n+1),Hne
      double precision Hval(nnzh)
      
      logical          warn,finalstg,accept,cstage0,init,do_diff
      logical          feasible
      !save             timef1,timeg1,times1,timef2,timeg2,times2
      include          'ctrl_get.inc'
!
      do i = 1,nnz
        storind(i)=i-1
      enddo

      delta_t = 0
      gelta_t = 0

!  remove below before releasing Version 1.0.0, or as many as possible,
!     starting with 2, 1, 6, 4, 3
!
!  resolved: 5,
!
!       if (finite_diff.or.check_grad)    stop ' illegal in sao_sparse 1'
!
!       if (force_converge.ne.0
!      &        .and.force_converge.ne.2) stop ' illegal in sao_sparse 2'
!
      if (approx_f.ne.4.and.approx_f.ne.1.and.approx_f.ne.14
     &    .and.approx_f.lt.100)
     &                                 stop ' illegal in sao_sparse 3'
!
      if (approx_c.ne.4.and.approx_c.ne.1.and.approx_c.ne.14
     &    .and.approx_c.lt.100)
     &                                 stop ' illegal in sao_sparse 4'
!
      if (subsolver.ne.1
     &        .and.subsolver.ne.21
     &        .and.subsolver.ne.25
     &        .and.subsolver.ne.26
     &        .and.subsolver.ne.27
     &        .and.subsolver.ne.28)    stop ' illegal in sao_sparse 6'
!


!  initialize counters and a few other thingies
      n1              = n+1
      message         = -200
      init            = .true.
      yhi             = 0.d0
      kloopmax        = outermax
      lloopmax        = innermax
      kloop           = 0
      lloop           = 0
      llooptot        = 0
      outerloop       = 0
      eps_locala      = 1.d-7
      eps_localb      = 1.d-7
      ksubitertot     = 0
      rfd             = big
      n5              = 5
      cstage0         = .false.
      finalstg        = .false.
      dml_inf0        = dml_infinity
      rho             = dml_infinity
      do i=1,n 
        x0(i)         = x(i)
        x_l(i)        = x_l0(i)
        x_u(i)        = x_u0(i)
      enddo
      current_tol     = 1.d0
      current_kkt     = 1.d0
      epsmch          = dpmeps()
      time_tmp        = time0
      amult           = 1.d0
      acurv           = 0.d0
      do j=1,ni
        bmult(j)      = 1.d0
        bcurv(j)      = 0.d0
      enddo
!
!  write a few headers
!       if (iname.eq.0) then
!         write(6,*)' '
!         write(6,*) 'Problem name: ',cname1 
!         write(9,*)' '
!         write(9,*) 'Problem name: ',cname1 
!       endif  
      if (iprint.ge.1) write(6,6010) subsolver,cname1
      write(8,5000)
      if (iprint.ge.1) write(9,6010) subsolver,cname1
      write(10,7500)
      write(11,7500)
      if (debug) write(12,6001)
      write(13,6002)
      if (nprob.ge.180.and.nprob.le.189) then
        open (15,file='Specialized.out')
        write(15,9500)
      endif
      if (check_grad) open (17,file='CheckSAOi_grad.out')
      open (55,file='Priv.out')

!  determine function and constraint values at starting point
      timef1 = SAOi_seconds()
      call SAOi_funcs (n,ni,ni-ne,ne,x,f,c,iuser,luser,cuser,ruser,
     &                 eqn,lin,ictrl,lctrl,rctrl,cctrl)
      timef2 = SAOi_seconds()
      nfe=nfe+1
      timef=timef+(timef2-timef1)

!  determine gradients of function and constraints at starting point
      timeg1 = SAOi_seconds()
      call formgrad (n,x,ni,ne,f,c,gf,gc,ictrl,lctrl,rctrl,cctrl,
     &               nfe,nge,do_diff,iuser,luser,cuser,ruser,nnz,
     &               Acol,Aptr,eqn,lin)
      timeg2 = SAOi_seconds()
      nge=nge+1
      timeg=timeg+(timeg2-timeg1)
      
!  calculate the active set on the sao level
      call form_act(ni,c,iact,nact,ictrl,lctrl,rctrl,cctrl,
     &              iuser,luser,cuser,ruser)

!  initialize the dual variables for the constraints
      do i=1,ni
        xlam_h(i)=xlam(i) ! xlam is initialized in SAOi_check for the time being
      enddo

!  initialize the dual variables for the simple bounds (box constraints)
      do i=1,n
        z_h=z(i) ! z is initialized in SAOi_check for the time being 
      enddo

!  determine if starting point is feasible
      if (unconstrained) then
        call SAOistatus(n,ni,iactblo,iactbhi,iactc,cplus_norm,xlam,
     &                  viol,c,x,x_l0,x_u0,5,1,hsum,ictrl,lctrl,
     &                  rctrl,cctrl,feasible,
     &                  iuser,luser,cuser,ruser)
      else
        call SAOistatus(n,ni,iactblo,iactbhi,iactc,cplus_norm,xlam,
     &                  viol,c,x,x_l0,x_u0,1,1,hsum,ictrl,lctrl,
     &                  rctrl,cctrl,feasible,
     &                  iuser,luser,cuser,ruser)
        call SAOistatusE(n,ni,iactblo,iactbhi,iactc,cplus_norm,xlam,
     &                   viol,c,x,x_l0,x_u0,1,1,hsum,ictrl,lctrl,
     &                   rctrl,cctrl,feasible,eqn,lin,
     &                   iuser,luser,cuser,ruser)
      endif
 
!  output some subproblem stuff
      ostring = ' the starting point '
      if (debug) write(12,8500) 0,0,0,0.d0,0.d0,f,viol,ostring

!  conditionally initialize cstage0
      if (viol.gt.feaslim) cstage0=.true.

!  construct the approximate diagonal hessian
      amult           = 1.d0
      do j=1,ni
        bmult(j)      = 1.d0
      enddo
!
      call diaHess_sq1s (n,ni,ne,acurv,bcurv,ccurv,x,x_h,
     &                   gf,gc,gh,f,c,h,f_h,c_h,h_h,
     &                   ictrl,lctrl,rctrl,cctrl,outerloop,
     &                   Acol,Aptr,nnz,
     &                   iuser,luser,cuser,ruser) 
     
!  form the sparse COO or CSR Hessian
      if (subsolver.ge.25.and.subsolver.le.28) then
      
        include './ProbDevelop.inc' ! Temporarily for xAx galahad test problems - hard coded
 
!  form the exact sparse COO Hessian if A is not diagonal - general case
        if (FormHessSparse) then
          call FormHessSparseN (n,ni,x,xlam,Hne,Hval,Hrow,Hcol,
     &                          ictrl,lctrl,rctrl,cctrl,iuser,luser,
     &                          cuser,ruser,nnzh,eqn,lin,acurv) 
        endif
!                
      endif

!  store first iterate
      ig=0
      if (approx_c.eq.2.or.approx_c.eq.3.or.approx_c.eq.5) ig=1
      call store_iters (n,ni,f_h,f,v_h,viol,x_h,x_h2,x,c_h,c,
     &                  cplus_norm,cplus_norm_h,gf_h,gf,gc_h,gc,
     &                  xlam,xlam_h,1,ig,z,z_h,
     &                  iuser,luser,cuser,ruser,ictrl,lctrl,
     &                  rctrl,cctrl)

!  write some output
      if (iprint.ge.1) write (*,7010) kloop,lloop,f,viol,feasible
      write (9,7010) kloop,lloop,f,viol,feasible

!  possibly do specialized things
      call SAOi_Special (n,ni,ne,x,f,c,h,iuser,luser,cuser,ruser,
     &                   eqn,lin,ictrl,lctrl,rctrl,cctrl,x_l0,x_u0,
     &                   kloop,.false.,xkkt,viol,iactblo,iactbhi,
     &                   iactc,delxnormc,delxnormi,0.d0)

!  prepare for the outer optimization loops
      kloop=1

!  flush a buffer the fortran way...
      close(9)
      open (9 ,file='History.out',status='old',position='append')
      if (debug) then
        close(12)
        open (12 ,file='Subproblems.out',status='old',position='append')
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                start the outer optimization loop                     !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do while (kloop.le.kloopmax)
!
        rctrl(34)=drandu1(iseed)
        outerloop=kloop
        lloop=0
!
        rho=dml_inf0
!
 110    continue
        do i=1,n
          drange=min(rangemax,rho*(x_u0(i)-x_l0(i)))
          x_l(i)=max(x(i)-drange,x_l0(i))
          x_u(i)=min(x(i)+drange,x_u0(i))

          !write(*,*) drange
          !write(*,*) x_l(1)
          !write(*,*) x_u(1)
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!               possibly start an inner conservative loop              !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  initialize and check for unconstrained problems


      if (GetFinDiffGradHessDense) then !From albert email
      i = int(n*(n+1)/2)!n*n !
        call FinDiffGradHessSparse (n,x,ni,ni-ne,ne,f,c,ictrl,lctrl,
     &                            rctrl,cctrl,do_diff,iuser,luser,
     &                            cuser,ruser,eqn,lin,sloop,acurv,
     &                            i,xlam)
!         stop
      endif
        ksubiter=0
        times1 = SAOi_seconds()
        if (unconstrained) then

!         call uncstr (n,x,x_h,acurvn,x_l,x_u,gf,ksubiter) ! no advantage over dense !

        else
          if (subsolver.eq.1) then
!  call the l-bfgs-b solver to solve the Falk dual for constrained PRIMAL problems (sparse)
            call drive_lbfgsb24f_l (n,x,ni,ne,f_a,c_a,h_a,iact,nact,
     &                         xlam,x_h,acurv,bcurv,ccurv,
     &                         acurvn,bcurvn,ccurvn,
     &                         x_l,x_u,f,c,h,gf,gc,gh,ksubiter,
     &                         xkktl,flam,cstage0,ostring,xlam_h,
     &                         ns,ictrl,lctrl,rctrl,cctrl,message,
     &                         nnz,Acol,Aptr,yhi,outerloop,eqn,lin,
     &                         amult,bmult,
     &                         iuser,luser,cuser,ruser)
!
            elseif (subsolver.eq.21) then
!
!  call the Cplex solver to solve the equivalent diagonal QP problem
            t_g1 = SAOi_seconds();
            call drive_cplexs_old  (n,x,ni,ne,f_a,c_a,
     &                         iact,nact,xlam,x_h,
     &                         acurv,bcurv,amult,bmult,x_l,x_u,
     &                         f,c,gf,gc,ksubiter,
     &                         xkkt,flam,cstage0,ostring,xlam_h,
     &                         nsx,ictrl,lctrl,rctrl,cctrl,message,yhi,
     &                         outerloop,lloop, nnz,Acol,Aptr,
     &                         storind,delta_t,eqn,
     &                         iuser,luser,cuser,ruser,nnzh)
!     
          elseif  (subsolver.ge.25.and.subsolver.le.28) then
!
! call the Galahad QP solvers to solve the equivalent diagonal QP problem
            call drive_galQPs (n,x,ni,ne,f_a,c_a,iact,nact,xlam,
     &                          x_h,acurv,bcurv,
     &                          x_l,x_u,f,c,gf,gc,ksubiter,
     &                          xkktl,flam,cstage0,ostring,
     &                          ictrl,lctrl,rctrl,cctrl,outerloop,
     &                          nnz,nnzh,Acol,Aptr,yhi,message,eqn,
     &                          lin,z,Hne,Hval,Hrow,Hcol,
     &                          iuser,luser,cuser,ruser)
!
            else
!            
            stop ' subsolver not yet implemented for sparse '
            ! but some are not very interesting anyway...
!            
          endif
        endif


! catch floating point exceptions
        call testNaN(f_a, ierr)
        if (ierr.lt.0) then
          !accept = .false.
          !trust_region = .true.
          !goto 1233
          write(*,*) ' NaN occured in sao_sparse '
          stop
        endif
        
!  do some timing
        times2 = SAOi_seconds()
        times=times+(times2-times1)
        ksubitertot=ksubitertot+ksubiter

!  determine true function and constraint values
        timef1 = SAOi_seconds()
        call SAOi_funcs (n,ni,ni-ne,ne,x,f,c,iuser,luser,cuser,ruser,
     &                   eqn,lin,ictrl,lctrl,rctrl,cctrl)
        timef2 = SAOi_seconds()
        nfe=nfe+1
        timef=timef+(timef2-timef1)

! catch floating point exceptions
!  determine if feasible
        if (unconstrained) then
          call SAOistatus(n,ni,iactblo,iactbhi,iactc,cplus_norm,xlam,
     &                   viol,c,x,x_l0,x_u0,5,-1,hsum,ictrl,lctrl,
     &                   rctrl,cctrl,feasible,
     &                   iuser,luser,cuser,ruser)
        else
          call SAOistatus(n,ni,iactblo,iactbhi,iactc,cplus_norm,xlam,
     &                   viol,c,x,x_l0,x_u0,1,-1,hsum,ictrl,lctrl,
     &                   rctrl,cctrl,feasible,
     &                   iuser,luser,cuser,ruser)
          call SAOistatusE(n,ni,iactblo,iactbhi,iactc,cplus_norm,xlam,
     &                    viol,c,x,x_l0,x_u0,1,1,hsum,ictrl,lctrl,
     &                    rctrl,cctrl,feasible,eqn,lin,
     &                    iuser,luser,cuser,ruser)
        endif

        accept=.true.

!  conditionally kill cstage0
        if (viol.lt.feaslim) cstage0=.false.                      ! agap-2008

!  output some subproblem stuff
        if (debug)
     &    write(12,8500) kloop,lloop,ksubiter,flam,xkktl,f,viol,ostring

!   escape if relative step size is too small
        delxnormc = norm2b(n,x_h,x)
        delxnormi = norm2f(n,x_h,x)
        if (delxnormc.lt.xtol.or.delxnormi.lt.xtol_inf) goto 1234 ! do not enforce conservatism

!  escape if a feasible descent step was made
        if (conservative.and.f.lt.f_h.and.unconstrained) then
          goto 1234                      ! do not enforce conservatism
        endif

!  escape if a feasible descent step was made
        if (conservative.and.f.lt.f_h.and.feasible.and.allow_f) then
          goto 1234                      ! do not enforce conservatism
        endif

!  escape if an infeasible restoration step was made
        if (conservative.and.viol.lt.v_h.and.cstage0.and.allow_c) then
          goto 1234                      ! do not enforce conservatism
        endif                                                          ! agap-2008

!  determine if a conservative inner loop is required
        if (conservative.and.kloop.gt.-1) then
          call conserve2 (n,ni,accept,f,f_a,c,c_a,eps_locala,
     &                   eps_localb,acurv,bcurv,
     &                   rho,ictrl,lctrl,rctrl,cctrl,lloop,amult,bmult,
     &                   iuser,luser,cuser,ruser)
        endif

!  determine if a trust region step is required
        if (trust_region) then
          call filter_ises (n,ni,accept,f,f_h,f_a,viol,v_h,
     &                      c,c_h,c_a,kloop,acurv,
     &                      acurvn,bcurv,bcurvn,
     &                      eps_locala,eps_localb,
     &                      filter_list,infilter,init,rho,
     &                      ictrl,lctrl,rctrl,cctrl,lloop,amult,bmult,
     &                      iuser,luser,cuser,ruser)
          if (rho.eq.rho_l_min) then
            goto 1234 ! no further reduction possible
          endif
        endif

!  increment inner loop counters
 1233   if (.not.accept.and.lloop.lt.lloopmax)  then
          lloop=lloop+1
          llooptot=llooptot+1

!  restore iterate
          call store_iters (n,ni,f_h,f,v_h,viol,x_h,x_h2,x,c_h,c,
     &                      cplus_norm,cplus_norm_h,gf_h,gf,gc_h,gc,
     &                      xlam,xlam_h,-1,0,z,z_h,
     &                      iuser,luser,cuser,ruser,ictrl,lctrl,
     &                      rctrl,cctrl)
          goto 110
        endif
!
 1234   continue ! escape from the inner loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                      end the inner conservative loop                 !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  determine if feasible
        if (unconstrained) then
          call SAOistatus(n,ni,iactblo,iactbhi,iactc,cplus_norm,xlam,
     &                    viol,c,x,x_l0,x_u0,5,1,hsum,ictrl,lctrl,
     &                    rctrl,cctrl,feasible,
     &                    iuser,luser,cuser,ruser)
        else
          call SAOistatus(n,ni,iactblo,iactbhi,iactc,cplus_norm,xlam,
     &                    viol,c,x,x_l0,x_u0,1,1,hsum,ictrl,lctrl,
     &                    rctrl,cctrl,feasible,
     &                    iuser,luser,cuser,ruser)
          call SAOistatusE(n,ni,iactblo,iactbhi,iactc,cplus_norm,xlam,
     &                    viol,c,x,x_l0,x_u0,1,1,hsum,ictrl,lctrl,
     &                    rctrl,cctrl,feasible,eqn,lin,
     &                    iuser,luser,cuser,ruser)
        endif
!  evaluate some filter thingies
        pred_f=f_h-f_a
        real_f=f_h-f
        if (pred_f.ne.0.d0) then
          qindicator=real_f/pred_f
        else
          qindicator=0.d0
        endif

!  calculate relative step sizes
        delxnormc = norm2b(n,x_h,x)
        delxnormi = norm2f(n,x_h,x)
        if (delxnormc.lt.xtol*1.d2.or.delxnormi.lt.xtol_inf*1.d2) then
          finalstg=.true.
          rctrl(17)=0.d0 ! btol1
        endif

!  calculate relative function difference
        rfd=(f_h-f)/(1.d0+dabs(f_h))

!  determine gradients
        timeg1 = SAOi_seconds()
        call formgrad (n,x,ni,ne,f,c,gf,gc,ictrl,lctrl,rctrl,cctrl,
     &                 nfe,nge,do_diff,iuser,luser,cuser,ruser,nnz,
     &                 Acol,Aptr,eqn,lin)
        timeg2 = SAOi_seconds()
        nge=nge+1
        timeg=timeg+(timeg2-timeg1)

!  construct the approximate diagonal hessian
        amult           = 1.d0
        do j=1,ni
          bmult(j)      = 1.d0
        enddo
!   
        call diaHess_sq1s (n,ni,ne,acurv,bcurv,ccurv,x,x_h,
     &                     gf,gc,gh,f,c,h,f_h,c_h,h_h,
     &                     ictrl,lctrl,rctrl,cctrl,outerloop,
     &                     Acol,Aptr,nnz,
     &                     iuser,luser,cuser,ruser)

!  form the sparse COO or CSR Hessian
        if (subsolver.ge.25.and.subsolver.le.28) then

          include './ProbDevelop.inc' ! Temporarily for xAx galahad test problems - hard coded
        
!  form the exact sparse COO Hessian if A is not diagonal - general case
          if (FormHessSparse) then
          
            write(*,*) 'FormHessSparse'
          
            call FormHessSparseN (n,ni,x,xlam,Hne,Hval,Hrow,Hcol,
     &                            ictrl,lctrl,rctrl,cctrl,iuser,luser,
     &                            cuser,ruser,nnzh,eqn,lin,acurv) 
          endif
!                
        endif
       
!  store iterate
        ig=0
        if (approx_c.eq.2.or.approx_c.eq.3.or.approx_c.eq.5) ig=1
        call store_iters (n,ni,f_h,f,v_h,viol,x_h,x_h2,x,c_h,c,
     &                    cplus_norm,cplus_norm_h,gf_h,gf,gc_h,gc,
     &                    xlam,xlam_h,1,ig,z,z_h,
     &                    iuser,luser,cuser,ruser,ictrl,lctrl,
     &                    rctrl,cctrl)

!  calculate the active set on the subproblem level
        call form_act(ni,c,iact,nact,ictrl,lctrl,rctrl,cctrl,
     &                iuser,luser,cuser,ruser)

!  calculate the true KKT conditions
        call form_kkts (xkkt,n,ni,ne,x,iact,nact,x_l0,x_u0,
     &                  gf,gc,gh,xlam,xlamsml,xlambig,nnz,
     &                  Acol,Aptr,ifree,
     &                  ictrl,lctrl,rctrl,cctrl,z,
     &                  iuser,luser,cuser,ruser)

! some timing
      timex = SAOi_seconds ()
      timey = timex - time_tmp

! write some output
        if (iprint.ge.1) write (*,7012) kloop,lloop,f,viol,feasible,
     &  finalstg,rho,qindicator,xkkt,delxnormc,rfd,
     &  iactblo,iactbhi,iactc,message,timey

! write some more output
        write (9,7012) kloop,lloop,f,viol,feasible,
     &  finalstg,rho,qindicator,xkkt,delxnormc,rfd,
     &  iactblo,iactbhi,iactc,message,timey

! some more timing
        time_tmp = timex

! write even more output
        if (delxnormc.lt.current_tol) then
          write(10,7600)current_tol,kloop,llooptot
          do i=1,15
            if (current_tol/10.d0.lt.delxnormc) then
              current_tol=current_tol/10.d0
              exit
            else
              current_tol=current_tol/10.d0
              write(10,7600)current_tol,kloop,llooptot
            endif
          enddo
        endif
        if (xkkt.lt.current_kkt) then
          write(11,7600)current_kkt,kloop,llooptot
          do i=1,15
            if (current_kkt/10.d0.lt.xkkt) then
              current_kkt=current_kkt/10.d0
              exit
            else
            current_kkt=current_kkt/10.d0
            write(11,7600)current_kkt,kloop,llooptot
            endif
          enddo
        endif
!  exit do loop here on final iteration
        if (kloop.eq.kloopmax) mayStop=.true.
        if (feasible) mayStop=.true.
        if (mayStop.and.(kloop.eq.kloopmax
     &      .or. (delxnormc.le.xtol.and.current_kkt.le.kkt_min)
     &      .or. (delxnormi.le.xtol_inf.and.current_kkt.le.kkt_min)
     &      .or. current_kkt.le.kkt_tol)) then
          if (iprint.ge.1) write(*,9600)
          write(9,*) ' '
          if (kloop.eq.kloopmax) then
            if (iprint.ge.1) write(*,6012)
            write(9,6012)
            call open_warn_file (warn)
            write (14,6013)
          end if
          if (delxnormc.le.xtol.and.current_kkt.le.kkt_min
     &                                             .and.mayStop) then
            if (iprint.ge.1) write(*,6014)
            write(9,6014)
          end if
          if (current_kkt.le.kkt_tol) then
            if (iprint.ge.1) write(*,6015)
            write(9,6015)
          end if
          if (delxnormi.le.xtol_inf.and.current_kkt.le.kkt_min
     &                                             .and.mayStop) then
            if (iprint.ge.1) write(*,6016)
            write(9,6016)
          end if
          if (iprint.ge.1) write(*,9600)
          write(9,*) ' '
          exit
        endif           ! stopping condition
!
        if (dabs(rfd).lt.ftol) then         ! ftol is normally negative !
          if (iprint.ge.1) write(*,9600)
          write(9,*) ' '
          if (iprint.ge.1) write(*,6017)
          write(9,6017)
          if (iprint.ge.1) write(*,9600)
          write(9,*) ' '
          exit
        endif          ! stopping condition on relative function difference

!  possibly do specialized things
        call SAOi_Special (n,ni,ne,x,f,c,h,iuser,luser,cuser,ruser,
     &                     eqn,lin,ictrl,lctrl,rctrl,cctrl,x_l0,x_u0,
     &                     kloop,.false.,xkkt,viol,iactblo,iactbhi,
     &                     iactc,delxnormc,delxnormi,timey)

        kloop=kloop+1    !increment outer loop counter
!

!  again flush some buffers the fortran way...
        close(9)
        open (9 ,file='History.out',status='old',position='append')
        if (debug) then
          close(12)
          open (12,file='Subproblems.out',status='old',
     &          position='append')
        endif
!
      end do             ! while kloop.lt.kloopmax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                   end of outer optimization loop                     !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  write a final summary
      do ipr=6,9,3
        if (ipr.eq.6.and.iprint.ge.2.or.ipr.eq.9) then
          write(ipr,*) ' number of variables n     : ',n
          write(ipr,*) ' number of inequalities ni : ',ni-ne
          if (ne.gt.0)
     &    write(ipr,*) ' number of equalities ne   : ',ne
          write(ipr,*) ' '
          write(ipr,*) ' Outer iterations k used   : ',kloop
          write(ipr,*) ' Inner iterations l used   : ',llooptot
          write(ipr,*) ' Subproblem evaluations    : ',ksubitertot
          write(ipr,*) ' '
          write(ipr,*) ' Machine precision         : ',epsmch
          write(ipr,*) ' '
          write(ipr,8101) f,viol
          write(ipr,8102) (x(i),i=1,min(n,n5))
          if (ni.gt.0) then
            write(ipr,8103) (c(i),i=1,min(ni,n5))
            write(ipr,8104) (xlam(i),i=1,min(ni,n5))
          endif
        endif
      enddo
      write(*,*) ' '

!  calculate the true KKT conditions
      call form_kkts (xkkt,n,ni,ne,x,iact,nact,x_l0,x_u0,
     &                  gf,gc,gh,xlam,xlamsml,xlambig,nnz,
     &                  Acol,Aptr,ifree,
     &                  ictrl,lctrl,rctrl,cctrl,z,
     &                  iuser,luser,cuser,ruser)

!  possibly do specialized things
        call SAOi_Special (n,ni,ne,x,f,c,h,iuser,luser,cuser,ruser,
     &                     eqn,lin,ictrl,lctrl,rctrl,cctrl,x_l0,x_u0,
     &                     kloop,.true.,xkkt,viol,iactblo,iactbhi,
     &                     iactc,delxnormc,delxnormi,timey)

      do ipr=6,9,3
        if (ipr.eq.6.and.iprint.ge.2.or.ipr.eq.9) then
          if (ni.gt.0) then
            write(ipr,*)' '
            write(ipr,8900) xlamsml
            write(ipr,8901) xlambig
          endif
          write(ipr,*)' '
          write(ipr,8902) xkkt
        endif
      enddo
!
      if (yhi.gt.1.d-8) then
        call open_warn_file (warn)
        write (14,8801)yhi
      endif
!
      if (nprob.ge.180.and.nprob.le.189) then
        phiBW=dble(iactblo+iactbhi)/dble(n)
        do ipr=6,9,3
          if (ipr.eq.6.and.iprint.ge.2.or.ipr.eq.9) then
            write(ipr,9000) f,kloop,llooptot,phiBW
          endif
        enddo
      endif
!
      do i=1,n
        write(8,5001) i,x0(i),x(i),x_l0(i),x_u0(i)
      enddo
!
      do j=1,ni
        write(13,8501) j,c(j),xlam(j)
        if (xlam(j).ge.biglam) then
          write(9,8600) j
          write(13,8601) j
          call open_warn_file (warn)
          write(14,8602) j,biglam
        endif
      enddo
!
!
!
!
      return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 5000 format (/,20x,'SAOi algorithm: primal variables',
     &     //,100('-'),/,1x,
     &    'Component  starting point     final point  ',
     &    '   lower bound     upper bound ',/,100('-'))
 5001 format (i10,6es16.6)
!
 6000 format (/,46x,'SAOi algorithm for Problem: ',a24,//,
     &    130('-'),/,1x,
     &    '  Outer Inner   Function val  Max constr Stats    ',
     &    ' dml       Qual      xnorm       Frel      ActLo      ActHi',
     &    '       ActC  Message',/,130('-'))
 6001 format (/,26x,'SAOi algorithm',//,137('-'),/,1x,
     & '  Outer Inner SubProb        gamma        ',
     & '  skkt             f             h',
     & '   comments',/,137('-'))
 6002 format (/,9x,'SAOi algorithm: constraints and ',
     & 'dual variables',//,77('-'),/,1x,
     & 'Number                  Constraint            ',
     & '   Dual variable  ',/,77('-'))
 6010 format (/,40x,'SAOi Algorithm with Subsolver',i3,
     &    ' for Problem: ',a24,//,146('-'),/,1x,
     &    '  Outer Inner   Function val  Max viol  S    ',
     &    '   dml      Qual       kkt     xnorm      Frel  ',
     &    '    ActLo      ActHi',
     &    '       ActC  Message      CPU(s)',/,146('-'))
!
 6012 format ('  Terminated on maximum number of steps ')
 6013 format ('  Convergence criteria not satisfied - terminated on ',
     &        'maximum number of steps ')
 6014 format ('  Terminated on x-tolerance (Euclidian norm) ')
 6015 format ('  Terminated on KKT-residual ')
 6016 format ('  Terminated on x-tolerance (infinity norm) ')
 6017 format ('  Terminated on relative f-value ')
!
 7000 format (2(1x,i6),1x,es14.7,1x,es11.4,1x,l1)
 7002 format (2(1x,i6),1x,es14.7,1x,es11.4,1x,
     &        2l1,1x,es10.3,1x,es10.3,1x,es10.3,1x,es10.3,
     &        1x,i10,1x,i10,1x,i10,5x,i4)
!
 7010 format (2(1x,i6),1x,es14.7,1x,es9.2,1x,l1)
 7012 format (2(1x,i6),1x,es14.7,1x,es9.2,1x,
     &        2l1,1x,es9.2,1x,es9.2,1x,es9.2,1x,es9.2,1x,es9.2,
     &        1x,i10,1x,i10,1x,i10,5x,i4,f12.2)
!
 7500 format (/,46x,'SAOi algorithm',//,109('-'),/,1x,
     &    '  Tolerance       k       l ',/,109('-'))
 7600 format (1es12.6,2i8)
!
 8101 format (' f, viol   : ',10es18.8)
 8102 format (' x1, ...   : ',10es18.8)
 8103 format (' c1, ...   : ',10es18.8)
 8104 format (' lam1, ... : ',10es18.8)
!
 8500 format (3i7,4es14.6,3x,'-',a90)
 8501 format (1i7,2es28.12)
!
 8600 format (/,' WARNING: Dual variable ',i10,
     &          ' is on its upper bound. Severity = 10 ',33('.'),'!')
 8601 format (  ' WARNING: Dual variable ',i10,
     &          ' is on its upper bound. Severity = 10 ')
 8602 format (  ' WARNING: Dual variable ',i10,
     &          ' is on its upper bound',1es18.6,' Severity = 10 ')
!
 8800 format (/,' WARNING: relaxation variable ',i10,
     &          ' is not equal to zero, but equals: ',1es18.6,
     &          ' Severity = 7 ',33('.'),'!')
 8801 format (  ' WARNING: the relaxation variable',
     &          ' is not equal to zero, but equals: ',1es18.6,
     &          ' Severity = 7 ')
!
 8900 format (' Smallest dual variable: ',1es18.8)
 8901 format (' Largest  dual variable: ',1es18.8)
 8902 format (' KKT residual: ',1es18.8)
!
 9000 format (/,/,'  Customized output: ',1f12.2,2i7,1f8.3,10f12.4,/)
 9500 format ('    It              f              h       b&w    p-SIMP'
     &        ,'     q-GSS      brns',//,77('-'))
 9600 format (5x)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end subroutine sao_sparse
!----------------------------------------------------------------------c

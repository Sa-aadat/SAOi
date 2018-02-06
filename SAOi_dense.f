! SAOi: 
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: The secondary SAO driver. Beware, some subsolvers are highly
! SAOi: experimental, and not for general use. Please refer to the
! SAOi: manual in case of doubt
! SAOi:


      subroutine sao_dense (n,ni,ne,x,x_l0,x_u0,timef,timeg,times,
     &                      f,viol,xkkt,kloop,llooptot,
     &                      ictrl,lctrl,rctrl,cctrl,warn,nfe,nge,
     &                      feasible,iuser,luser,cuser,ruser,shift,
     &                      nnz,nnzh,eqn,lin,nw,nm,mxloop,ierr,time0,
     &                      iseed,xlam,iname,sloop,z)
      implicit         none
      include          'ctrl.h'
      character*90     ostring
      integer          kloopmax,lloopmax,kloop,llooptot,lloop,ig,ns
      integer          iactblo,iactbhi,iactc,ipr,n5,i,j,nact,nfe,nge
      integer          n,ni,ne,iact(ni),ksubitertot,infilter,ksubiter
      integer          outerloop,message,ifree,nw,mxloop,ierr,iname
      integer          Acol(1),ncol(1),n1,nm,nnz,nnzh,sloop
      integer          ivec,lenr,iseed,maxit,maxeval,maxLS
      double precision x(n),x0(n),xlam(ni),c_h(ni),h_h(ne),v_h,pen
      double precision h_a(ne),c_a(ni),x_h(n),x_h2(n),f_h,drandu1
      double precision gf_h(n),gc_h(nnz),gh_h(ne,n),x_l0(n),x_u0(n)
      double precision acurv,bcurv(ni),ccurv(ne),x_l(n),x_u(n)
      double precision acurvn(n),bcurvn(ni,n),ccurvn(ne,n)
      double precision f,c(ni),h(ne),gf(n),filter_list(2,nmax)
      double precision gc(nnz),gh(ne,n),yhi,cplus_norm_h,timey
      double precision timef1,timeg1,times1,timef2,timeg2,times2,ss(ni)
      double precision eps_locala,eps_localb,rfd,dml_inf0,rho,s_h(ni)
      double precision current_tol,current_kkt,f_a,cplus_norm,norm2b
      double precision SAOi_seconds,epsmch,dpmeps,viol,drange,timex
      double precision pred_f,real_f,qindicator,delxnormc,xkkt,flam
      double precision xlamsml,xlambig,phibw,timef,timeg,times,xkktl
      double precision Li(n),Ui(n),Li_h(n),Ui_h(n),delxnormi,norm2f
      double precision hsum,xlam_h(ni),s(n)
      double precision time0,time_tmp 
      double precision y(ni,n),shift(*),xkkt_mp,z(n),z_h(n)
      logical          warn,finalstg,accept,cstage0,init,do_diff
      logical          feasible,eqn(*),lin(*)
      integer          status
      !save             timef1,timeg1,times1,timef2,timeg2,times2
      include          'ctrl_get.inc'
     

!  initialize counters and a few other thingies
      ierr            = 0
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
      acurv           = 0.d0
      do j=1,ni
        bcurv(j)      = 0.d0
      enddo
      do j=1,ne
        ccurv(j)      = 0.d0
      enddo
      do i=1,n
        acurvn(i)     = 0.d0
      enddo
      do i=1,n
        do j=1,ni
          bcurvn(j,i) = 0.d0
        enddo
      enddo
      do i=1,n
        do j=1,ne
          ccurvn(j,i) = 0.d0
        enddo
      enddo
      do j=1,ni
        ss(j)         = 0.d0
        s_h(j)        = 0.d0
      enddo
      time_tmp=time0

!  write a few headers
!       if (iname.eq.0) then
!         write(6,*)' '
!         write(6,*) 'Problem name: ',cname1 
!         write(9,*)' '
!         write(9,*) 'Problem name: ',cname1 
!       endif  
      if (iprint.ge.1) write(6,6010) subsolver,cname1
      if (iprint.ge.3.or.multimodal) write(8,5000)
      if (iprint.ge.1) write(9,6010) subsolver,cname1
      write(10,7500)
      write(11,7500)
      if (debug) write(12,6001)
      if (iprint.ge.3.or.multimodal) write(13,6002)
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
     &               Acol,ncol,eqn,lin)
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
     &                  viol,c,x,x_l0,x_u0,1,1,hsum,ictrl,lctrl,
     &                  rctrl,cctrl,feasible,eqn,lin,
     &                  iuser,luser,cuser,ruser)
      endif

!  adjust the relaxation variable upper bound
      ymax = max(0.d0,viol)
      rctrl(25) = ymax

!  output some subproblem stuff
      ostring = ' the starting point '
      if (debug) write(12,8500) 0,0,0,0.d0,0.d0,f,viol,ostring

!  conditionally initialize cstage0
      if (viol.gt.feaslim) cstage0=.true.

!  construct approximate diagonal hessian
      call diaHess(n,x,ni,ne,x_h,x_h2,acurv,bcurv,ccurv,shift,
     &             acurvn,bcurvn,ccurvn,x_l,x_u,f,c,h,gf,gc,gh,
     &             gf_h,gc_h,gh_h,Li,Ui,Li_h,Ui_h,f_h,c_h,h_h,s,
     &             ictrl,lctrl,rctrl,cctrl,outerloop,iuser,luser,
     &             cuser,ruser)

!  optionally determine first and second order derivatives here - for debugging only... ! 
      if (GetFinDiffGradHessDense) then
        call FinDiffGradHessDense (n,x,ni,ni-ne,ne,f,c,ictrl,lctrl,
     &                             rctrl,cctrl,do_diff,iuser,luser,
     &                             cuser,ruser,eqn,lin,sloop) 
      endif

!  store first iterate
      ig=0
      if (approx_c.eq.2.or.approx_c.eq.3.or.approx_c.eq.5
     &                                  .or.approx_c.eq.8) ig=1
      call store_iter (n,ni,ne,f_h,f,v_h,viol,x_h,x_h2,x,c_h,c,
     &                 cplus_norm,cplus_norm_h,gf_h,gf,gc_h,gc,
     &                 h_h,h,gh_h,gh,xlam,xlam_h,1,ig,z,z_h,
     &                 iuser,luser,cuser,ruser,ictrl,lctrl,
     &                 rctrl,cctrl)

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

      if (outermax.eq.0) then
        write(*,*) ' '
        write(*,*) ' STOP: Optimization not requested; terminating at',
     &             ' starting point.'
        write(*,*) ' '
       goto 500
      endif
!
      do while (kloop.le.kloopmax)
!
        rctrl(34)=drandu1(iseed) ! generate a new random number (used in some test problems)
!        
        outerloop=kloop
        lloop=0
        rho=dml_inf0
!
 110    continue
        do i=1,n
          drange=min(rangemax,rho*(x_u0(i)-x_l0(i)))
          x_l(i)=max(x(i)-drange,x_l0(i))
          x_u(i)=min(x(i)+drange,x_u0(i))
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!               possibly start an inner conservative loop              !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  determine the sparsity of the Jacobian
        call sparser (n,ni,gf,gc,ns,ictrl,lctrl,rctrl,cctrl,warn,
     &                iuser,luser,cuser,ruser)
     
!  initialize and check for unconstrained problems
        ksubiter=0
        times1 = SAOi_seconds()
!
        if (unconstrained) then 
!
          if (outerloop.eq.1) message = -50
!
          if (alg_unconstr.eq.2) then
!          
!  call the pure spherical quadratic approximation - not recommended ...             
            call uncstr (n,x,f_a,x_h,acurv,acurvn,x_l,x_u,f,gf,ksubiter,
     &                  iuser,luser,cuser,ruser,ictrl,lctrl,rctrl,cctrl)
!            
          elseif (alg_unconstr.eq.1) then
!          
!  call l-bfgs-b-2.4           
            call drive_lbfgsb24u (n,x,ni,ne,x_h,
     &                           acurv,bcurv,ccurv,acurvn,
     &                           bcurvn,ccurvn,x_l,x_u,
     &                           f,c,h,gf,gc,gh,ksubiter,
     &                           xkkt,flam,cstage0,ostring,
     &                           ns,ictrl,lctrl,rctrl,cctrl,message,yhi,
     &                           outerloop,eqn,lin,acol,nnz,do_diff,
     &                           nfe,nge,iuser,luser,cuser,ruser)
          endif ! unconstrained
!
        else    !   constrained
!        
         if (subsolver.eq.1) then
!         
!  call the l-bfgs-b-2.4 solver to solve the inequality / equality Falk dual (iFalk = 0 / 1)
            call drive_lbfgsb24f (n,x,ni,ne,f_a,c_a,h_a,iact,nact,xlam,
     &                           x_h,acurv,bcurv,ccurv,
     &                           acurvn,bcurvn,ccurvn,
     &                           x_l,x_u,f,c,h,gf,gc,gh,ksubiter,
     &                           xkktl,flam,cstage0,ostring,xlam_h,
     &                           ns,ictrl,lctrl,rctrl,cctrl,message,yhi,
     &                           outerloop,eqn,lin,
     &                           iuser,luser,cuser,ruser)
!
          elseif (subsolver.eq.21) then
!
!  call the cplex solver to solve the equivalent diagonal QP problem
            call drive_cplex (n,x,ni,ne,f_a,c_a,h_a,iact,nact,xlam,
     &                          x_h,acurv,bcurv,ccurv,
     &                          acurvn,bcurvn,ccurvn,
     &                          x_l,x_u,f,c,h,gf,gc,gh,ksubiter,
     &                          xkktl,flam,cstage0,ostring,xlam_h,
     &                          ns,ictrl,lctrl,rctrl,cctrl,message,yhi,
     &                          outerloop,lloop,eqn,
     &                          iuser,luser,cuser,ruser)      
!
          elseif (subsolver.ge.25.and.subsolver.le.28) then
!          
!  call the Galahad QP solvers to solve the equivalent diagonal QP problem
            call drive_galQP (n,x,ni,ne,f_a,c_a,h_a,xlam,x_h,
     &                        acurvn,bcurvn,ccurvn,x_l,x_u,f,
     &                        c,h,gf,gc,gh,ksubiter,
     &                        xkktl,flam,cstage0,ostring,xlam_h,
     &                        ns,ictrl,lctrl,rctrl,cctrl,message,yhi,
     &                        outerloop,(n*n+n)/2,ierr,eqn,lin,
     &                        iuser,luser,cuser,ruser,sloop,z)
!
          else
!          
            stop ' subsolver not defined '
!            
          endif    !   constrained
!
        endif

! catch floating point exceptions
        call testNaN(f_a, ierr)
        if (ierr.lt.0) then
          !accept = .false.
          !trust_region = .true.
          !goto 1233
          write(*,*) ' NaN occured in sao_dense '
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

!  determine if feasible
        if (unconstrained) then
          call SAOistatus(n,ni,iactblo,iactbhi,iactc,cplus_norm,xlam,
     &                   viol,c,x,x_l0,x_u0,5,-1,hsum,ictrl,lctrl,
     &                   rctrl,cctrl,feasible,iuser,luser,cuser,ruser)
        else
          call SAOistatus(n,ni,iactblo,iactbhi,iactc,cplus_norm,xlam,
     &                   viol,c,x,x_l0,x_u0,1,-1,hsum,ictrl,lctrl,
     &                   rctrl,cctrl,feasible,iuser,luser,cuser,ruser)
          call SAOistatusE(n,ni,iactblo,iactbhi,iactc,cplus_norm,xlam,
     &                    viol,c,x,x_l0,x_u0,1,1,hsum,ictrl,lctrl,
     &                    rctrl,cctrl,feasible,eqn,lin,
     &                    iuser,luser,cuser,ruser)
        endif

        accept=.true.

!  adjust the relaxation variable upper bound
      ymax = max(0.d0,viol)
      rctrl(25) = ymax

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
!
        if (trust_region.and.f.lt.f_h.and.unconstrained) then
          goto 1234                      ! do not enforce conservatism
        endif

!  escape if a feasible descent step was made
        if (conservative.and.f.lt.f_h.and.feasible.and.allow_f) then
          goto 1234                      ! do not enforce conservatism
        endif

!  escape if an infeasible restoration step was made
!         if (conservative.and.viol.lt.v_h.and.cstage0.and.allow_c) then
!           goto 1234                      ! do not enforce conservatism
!         endif                                                 ! agap-2008

!  determine if a conservative inner loop is required
        if (conservative.and.kloop.gt.-1) then
          call conserve1 (n,ni,accept,f,f_a,c,c_a,eps_locala,
     &                    eps_localb,acurv,acurvn,bcurv,
     &                    bcurvn,rho,ictrl,lctrl,rctrl,cctrl,lloop,
     &                    x_h,f_h,c_h,gf_h,gc_h,x,gf,gc,
     &                    iuser,luser,cuser,ruser)
        endif

!  determine if a trust region step is required
        if (trust_region) then
          call filter_ise (n,ni,accept,f,f_h,f_a,viol,v_h,
     &                     c,c_h,c_a,kloop,acurv,
     &                     acurvn,bcurv,bcurvn,
     &                     eps_locala,eps_localb,
     &                     filter_list,infilter,init,rho,
     &                     ictrl,lctrl,rctrl,cctrl,lloop,
     &                     x_h,gf_h,gc_h,x,gf,gc,
     &                     iuser,luser,cuser,ruser)
          if (rho.eq.rho_l_min) then
            goto 1234 ! no further reduction possible
          endif
        endif

!  increment inner loop counters
 1233   if (.not.accept.and.lloop.lt.lloopmax)  then
          lloop=lloop+1
          llooptot=llooptot+1

!  restore iterate
          call store_iter (n,ni,ne,f_h,f,v_h,viol,x_h,x_h2,x,c_h,c,
     &                     cplus_norm,cplus_norm_h,gf_h,gf,gc_h,gc,
     &                     h_h,h,gh_h,gh,xlam,xlam_h,-1,0,z,z_h,
     &                     iuser,luser,cuser,ruser,ictrl,lctrl,
     &                     rctrl,cctrl)
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
     &                     viol,c,x,x_l0,x_u0,1,1,hsum,ictrl,lctrl,
     &                     rctrl,cctrl,feasible,eqn,lin,
     &                     iuser,luser,cuser,ruser)
        endif

!  evaluate some filter thingies
        pred_f=f_h-f_a
        real_f=f_h-f
        if (dabs(pred_f).gt.1.d-8) then
          qindicator=real_f/pred_f
        else
          qindicator=0.d0
        endif

!  calculate relative step sizes
        delxnormc = norm2b(n,x_h,x)
        delxnormi = norm2f(n,x_h,x)
        if (delxnormc.lt.xtol*1.d1.or.delxnormi.lt.xtol_inf*1.d1)
     &       finalstg=.true.

!  calculate relative function difference
        rfd=(f_h-f)/(1.d0+dabs(f_h))

!  determine gradients
        timeg1 = SAOi_seconds()
        call formgrad (n,x,ni,ne,f,c,gf,gc,ictrl,lctrl,rctrl,cctrl,
     &                 nfe,nge,do_diff,iuser,luser,cuser,ruser,nnz,
     &                 Acol,ncol,eqn,lin)
        timeg2 = SAOi_seconds()
        nge=nge+1
        timeg=timeg+(timeg2-timeg1)

!  construct the approximate diagonal hessian
        call diaHess(n,x,ni,ne,x_h,x_h2,acurv,bcurv,ccurv,shift,
     &               acurvn,bcurvn,ccurvn,x_l,x_u,f,c,h,gf,gc,gh,
     &               gf_h,gc_h,gh_h,Li,Ui,Li_h,Ui_h,f_h,c_h,h_h,s,
     &               ictrl,lctrl,rctrl,cctrl,outerloop,iuser,luser,
     &               cuser,ruser)

!  optionally determine first and second order derivatives here - for debugging only... ! 
      if (GetFinDiffGradHessDense) then
        call FinDiffGradHessDense (n,x,ni,ni-ne,ne,f,c,ictrl,lctrl,
     &                             rctrl,cctrl,do_diff,iuser,luser,
     &                             cuser,ruser,eqn,lin,sloop) 
      endif

!  store iterate
        ig=0
        if (approx_c.eq.2.or.approx_c.eq.3.or.approx_c.eq.5
     &                                    .or.approx_c.eq.8) ig=1
        call store_iter (n,ni,ne,f_h,f,v_h,viol,x_h,x_h2,x,c_h,c,
     &                   cplus_norm,cplus_norm_h,gf_h,gf,gc_h,gc,
     &                   h_h,h,gh_h,gh,xlam,xlam_h,1,ig,z,z_h,
     &                   iuser,luser,cuser,ruser,ictrl,lctrl,
     &                   rctrl,cctrl)

!  calculate the active set on the subproblem level
        call form_act(ni,c,iact,nact,ictrl,lctrl,rctrl,cctrl,
     &                iuser,luser,cuser,ruser)

!  calculate the true KKT conditions
        if (calc_kkt.and.subsolver.ne.5) then
          call form_kkt (xkkt,n,ni,ne,x,iact,nact,x_l0,x_u0,
     &                   gf,gc,gh,xlam,xlamsml,xlambig,ifree,
     &                   ictrl,lctrl,rctrl,cctrl,z,
     &                   iuser,luser,cuser,ruser)
        elseif (calc_kkt.and.subsolver.eq.5) then
          xkkt=xkkt_mp
        else
          xkkt=-1.d0
        endif

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
!

!  again flush some buffers the fortran way...
        close(9)
        open (9 ,file='History.out',status='old',position='append')
        if (debug) then
          close(12)
          open (12,file='Subproblems.out',status='old',
     &          position='append')
        endif

!  possibly do specialized things
        call SAOi_Special (n,ni,ne,x,f,c,h,iuser,luser,cuser,ruser,
     &                     eqn,lin,ictrl,lctrl,rctrl,cctrl,x_l0,x_u0,
     &                     kloop,.false.,xkkt,viol,iactblo,iactbhi,
     &                     iactc,delxnormc,delxnormi,timey)

! increment outer loop counter
        kloop=kloop+1

! a heuristic
!        rctrl(17)=0.d0

!
      end do             ! while kloop.lt.kloopmax
 500  continue
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                   end of outer optimization loop                     !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  write a final summary
      do ipr=6,9,3
        if ((ipr.eq.6.and.iprint.ge.2).or.ipr.eq.9) then
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
!
      write (55,*) priv1,kloop,ksubitertot,f,viol

!  calculate the true KKT conditions
      if (calc_kkt) then
        call form_kkt (xkkt,n,ni,ne,x,iact,nact,x_l0,x_u0,
     &                 gf,gc,gh,xlam,xlamsml,xlambig,ifree,
     &                 ictrl,lctrl,rctrl,cctrl,z,
     &                 iuser,luser,cuser,ruser)
      elseif (calc_kkt.and.subsolver.eq.5) then
        xkkt=xkkt_mp
      else
        xkkt=-1.d0
      endif

!  possibly do specialized things
        call SAOi_Special (n,ni,ne,x,f,c,h,iuser,luser,cuser,ruser,
     &                     eqn,lin,ictrl,lctrl,rctrl,cctrl,x_l0,x_u0,
     &                     kloop,.true.,xkkt,viol,iactblo,iactbhi,
     &                     iactc,delxnormc,delxnormi,timey)

!  print a few things
 4000 do ipr=6,9,3
        if ((ipr.eq.6.and.iprint.ge.2).or.ipr.eq.9) then
          if (ni+ne.gt.0) then
            write(ipr,*)' '
            write(ipr,8900) xlamsml
            write(ipr,8901) xlambig
          endif
          write(ipr,*)' '
          write(ipr,8902) xkkt
!           if (fapriori.lt.big/10.d0) then
!             write(ipr,*)' '
!             write(ipr,8903) f-fapriori
!             if (ipr.eq.6) then
!               if (ne.ne.0) then
!                 write(47,*) mxloop,n,ni,ne
!               else
!                 write(47,*) mxloop,n,ni,ne,f-fapriori,
!      &                      xkkt,kloop,llooptot
!               endif
!             endif
!           endif
        endif
      enddo
!
      if (yhi.gt.1.d-8) then
        call open_warn_file (warn)
        write (14,8801)yhi
      endif
!
      if (iprint.ge.3.or.multimodal) then
        do i=1,n
          write(8,5001) i,x0(i),x(i),x_l0(i),x_u0(i)
        enddo
      endif
!
      if (iprint.ge.3.or.multimodal) then
        do j=1,ni
          write(13,8501) j,c(j),xlam(j)
          if (xlam(j).ge.biglam) then
            write(9,8600) j
            write(13,8601) j
            call open_warn_file (warn)
            write(14,8602) j,biglam
          endif
        enddo
      endif
!
      if (debug) then
        close(12)
        open (12 ,file='Subproblems.out',status='old',position='append')
      endif
!
        close(28)
!       open (28,position='append')
!       write(28,*) n,f,iactblo+iactbhi,xlam(1),xlam(2),kloop
!       close(28)
!
!
!      write(*,*) f -f_a
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
 6000 format (/,46x,'SAOi algorithm',//,130('-'),/,1x,
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
 8903 format (' f-f_optimal : ',1es18.8)
!
 9000 format (/,/,'  Customized output: ',1f12.2,2i7,1f8.3,10f12.4,/)
 9500 format ('    It              f              h       b&w    p-SIMP'
     &        ,'     q-GSS      brns',//,77('-'))
 9600 format (5x)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end subroutine sao_dense

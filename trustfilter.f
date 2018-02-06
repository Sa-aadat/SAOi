! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Enforce convergence using a trust region. Experimental
! SAOi: implementation; improve the coding
! SAOi:


      subroutine filter_ise (n,ni,accept,f,f_h,f_a,v,v_h,c,c_h,c_a,
     &                       kloop,acurv,acurvn,bcurv,
     &                       bcurvn,eps_locala,eps_localb,
     &                       filter_list,infilter,init,rho,
     &                       ictrl,lctrl,rctrl,cctrl,lloop,
     &                       x_h,gf_h,gc_h,x,gf,gc,
     &                       iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          kloop,infilter,i,n,ni,il,lloop
      double precision c(ni),c_h(ni),c_a(ni),filter_list(2,nmax)
      double precision acurv,acurvn(n),bcurv(ni),bcurvn(ni,n)
      double precision rho,beta,gamma,delta_f
      double precision sigma,fact_red,delta_q
      double precision test1left,test2left,v,f,v_h,f_h,gf(n)
      double precision test1right,test2right,f_a,eps_locala,eps_localb
      double precision x_h(n),gf_h(n),gc_h(ni,n),x(n),gc(ni,n)
      logical          accept,acceptable,update_filt,init,iprntl
!
      data             gamma                /1.00d-3/     !  /1.00d-5/
      data             sigma                /1.00d-1/     !  /1.00d-1/
      data             fact_red             /2.00d00/     !  /2.00d00/
      data             iprntl               /.false./     !  /.false./
!
      include          'ctrl_get.inc'
!      
      if (iprntl) write(*,*) ' entering trust-filter'
!      
      beta=1.d0-gamma
!
      if (kloop.eq.1.and.init) then  ! initialize the filter and check
        do i=1,outermax
          filter_list(1,i)=0.d0
          filter_list(2,i)=0.d0
       enddo
        if (sigma.lt.gamma) stop ' sigma < gamma in filter_ise'
        if (v.gt.beta*filter_hi) then
          write(*,*) ' '
          write(*,*) ' Initial filter behavior depends on filter_hi... '
          write(*,*) ' Press Enter to continue '
          read (*,*)
          ! TO-DO: automate reasonable filter upper bounds
        endif
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
          else
            call conserve1(n,ni,accept,f,f_a,c,c_a,eps_locala,
     &                     eps_localb,acurv,acurvn,bcurv,bcurvn,
     &                     rho,ictrl,lctrl,rctrl,cctrl,lloop,
     &                     x_h,f_h,c_h,gf_h,gc_h,x,gf,gc,
     &                     iuser,luser,cuser,ruser)
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
        if (iprntl) write(*,*)' case a: do nothing'
      else
        acceptable=.false.
        accept=.false.
        if (classic_trust) then
          rho=rho/fact_red
        else
          call conserve1(n,ni,accept,f,f_a,c,c_a,eps_locala,
     &                   eps_localb,acurv,acurvn,bcurv,bcurvn,
     &                   rho,ictrl,lctrl,rctrl,cctrl,lloop,
     &                   x_h,f_h,c_h,gf_h,gc_h,x,gf,gc,
     &                   iuser,luser,cuser,ruser)
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
          if (rho.le.xtol*10.d0) goto 1000
        else
          call conserve1(n,ni,accept,f,f_a,c,c_a,eps_locala,
     &                    eps_localb,acurv,acurvn,bcurv,bcurvn,
     &                    rho,ictrl,lctrl,rctrl,cctrl,lloop,
     &                    x_h,f_h,c_h,gf_h,gc_h,x,gf,gc,
     &                    iuser,luser,cuser,ruser)
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
      subroutine update_filter_ise (filter_list,infilter,v,f,beta,gamma,
     &                              ictrl,lctrl,rctrl,cctrl,
     &                              iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,infilter,infilter_del
      double precision filter_list(2,nmax),v,f,beta,gamma
      logical          delpoint
      include          'ctrl_get.inc'
!
      if (v.lt.0.d0) return

! expand the filter by 1, add point x^k <---> (h^k, f^k)
      infilter=infilter+1
      if (infilter.gt.outermax+1) then
        stop ' infilter > outermax+1 in conserve1.f'
      endif
      filter_list(1,infilter)=v
      filter_list(2,infilter)=f
      infilter_del=0
!
      delpoint=.false.
      do i=1,infilter-1
        if (filter_list(1,i) .ge. v .and. filter_list(2,i) .ge. f) then
          delpoint=.true.
          filter_list(1,i)=0.d0
          filter_list(2,i)=0.d0
          infilter=infilter-1
        endif
      enddo
!
      call bubsort3a(filter_list,2,nmax)
!
      return
      end
!----------------------------------------------------------------------!

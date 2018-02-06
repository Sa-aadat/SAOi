! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Approximate diagonal `Hessian' terms
! SAOi: To do: improve coding & blassify
! SAOi:


      subroutine diaHess(n,x,ni,ne,x_h,x_h2,acurv,bcurv,ccurv,shift,
     &                   acurvn,bcurvn,ccurvn,x_l,x_u,f,c,h,gf,gc,gh,
     &                   gf_h,gc_h,gh_h,Li,Ui,Li_h,Ui_h,f_h,c_h,h_h,s,
     &                   ictrl,lctrl,rctrl,cctrl,outerloop,iuser,luser,
     &                   cuser,ruser)
!
      implicit         none
      include          'ctrl.h'
      integer          n,ni,ne,i,j,outerloop
      double precision x(n),x_h(n),x_h2(n),s(n),shift(*)
      double precision x_l(n),x_u(n),acurv,bcurv(ni),ccurv(ne)
      double precision acurvn(n),bcurvn(ni,n),ccurvn(ne,n)
      double precision f,c(ni),h(ne),gf(n),gc(ni,n),gh(ne,n)
      double precision f_h,c_h(ni),h_h(ne)
      double precision gf_h(n),gc_h(ni,n),gh_h(ne,n)
      double precision Li(n),Ui(n),Li_h(n),Ui_h(n)
      include          'ctrl_get.inc'
!=================================================================================================!
!                                                                                                 !
!  This soubroutine allows for the selection of various diagonal                                  !
!    approximations. If exact DIAGONAL information is available,                                  !
!    this may be coded in diaHessUser by setting approx_f and/or                                  !
!    approx_c = 100. This is of course only possible if the problem                               !
!    is separable (the objective and/or all the constraints).                                     !
!                                                                                                 !
!=================================================================================================!
!
!
!
! notation: "T2:*" = "quadratic Taylor series expansion to approximation *"
!
!  0      : automated (incomplete)
!  1      : classical spherical quadratic
!  2      : non-classical spherical quadratic (error-norm)
!  3      : nabla2
!  4      : T2 to the reciprocal approximation R, denoted T2:R (with abs operator)
!  5      : T2 to the exponential approximation E, denoted T2:E (with abs operator, ai < 0)
!  6      : T2 to the MMA approximations, denoted T2:MMA
!  7      : T2 to the CONLIN approximations, denoted T2:CONLIN
!  8      : DQA of Park and Choi
! 14      : T2 to the reciprocal approximation R, denoted T2:R (nonconvex, without abs operator)
! >= 100  : user approximations
! else    : nothing was defined; terminal
!
!
! x    = x^(k)
! x_h  = x^(k-1)
! x_h2 = x^(k-2)
!

! check if the requested approximations have indeed been implemented
      if ((approx_f.lt.0.or.approx_f.gt.8).and.
     &     approx_f.ne.14.and.
     &     approx_f.ne.20.and.
     &     approx_f.ne.21.and.
     &     approx_f.ne.50.and.
     &     approx_f.ne.100) then
        stop 'approx_f not defined '
      endif
      if ((approx_c.lt.0.or.approx_c.gt.8).and.
     &     approx_c.ne.14.and.
     &     approx_c.ne.20.and.
     &     approx_c.ne.21.and.
     &     approx_c.ne.50.and.
     &     approx_c.ne.100) then
        stop 'approx_c not defined '
      endif

! best fit
      if     (approx_f.eq.0 .or. approx_c.eq.0) then
        stop ': approx_f .eq. 0 .and. approx_c .eq. 0 are incomplete'
      endif

! classical spherical quadratic
      if (approx_f.eq.1 .or. approx_c.eq.1) then
        call diaHess_sq1 (n,ni,ne,acurv,bcurv,ccurv,x,x_h,
     &                    gf,gc,gh,f,c,h,f_h,c_h,h_h,
     &                    acurvn,bcurvn,ccurvn,
     &                    ictrl,lctrl,rctrl,cctrl,outerloop,
     &                    iuser,luser,cuser,ruser)
      endif

! non-classical spherical quadratic
      if (approx_f.eq.2 .or. approx_c.eq.2) then
        call diaHess_sq2 (n,ni,ne,acurv,bcurv,ccurv,x,x_h,
     &                    gf,gc,gh,f,c,h,f_h,c_h,h_h,gf_h,gc_h,gh_h,
     &                    acurvn,bcurvn,ccurvn,
     &                    ictrl,lctrl,rctrl,cctrl,outerloop,
     &                    iuser,luser,cuser,ruser)
      endif

! nabla2
      if (approx_f.eq.3 .or. approx_c.eq.3) then
        call diaHess_n2 (n,ni,ne,acurvn,bcurvn,ccurvn,
     &                   x,x_h,gf,gc,gh,gf_h,gc_h,gh_h,
     &                   ictrl,lctrl,rctrl,cctrl,outerloop,
     &                   iuser,luser,cuser,ruser)
      endif

! quadratic Taylor series of the reciprocal approximation (with abs operator)
      if (approx_f.eq.4 .or. approx_c.eq.4) then
        call diaHess_t2r (n,ni,ne,acurvn,bcurvn,ccurvn,x,gf,gc,gh,
     &                    ictrl,lctrl,rctrl,cctrl,outerloop,
     &                    iuser,luser,cuser,ruser)
      endif

! quadratic Taylor series of the exponential approximation (with abs operator, ai < 0)
      if (approx_f.eq.5 .or. approx_c.eq.5) then
        call diaHess_t2e (n,ni,ne,x,x_h,gf,gc,gh,
     &                    gf_h,gc_h,gh_h,acurvn,bcurvn,ccurvn,
     &                    ictrl,lctrl,rctrl,cctrl,outerloop,
     &                    iuser,luser,cuser,ruser)
      endif

! quadratic Taylor series of the MMA approximations
      if (approx_f.eq.6 .or. approx_c.eq.6) then
        call diaHess_t2mma0 (n,ni,ne,x,x_h,gf,gc,gh,x_h2,Li,Ui,
     &                       x_u,x_l,Li_h,Ui_h,acurvn,bcurvn,ccurvn,
     &                       ictrl,lctrl,rctrl,cctrl,outerloop,
     &                       iuser,luser,cuser,ruser)
      endif

! quadratic Taylor series of the reciprocal approximation (with abs operator)
      if (approx_f.eq.7 .or. approx_c.eq.7) then
        call diaHess_t2conlin (n,ni,ne,acurvn,bcurvn,ccurvn,x,gf,gc,gh,
     &                         ictrl,lctrl,rctrl,cctrl,outerloop,
     &                         iuser,luser,cuser,ruser)
      endif

! quadratic approximation DQA of Park and Choi
      if (approx_f.eq.8 .or. approx_c.eq.8) then
         call diaHess_DQA (n,ni,ne,x,x_l,x_h,gf,gc,gh,f,c,h,f_h,c_h,
     &                     h_h,gf_h,gc_h,gh_h,acurvn,bcurvn,ccurvn,
     &                     shift,ictrl,lctrl,rctrl,cctrl,outerloop,
     &                     iuser,luser,cuser,ruser)
      endif

! quadratic Taylor series of the reciprocal approximation (without abs operator)
      if (approx_f.eq.14 .or. approx_c.eq.14) then
        call diaHess_t2r_nc (n,ni,ne,acurvn,bcurvn,ccurvn,x,gf,gc,gh,
     &                       ictrl,lctrl,rctrl,cctrl,outerloop,
     &                       iuser,luser,cuser,ruser)
      endif

! modified classical spherical quadratic
      if (approx_f.eq.20 .or. approx_c.eq.20) then
        call diaHess_sq1sprse (n,ni,ne,acurv,bcurv,ccurv,x,x_h,
     &                         gf,gc,gh,f,c,h,f_h,c_h,h_h,
     &                         acurvn,bcurvn,ccurvn,
     &                         ictrl,lctrl,rctrl,cctrl,outerloop,
     &                         iuser,luser,cuser,ruser)
      endif

! another non-classical spherical quadratic
      if (approx_f.eq.21 .or. approx_c.eq.21) then
        call diaHess_sq3 (n,ni,ne,acurv,bcurv,ccurv,x,x_h,x_h2,
     &                    gf,gc,gh,f,c,h,f_h,c_h,h_h,x_l,x_u,
     &                    acurvn,bcurvn,ccurvn,s,
     &                    ictrl,lctrl,rctrl,cctrl,outerloop,
     &                    iuser,luser,cuser,ruser)
      endif

! quasi-Cauchy updates
      if (approx_f.eq.50 .or. approx_c.eq.50) then
        call diaHess_qcauchy (n,ni,ne,x,x_h,gf,gc,gh,
     &                        gf_h,gc_h,gh_h,acurvn,bcurvn,ccurvn,
     &                        ictrl,lctrl,rctrl,cctrl,outerloop,
     &                        iuser,luser,cuser,ruser)
      endif

! some user approximations; many more may be added
      if (approx_f.eq.100 .or. approx_c.eq.100) then
        call diaHessUser (n,ni,ne,acurv,bcurv,ccurv,f,c,h,
     &                    x,x_h,x_h2,gf,gc,gh,gf_h,gc_h,gh_h,
     &                    acurvn,bcurvn,ccurvn,ictrl,lctrl,
     &                    rctrl,cctrl,outerloop,iuser,luser,
     &                    cuser,ruser,shift)
      endif
!
      return
      end subroutine diaHess
!----------------------------------------------------------------------c
      subroutine diaHess_sq1 (n,ni,ne,acurv,bcurv,ccurv,x,x_h,
     &                        gf,gc,gh,f,c,h,f_h,c_h,h_h,
     &                        acurvn,bcurvn,ccurvn,
     &                        ictrl,lctrl,rctrl,cctrl,outerloop,
     &                        iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,outerloop
      double precision acurv,bcurv(ni),ccurv(ne),x(n),x_h(n),delx(n)
      double precision acurvn(n),bcurvn(ni,n),ccurvn(ne,n)
      double precision gf(n),gc(ni,n),gh(ne,n),dp,zero,one,delxnorm2
      double precision f,c(ni),h(ne),f_h,c_h(ni),h_h(ne)
      data             zero /0.d0/
      data             one  /1.d0/
      include          'ctrl_get.inc'
!
! the classical spherical quadratic approximation
!
      if (outerloop.eq.0) then    ! initialize the approximation

        if (approx_f.eq.1) then
          do i=1,n
            acurvn(i) = one
          end do
        endif
!
        if (approx_c.eq.1) then
         do i=1,n
           do j=1,ni
              bcurvn(j,i) = one ! btol1 ! zero 
            end do
          end do
         do j=1,ni
            bcurv(j) = one ! btol1 ! zero 
          end do
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
          do i=1,n
            acurvn(i)=acurv
          enddo
        endif
        if (approx_c.eq.1) then
          do j=1,ni
            dp=0.d0
            do i=1,n
              dp=dp+gc(j,i)*delx(i) ! recode
            end do
            bcurv(j)=2.d0*(c_h(j)-c(j)-dp)/max(1.d-10,delxnorm2)
            bcurv(j)=max(btol1,bcurv(j))
          end do
          do i=1,n
            do j=1,ni
              bcurvn(j,i)=bcurv(j)
            enddo
          enddo
        endif
      endif
!
      return
      end subroutine diaHess_sq1
!----------------------------------------------------------------------c
      subroutine diaHess_sq3 (n,ni,ne,acurv,bcurv,ccurv,x,x_h,x_h2,
     &                        gf,gc,gh,f,c,h,f_h,c_h,h_h,x_l,x_u,
     &                        acurvn,bcurvn,ccurvn,s,
     &                        ictrl,lctrl,rctrl,cctrl,outerloop,
     &                        iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,outerloop
      double precision acurv,bcurv(ni),ccurv(ne),x(n),x_h(n),delx(n)
      double precision acurvn(n),bcurvn(ni,n),ccurvn(ne,n)
      double precision x_h2(n),s(n),x_l(n),x_u(n)
      double precision gf(n),gc(ni,n),gh(ne,n),dp,zero,one,delxnorm2
      double precision f,c(ni),h(ne),f_h,c_h(ni),h_h(ne)
      data             zero /0.d0/
      data             one  /1.d0/
      include          'ctrl_get.inc'
!
! another non-classical spherical quadratic approximation
!
      call smma (n,ni,ne,x,x_h,x_h2,x_u,x_l,s,
     &           ictrl,lctrl,rctrl,cctrl,outerloop,
     &           iuser,luser,cuser,ruser)
!
      if (outerloop.eq.0) then    ! initialize the approximation

        if (approx_f.eq.21) then
          do i=1,n
            acurvn(i) = one
          end do
        endif
!
        if (approx_c.eq.21) then
         do i=1,n
           do j=1,ni
              bcurvn(j,i) = one ! btol1
            end do
          end do
        endif
!
      else     ! the iteration number > 0
!
        delxnorm2=0.d0
        do i=1,n
          delx(i)=x_h(i)-x(i)
          delxnorm2=delxnorm2+(delx(i)/s(i))**2
        end do
        if (approx_f.eq.21) then
          dp=0.d0
          do i=1,n
            dp=dp+gf(i)*delx(i)
          end do
          acurv=2.d0*(f_h-f-dp)/max(1.d-10,delxnorm2)
          do i=1,n
            acurvn(i)=max(atol1,acurv/s(i)**2)
          enddo
        endif
        if (approx_c.eq.21) then
          do j=1,ni
            dp=0.d0
            do i=1,n
              dp=dp+gc(j,i)*delx(i) ! recode
            end do
            bcurv(j)=2.d0*(c_h(j)-c(j)-dp)/max(1.d-10,delxnorm2)
          end do
          do i=1,n
            do j=1,ni
              bcurvn(j,i)=max(btol1,bcurv(j)/s(i)**2)
            enddo
          enddo
        endif
      endif
!
      return
      end subroutine diaHess_sq3
!----------------------------------------------------------------------c
      subroutine diaHess_sq1sprse (n,ni,ne,acurv,bcurv,ccurv,x,x_h,
     &                             gf,gc,gh,f,c,h,f_h,c_h,h_h,
     &                             acurvn,bcurvn,ccurvn,
     &                             ictrl,lctrl,rctrl,cctrl,outerloop,
     &                             iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,outerloop
      double precision acurv,bcurv(ni),ccurv(ne),x(n),x_h(n),delx(n)
      double precision acurvn(n),bcurvn(ni,n),ccurvn(ne,n)
      double precision gf(n),gc(ni,n),gh(ne,n),dp,zero,one,delxnorm2
      double precision f,c(ni),h(ne),f_h,c_h(ni),h_h(ne),delxnormj
      data             zero /0.d0/
      data             one  /1.d0/
      include          'ctrl_get.inc'
!
! the classical spherical quadratic approximation
!
      if (outerloop.eq.0) then    ! initialize the approximation

        if (approx_f.eq.20) then
          do i=1,n
            acurvn(i) = one
          end do
        endif
!
        if (approx_c.eq.20) then
         do i=1,n
           do j=1,ni
              if (gc(j,i).eq.0.d0) then
                bcurvn(j,i) = 0.d0 ! btol1
              else
                bcurvn(j,i) = one ! btol1
              endif
            end do
          end do
        endif
!
      else     ! the iteration number > 0
!
        delxnorm2=0.d0
        do i=1,n
          delx(i)=x_h(i)-x(i)
          delxnorm2=delxnorm2+delx(i)**2
        end do
        if (approx_f.eq.20) then
          dp=0.d0
          do i=1,n
            dp=dp+gf(i)*delx(i)
          end do
          acurv=2.d0*(f_h-f-dp)/max(1.d-10,delxnorm2)
          do i=1,n
            acurvn(i)=max(atol1,acurv)
          enddo
        endif
!
        if (approx_c.eq.20) then
          do j=1,ni
            delxnormj=0.d0
            dp=0.d0
            do i=1,n
              if (gc(j,i).eq.0.d0) then
                ! do nothing
              else
                delx(i)=x_h(i)-x(i)
                delxnormj=delxnormj+delx(i)**2
                dp=dp+gc(j,i)*delx(i)             ! recode
              endif
            enddo
            bcurv(j)=2.d0*(c_h(j)-c(j)-dp)/max(1.d-10,delxnormj)
            do i=1,n
              if (gc(j,i).eq.0.d0) then
                bcurvn(j,i)=0.d0
              else
                bcurvn(j,i)=max(btol1,bcurv(j))   ! no need to store n components of bcurv !
              endif
            enddo
!
          end do       ! for j=1,ni
        endif          ! if approx_c = 20
      endif            ! if outerloop > 0
!
      return
      end subroutine diaHess_sq1sprse
!----------------------------------------------------------------------c
      subroutine diaHess_sq2 (n,ni,ne,acurv,bcurv,ccurv,x,x_h,
     &                        gf,gc,gh,f,c,h,f_h,c_h,h_h,gf_h,gc_h,gh_h,
     &                        acurvn,bcurvn,ccurvn,
     &                        ictrl,lctrl,rctrl,cctrl,outerloop,
     &                        iuser,luser,cuser,ruser)

      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,outerloop
      double precision acurv,bcurv(ni),ccurv(ne),x(n),x_h(n),delx(n)
      double precision acurvn(n),bcurvn(ni,n),ccurvn(ne,n)
      double precision gf(n),gc(ni,n),gh(ne,n),dp,zero,one,delxnorm2
      double precision gf_h(n),gc_h(ni,n),gh_h(ne,n)
      double precision f,c(ni),h(ne),f_h,c_h(ni),h_h(ne)
      data             zero /0.d0/
      data             one  /1.d0/
      include          'ctrl_get.inc'
!
! the non-classical spherical quadratic approximation
!
      if (outerloop.eq.0) then    ! initialize the approximation
!
        if (approx_f.eq.2) then
          do i=1,n
            acurvn(i) = one
          end do
        endif
        if (approx_c.eq.2) then
          do i=1,n
            do j=1,ni
              bcurvn(j,i) = one ! btol1
            end do
          end do
        endif
!
      else     ! the iteration number > 0
!
        delxnorm2=0.d0
        do i=1,n
          delx(i)=x_h(i)-x(i)
          delxnorm2=delxnorm2+delx(i)**2
        end do
        if (approx_f.eq.2) then
          dp=0.d0
          do i=1,n
            dp=dp+(gf_h(i)-gf(i))*delx(i)
          end do
          acurv=dp/max(1.d-10,delxnorm2)
          do i=1,n
            acurvn(i)=max(atol1,acurv)
          enddo
        endif
        if (approx_c.eq.2) then
          do j=1,ni
            dp=0.d0
            do i=1,n
              dp=dp+(gc_h(j,i)-gc(j,i))*delx(i) ! recode
            end do
            bcurv(j)=dp/max(1.d-10,delxnorm2)
          end do
          do i=1,n
            do j=1,ni
              bcurvn(j,i)=max(btol1,bcurv(j))
            enddo
          enddo
        endif
      endif
!
      return
      end subroutine diaHess_sq2
!----------------------------------------------------------------------c
      subroutine diaHess_n2 (n,ni,ne,acurvn,bcurvn,ccurvn,
     &                       x,x_h,gf,gc,gh,gf_h,gc_h,gh_h,
     &                       ictrl,lctrl,rctrl,cctrl,outerloop,
     &                       iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,outerloop
      double precision acurvn(n),bcurvn(ni,n),ccurvn(ne,n)
      double precision x(n),x_h(n),delx(n),zero,one,xtollocal
      double precision gf(n),gc(ni,n),gh(ne,n),fcurvres,ccurvres
      double precision gf_h(n),gc_h(ni,n),gh_h(ne,n)
      data             zero        /0.d00/
      data             one         /1.d00/
      data             xtollocal   /1.d-6/
      data             fcurvres    /1.d-3/
      data             ccurvres    /1.d-6/
      include          'ctrl_get.inc'
!
!   nonspherical quadratic based on gradients and direct variables
!
      if (outerloop.eq.0) then    ! initialize the approximation
!
        if (approx_f.eq.3) then
          do i=1,n
            acurvn(i) = one
          end do
        endif
        if (approx_c.eq.3) then
          do i=1,n
            do j=1,ni
              bcurvn(j,i) = btol1
            end do
          end do
        endif
!
      else     ! the iteration number > 0
!
        do i=1,n
          delx(i)=x_h(i)-x(i)
        end do
        if (approx_f.eq.3) then
          do i=1,n
            if (dabs(delx(i)).gt.xtollocal) then
              acurvn(i)=max(atol1,(gf_h(i)-gf(i))/delx(i))
            else
              acurvn(i)=max(atol1,fcurvres)
            endif
          enddo
        endif
        if (approx_c.eq.3) then
          do i=1,n
            do j=1,ni
              if (dabs(delx(i)).gt.xtollocal) then
                bcurvn(j,i)=max(btol1,(gc_h(j,i)-gc(j,i))/delx(i),
     &                          ccurvres)
              else
                bcurvn(j,i)=max(btol1,ccurvres)
              endif
            enddo
          enddo
        endif
      endif
!
      return
      end subroutine diaHess_n2
!-----------------------------------------------------------------------
      subroutine diaHess_t2r (n,ni,ne,acurvn,bcurvn,ccurvn,x,gf,gc,gh,
     &                        ictrl,lctrl,rctrl,cctrl,outerloop,
     &                        iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,outerloop
      double precision x(n),acurvn(n),bcurvn(ni,n),ccurvn(ne,n)
      double precision gf(n),gc(ni,n),gh(ne,n),zero,one
      data             zero /0.d0/
      data             one  /1.d0/
      include          'ctrl_get.inc'
!
! T2:R (may be constructed even when outerloop = 0)
!
      if (approx_f.eq.4) then
        do i=1,n
          acurvn(i) = max(atol1,2.d0/x(i)*dabs(gf(i)))
        enddo
      endif
      if (approx_c.eq.4) then
        call dgerx1(ni,n,2.d0,x,1,gc,1,bcurvn,ni,btol1)
      endif
!
      return
      end subroutine diaHess_t2r
!----------------------------------------------------------------------c
      subroutine diaHess_t2e (n,ni,ne,x,x_h,gf,gc,gh,
     &                        gf_h,gc_h,gh_h,acurvn,bcurvn,ccurvn,
     &                        ictrl,lctrl,rctrl,cctrl,outerloop,
     &                        iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,outerloop
      double precision ap(n),bp(ni,n),cp(ne,n)
      double precision acurvn(n),bcurvn(ni,n),ccurvn(ne,n)
      double precision x(n),x_h(n),delx(n),zero,one,exphi,explo
      double precision gf(n),gc(ni,n),gh(ne,n),exptol,a,b,anum
      double precision gf_h(n),gc_h(ni,n),gh_h(ne,n)
      data             zero   /0.d0/
      data             one    /1.d0/
      data             exphi  /-0.1d0/
      data             explo  /-4.d0/
      data             exptol /1.d-8/
      include          'ctrl_get.inc'
!
!  T2:E
!
      if (outerloop.eq.0) then    ! initialize the approximation
!
        if (approx_f.eq.5) then
          do i=1,n
            acurvn(i) = max(atol1,2.d0/x(i)*dabs(gf(i)))
          enddo
        endif
        if (approx_c.eq.5) then
          call dgerx1(ni,n,2.d0,x,1,gc,1,bcurvn,ni,btol1)
        endif
!
      else     ! the iteration number > 0
!
        if (approx_f.eq.5) then
          do i=1,n
            anum =-one
            b=x_h(i)/x(i)
            if (b.ne.1.d0.and.gf(i).ne.0.d0.and.
     &                                  gf_h(i).ne.0.d0) then
              a = gf_h(i)/gf(i)
              anum = one + dlog(a)/dlog(b)
            endif
            ap(i)= max(explo,min(anum,exphi))
          enddo
          do i = 1,n
            acurvn(i) = -(ap(i)-1.d0)/x(i)*dabs(gf(i))
            acurvn(i) = max(atol1,acurvn(i))
          enddo
        endif
!
        if (approx_c.eq.5) then
          do i=1,n
            do j=1,ni
              anum =-one
              b=x_h(i)/x(i)
              if (b.ne.1.d0.and.gc(j,i).ne.0.d0.and.
     &                                  gc_h(j,i).ne.0.d0) then
                a = gc_h(j,i)/gc(j,i)
                anum = one + dlog(a)/dlog(b)
              endif
              bp(j,i) = max(explo,min(anum,exphi))
            enddo
          enddo
          call dgerx2(ni,n,1.d0,x,1,gc,1,bcurvn,ni,bp,btol1)
        endif
      endif
!
      return
      end subroutine diaHess_t2e
!----------------------------------------------------------------------c
      subroutine diaHess_t2mma0 (n,ni,ne,x,x_h,gf,gc,gh,x_h2,Li,Ui,
     &                           x_u,x_l,Li_h,Ui_h,acurvn,bcurvn,ccurvn,
     &                           ictrl,lctrl,rctrl,cctrl,outerloop,
     &                           iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,outerloop
      double precision ap(n),bp(ni,n),cp(ne,n)
      double precision acurvn(n),bcurvn(ni,n),ccurvn(ne,n)
      double precision x(n),x_h(n),x_h2(n),delx(n),x_u(n),x_l(n)
      double precision gf(n),gc(ni,n),gh(ne,n),gammai
      double precision Li(n),Ui(n),Li_h(n),Ui_h(n),diff,osc
      double precision zero,one,two
      data             zero /0.d0/
      data             one  /1.d0/
      data             two  /2.d0/
      include          'ctrl_get.inc'
!
!  T2:MMA
!
      do i = 1,n ! recode
        if (outerloop.lt.2) then
          Li(i) = x(i)-(x_u(i)-x_l(i))/two
          Ui(i) = x(i)+(x_u(i)-x_l(i))/two
        else
          osc = (x(i)-x_h(i))*(x_h(i)-x_h2(i))
          if (osc.lt.zero) then
            gammai = mma_lo
          elseif (osc.gt.zero) then
            gammai = mma_hi
          else
            gammai = one
          endif
          Li(i) = x(i)-gammai*(x_h(i)-Li_h(i))
          Ui(i) = x(i)+gammai*(Ui_h(i)-x_h(i))
        endif
        Li_h(i) = Li(i)
        Ui_h(i) = Ui(i)
        if (approx_f.eq.6) then
          if (gf(i).le.0.d0) then
            acurvn(i) = -2.d0*gf(i)/(x(i)-Li(i))
          else
            acurvn(i) =  2.d0*gf(i)/(Ui(i)-x(i))
          endif
          acurvn(i) = max(atol1,acurvn(i))
        endif
        if (approx_c.eq.6) then
          do j=1,ni
            if (gc(j,i).le.0.d0) then
              bcurvn(j,i) = -2.d0*gc(j,i)/(x(i)-Li(i))
            else
              bcurvn(j,i) =  2.d0*gc(j,i)/(Ui(i)-x(i))
            endif
            bcurvn(j,i) = max(btol1,bcurvn(j,i))
          enddo
        endif
      enddo
!
      return
      end subroutine diaHess_t2mma0
!----------------------------------------------------------------------c
      subroutine diaHess_t2mma (n,ni,ne,x,x_h,gf,gc,gh,x_h2,Li,Ui,
     &                          x_u,x_l,Li_h,Ui_h,acurvn,bcurvn,ccurvn,
     &                          s,ictrl,lctrl,rctrl,cctrl,outerloop,
     &                          iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,outerloop
      double precision ap(n),bp(ni,n),cp(ne,n)
      double precision acurvn(n),bcurvn(ni,n),ccurvn(ne,n)
      double precision x(n),x_h(n),x_h2(n),delx(n),x_u(n),x_l(n),s(n)
      double precision gf(n),gc(ni,n),gh(ne,n),gammai
      double precision Li(n),Ui(n),Li_h(n),Ui_h(n),diff,osc
      double precision zero,one,two
      data             zero /0.d0/
      data             one  /1.d0/
      data             two  /2.d0/
      include          'ctrl_get.inc'
!
!  T2:MMA
!
      call smma (n,ni,ne,x,x_h,x_h2,x_u,x_l,s,
     &           ictrl,lctrl,rctrl,cctrl,outerloop,
     &           iuser,luser,cuser,ruser)
!
      do i = 1,n
        Li(i) = x(i) - s(i)
        Ui(i) = x(i) + s(i)
        if (approx_f.eq.6) then
          if (gf(i).le.0.d0) then
            acurvn(i) = -2.d0*gf(i)/(x(i)-Li(i))
          else
            acurvn(i) =  2.d0*gf(i)/(Ui(i)-x(i))
          endif
          acurvn(i) = max(atol1,acurvn(i))
        endif
        if (approx_c.eq.6) then
          do j=1,ni
            if (gc(j,i).le.0.d0) then
              bcurvn(j,i) = -2.d0*gc(j,i)/(x(i)-Li(i))
            else
              bcurvn(j,i) =  2.d0*gc(j,i)/(Ui(i)-x(i))
            endif
            bcurvn(j,i) = max(btol1,bcurvn(j,i))
          enddo
        endif
      enddo
!
      return
      end subroutine diaHess_t2mma
!----------------------------------------------------------------------c
      subroutine diaHess_t2conlin (n,ni,ne,acurvn,bcurvn,ccurvn,x,gf,
     &                        gc,gh,ictrl,lctrl,rctrl,cctrl,outerloop,
     &                        iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,outerloop
      double precision x(n),acurvn(n),bcurvn(ni,n),ccurvn(ne,n)
      double precision gf(n),gc(ni,n),gh(ne,n),zero,one
      data             zero /0.d0/
      data             one  /1.d0/
      include          'ctrl_get.inc'
!
! T2:CONLIN (may be constructed even when outerloop = 0)
!
      if (approx_f.eq.7) then
        do i = 1,n
          if (gf(i).lt.0.d0) then
            acurvn(i) = -2.d0/x(i)*gf(i)
          else
            acurvn(i) =  0.d0 ! 2.d0/x(i)*gf(i)
         endif
          acurvn(i) = max(atol1,acurvn(i))
        enddo
      endif
      if (approx_c.eq.7) then
        do j = 1,ni
          do i = 1,n
            if (gc(j,i).lt.0.d0) then
              bcurvn(j,i) = -2.d0/x(i)*gc(j,i) 
            else
              bcurvn(j,i) =  0.d0 ! 2.d0/x(i)*gc(j,i)
            endif
            bcurvn(j,i) = max(btol1,bcurvn(j,i))
          enddo
        enddo
      endif
!
      return
      end subroutine diaHess_t2conlin
!-----------------------------------------------------------------------
      subroutine diaHess_t2r_nc (n,ni,ne,acurvn,bcurvn,ccurvn,x,gf,gc, 
     &                           gh,ictrl,lctrl,rctrl,cctrl,outerloop,
     &                           iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,outerloop
      double precision x(n),acurvn(n),bcurvn(ni,n),ccurvn(ne,n)
      double precision gf(n),gc(ni,n),gh(ne,n),zero,one
      data             zero /0.d0/
      data             one  /1.d0/
      include          'ctrl_get.inc'
!
! T2:R (may be constructed even when outerloop = 0)
!
      if (approx_f.eq.14) then
        do i=1,n
          acurvn(i) = max(atol1,-2.d0/x(i)*gf(i))
        enddo
      endif
      if (approx_c.eq.14) then
        call dgerx1_nc(ni,n,2.d0,x,1,gc,1,bcurvn,ni,btol1)
      endif
!
      return
      end subroutine diaHess_t2r_nc
!----------------------------------------------------------------------c
      subroutine diaHess_DQA (n,ni,ne,x,x_l,x_h,gf,gc,gh,f,c,h,f_h,c_h,
     &                        h_h,gf_h,gc_h,gh_h,acurvn,bcurvn,ccurvn,
     &                        shift,ictrl,lctrl,rctrl,cctrl,outerloop,
     &                        iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,outerloop
      double precision ap(n),bp(ni,n),cp(ne,n),f,c(ni),h(ne),shift(n)
      double precision acurvn(n),bcurvn(ni,n),ccurvn(ne,n),c_h(ni),f_h
      double precision x(n),x_h(n),delx(n),zero,one,exphi,explo,h_h(ne)
      double precision gf(n),gc(ni,n),gh(ne,n),exptol,a,b,anum,dsum
      double precision gf_h(n),gc_h(ni,n),gh_h(ne,n),Gfi(n),Hfi(n),two
      double precision eta_f,eta_c(ni),deltaf,Cf(n),Cc(ni,n),deltac(ni)
      double precision Gci(ni,n),Hci(ni,n),yc(ni,n),yc_h(ni,n),three
      double precision yf(n),yf_h(n),temp1,temp2,temp3(ni),temp4(ni)
      double precision x_l(n)
      data             zero   /0.d0/
      data             one    /1.d0/
      data             two    /2.d0/
      data             three  /3.d0/
      data             exphi  /4.d0/
      data             explo  /-2.d0/
      data             exptol /1.d-8/
      include          'ctrl_get.inc'
!
!  DQA of Park and Choi
!  Coded by Seonho Park, 
!  errors or comments, contact to : Park.Seonho@gmail.com
!

! Dear Seanho, I have moved the code below to defaults.f.
!      You also have the option to specify shift_tol in 
!      Initialize.f (the default is 1.d-3)   

!
!       do i=1,n
!        if (x_l(i).lt.1.d-3) then
!          shift(i)=dabs(x_l(i))+1.d-3
!        else
!         shift(i)=0.d0
!        endif
!       enddo
!
      if (outerloop.eq.0) then    ! initialize the approximation
      !
        if (approx_f.eq.8) then
          do i=1,n
            if (gf(i).ge.0.d0) then
            ap(i) = one+exptol
            else
            ap(i) = -one
            endif
            acurvn(i) = (ap(i)-1.d0)/(x(i)+shift(i))*gf(i)
          enddo
        endif
        if (approx_c.eq.8) then
          do j=1,ni
            do i=1,n
                if (gc(j,i).ge.0.d0) then
                    bp(j,i) = one+exptol
                else
                    bp(j,i) = -one
                endif
            enddo
          enddo
          call dgerx3(ni,n,1.d0,x,1,gc,1,bcurvn,ni,bp)
        endif
!        
      else     ! the iteration number > 0
!
        if (approx_f.eq.8) then
          do i=1,n
            if (gf_h(i)*gf(i).le.0.d0) then
                anum=one
            else
                 a = gf_h(i)/gf(i)
                 b = (x_h(i)+shift(i))/(x(i)+shift(i))
                if (gf(i).lt.0.d0.and.(a-1.d0)*(b-1.d0).ge.0.d0) then
                      anum = -one
                elseif (gf(i).ge.0.d0.
     &          and.(a-1.d0)*(b-1.d0).le.0.d0) then
                    anum = three
                else
                if (gf(i).ge.0.d0) then
                 anum = min(exphi,max(two,one + dlog(a)/dlog(b)))
                else
                 anum = min(zero,max(explo,one + dlog(a)/dlog(b)))
                endif
                endif
            endif
            ap(i)=anum
            yf(i) = (x(i)+shift(i))**ap(i)
            yf_h(i) = (x_h(i)+shift(i))**ap(i)
            if (dabs(yf_h(i)-yf(i)).lt.exptol) then
                Gfi(i) = zero
            else
                Gfi(i) = (1.d0/(yf_h(i)-yf(i)))*(gf_h(i)
     &          *(x_h(i)+shift(i))**(1.d0-ap(i))/ap(i)-gf(i)
     &          *(x(i)+shift(i))**(1.d0-ap(i))/ap(i))
            endif
            if (gf_h(i)*gf(i).le.0.d0) then
            Hfi(i) = Gfi(i)
            else
            Hfi(i) = one
            endif
          enddo
          temp1 = zero
          temp2 = zero
          deltaf = zero
          do i=1,n
          temp1 = temp1+gf(i)*(x(i)+shift(i))**(1.d0-ap(i))
     &    /ap(i)*(yf_h(i)-yf(i))
          temp2 = temp2+Gfi(i)*(yf_h(i)-yf(i))**2.d0
          deltaf = deltaf+Hfi(i)*(yf_h(i)-yf(i))**2.d0
          enddo
          eta_f = 2*(f_h-f-temp1-0.5d0*temp2)
          do i=1,n
          if (deltaf.ne.0.d0) then
          Cf(i) = (eta_f*Hfi(i)/deltaf+
     &     Gfi(i))*(ap(i)*((x(i)+shift(i))**(ap(i)-1.d0)))**2        
          else
          Cf(i) = (eta_f*Hfi(i)/exptol+
     &    Gfi(i))*(ap(i)*((x(i)+shift(i))**(ap(i)-1.d0)))**2         
          endif
          Cf(i) = max(atol1,Cf(i))
          !Cf(i) = 0
          acurvn(i) = (ap(i)-1.d0)/(x(i)+shift(i))*gf(i)+Cf(i)
 !         acurvn(i) = max(atol1,acurvn(i))
          enddo
        endif
!
        if (approx_c.eq.8) then
          do i=1,n
            do j=1,ni
                if (gc_h(j,i)*gc(j,i).le.0.d0) then
                anum = one
                else
                 a = gc_h(j,i)/gc(j,i)
                 b = (x_h(i)+shift(i))/(x(i)+shift(i))
                 if (gc(j,i).lt.0.d0.and.(a-1.d0)*(b-1.d0).ge.0.d0) then
                    anum = -one
                 elseif (gc(j,i).ge.0.d0.and.
     &              (a-1.d0)*(b-1.d0).le.0.d0) then
                    anum = three
                 else
                 if (gc(j,i).ge.0.d0) then
                    anum = min(exphi,max(two,one + dlog(a)/dlog(b)))
                 else
                    anum = min(zero,max(explo,one + dlog(a)/dlog(b)))
                 endif
                 endif
                endif
                bp(j,i) = anum
                yc(j,i) = (x(i)+shift(i))**bp(j,i)
                yc_h(j,i) = (x_h(i)+shift(i))**bp(j,i)
                if (dabs(yc_h(j,i)-yc(j,i)).lt.exptol) then
                    Gci(j,i) = zero
                else
                  Gci(j,i)= (1.d0/(yc_h(j,i)-yc(j,i)))*(gc_h(j,i)
     &               *(x_h(i)+shift(i))**(1.d0-bp(j,i))/bp(j,i)-gc(j,i)
     &               *(x(i)+shift(i))**(1.d0-bp(j,i))/bp(j,i))
                endif
                if (gc_h(j,i)*gc(j,i).le.0.d0) then
                Hci(j,i) = Gci(j,i)
                else
                Hci(j,i) = one
                endif
             enddo
          enddo
          temp3 = zero
          temp4 = zero
          deltac = zero
          do j=1,ni
            do i=1,n
             temp3(j) = temp3(j)+gc(j,i)*(x(i)+shift(i))**(1.d0-bp(j,i))
     &        /bp(j,i)*(yc_h(j,i)-yc(j,i))
             temp4(j) = temp4(j)+Gci(j,i)*(yc_h(j,i)-yc(j,i))**2.d0
             deltac(j) = deltac(j)+Hci(j,i)*(yc_h(j,i)-yc(j,i))**2.d0  
            enddo
            eta_c(j) = 2*(c_h(j)-c(j)-temp3(j)-0.5d0*temp4(j))
          enddo
          do j=1,ni
            do i=1,n
             if (deltac(j).ne.0.d0) then
                Cc(j,i) = (eta_c(j)*Hci(j,i)/deltac(j)+
     &           Gci(j,i))*(bp(j,i)*((x(i)+shift(i))
     &           **(bp(j,i)-1.d0)))**2
             else
                Cc(j,i) = (eta_c(j)*Hci(j,i)/exptol+
     &           Gci(j,i))*(bp(j,i)*((x(i)+shift(i))
     &           **(bp(j,i)-1.d0)))**2
             endif
                Cc(j,i) = max(btol1,Cc(j,i))
                !Cc(j,i) = 0
            enddo
          enddo
          call dgerx3(ni,n,1.d0,x,1,gc,1,bcurvn,ni,bp)
          do j=1,ni
            do i=1,n
             bcurvn(j,i)=bcurvn(j,i)+Cc(j,i)
   !          bcurvn(j,i) = max(zero,bcurvn(j,i))
            enddo
          enddo
        endif
      endif
!
      return
      end subroutine diaHess_DQA
!----------------------------------------------------------------------c
      subroutine diaHess_qcauchy (n,ni,ne,x,x_h,gf,gc,gh,
     &                            gf_h,gc_h,gh_h,acurvn,bcurvn,ccurvn,
     &                            ictrl,lctrl,rctrl,cctrl,outerloop,
     &                            iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,ii,outerloop
      double precision ap(n),bp(ni,n),cp(ne,n)
      double precision acurvn(n),bcurvn(ni,n),ccurvn(ne,n)
      double precision x(n),x_h(n),delx(n),zero,one,exphi,explo
      double precision gf(n),gc(ni,n),gh(ne,n),exptol,a,b,anum
      double precision gf_h(n),gc_h(ni,n),gh_h(ne,n),s(n),y(n)
      double precision gcv(n),gcv_h(n),bcv(n)
      data             zero   /0.d0/
      data             one    /1.d0/
      data ii / -1/
      save ii
      include          'ctrl_get.inc'

      write(*,*) ' quasi-cauchy'

!
!      if (ii.lt.0) then
!        write(*,*) ' Enter ii '
!        read(*,*) ii
         ii=1
!      endif

!  quazi-Cauchy
      do i=1,n
        s(i)=x(i)- x_h(i)
        y(i)=gf_h(i)
      enddo
!
      if (approx_f.eq.50) then
        if (ii.eq.0) then
          call cauchy0 (outerloop+1,x,s,acurvn,gf,gf_h,y,n,
     &                  iuser,luser,cuser,ruser)
        elseif (ii.eq.1.or.ii.eq.2) then
          call qcauchy (outerloop+1,ii,x,s,acurvn,gf,gf_h,y,n,
     &                  iuser,luser,cuser,ruser)
        else
          stop ' ii not defined'
        endif
      endif
!
      if (approx_c.eq.50) then
        do j=1,ni
          do i=1,n
            gcv(i)   = gc(j,i)
            gcv_h(i) = gc_h(j,i)
            y(i)     = gc_h(j,i)
            bcv(i)   = bcurvn(j,i)
          enddo
          if (ii.eq.0) then
            call cauchy0 (outerloop+1,x,s,bcv,gcv,gcv_h,y,n,
     &                    iuser,luser,cuser,ruser)
          elseif (ii.eq.1.or.ii.eq.2) then
            call qcauchy (outerloop+1,ii,x,s,bcv,gcv,gcv_h,y,n,
     &                    iuser,luser,cuser,ruser)
            do i=1,n
              bcurvn(j,i) = bcv(i)
            enddo
          else
            stop ' ii not defined'
          endif
        enddo
      endif
!
      return
      end subroutine diaHess_qcauchy
!----------------------------------------------------------------------c

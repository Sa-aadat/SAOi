! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Enforce convergence using conservatism. Note to the developers:
! SAOi: gradients lag functions by one, so f_h and gf are complimentary,
! SAOi: etc.


      subroutine conserve1 (n,ni,accept,f,f_a,c,c_a,eps_locala,
     &                      eps_localb,acurv,acurvn,bcurv,bcurvn,rho,
     &                      ictrl,lctrl,rctrl,cctrl,lloop,
     &                      x_h,f_h,c_h,gf_h,gc_h,x,gf,gc,
     &                      iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,n,ni,lloop
      double precision bcurv(ni),bcurvn(ni,n),fact_red,fact_inc,factx
      double precision rho,c(ni),c_a(ni),acurv,acurvn(n),zero,one,two
      double precision f,f_a,eps_locala,eps_localb,sumt,sumb,dx,cbig
      double precision x_h(n),f_h,c_h(ni),gf_h(n),gc_h(ni,n),gf(n)
      double precision delxnorm2,delx(n),x(n),gc(ni,n),dp,alpha,fbig
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
        sumt=zero
        sumb=zero
        do i=1,n
          dx=x(i)-x_h(i)
          sumt=sumt+gf(i)*dx
          sumb=sumb+acurvn(i)*dx**2
        enddo
        alpha=min(2.d0*(lloop+1),2.d0*(f-f_h-sumt)/sumb)
        do i=1,n
          if (lloop.gt.0) then
            acurvn(i)=min(fbig,max(acurvn(i),atol2)*fact_inc*factx)
          else 
            acurvn(i)=min(fbig,max(acurvn(i),atol2)*alpha)
          endif
        enddo
        !write(*,*) lloop,acurvn(1)
      endif
!
      if (.not.unconstrained) then
        do j=1,ni
          if (c_a(j).lt.c(j)-eps_localb.and.c(j).gt.-1.d-4) then
            accept=.false.
            sumt=zero
            sumb=zero
            do i=1,n
              dx=x(i)-x_h(i)
              sumt=sumt+gc(j,i)*dx
              sumb=sumb+bcurvn(j,i)*dx**2
            enddo
            alpha=min(2.d0*(lloop+1),2.d0*(c(j)-c_h(j)-sumt)/sumb)
            do i=1,n
              if (lloop.gt.0) then
                bcurvn(j,i)=min(cbig,max(bcurvn(j,i),btol2)
     &                                         *fact_inc*factx)
              else
                bcurvn(j,i)=min(cbig,max(bcurvn(j,i),btol2)*alpha)
              endif  
            enddo
            !write(*,*) lloop,bcurvn(j,1)
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
      end subroutine conserve1

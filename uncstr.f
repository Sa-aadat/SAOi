! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Unconstrained `solver'
! SAOi:


      subroutine uncstr (n, x, f_a, x_h, acurv, acurvn, l, u, f, gf,
     &                   ksubiter, iuser, luser, cuser, ruser, 
     &                   ictrl, lctrl, rctrl, cctrl)
      implicit         none
      include          'ctrl.h'
      integer          n, i, ksubiter
      double precision flambda
      double precision x(n), x_h(n), l(n), u(n), f, f_a
      double precision acurv, acurvn(n), gf(n)
      include          'ctrl_get.inc'
!
      ksubiter=ksubiter+1
!
      do i=1,n
        x(i) = x_h(i)-gf(i)/acurvn(i)
        x(i) = max(l(i),x(i))
        x(i) = min(u(i),x(i))
      enddo
!
      call fun_a (n, x, f_a, x_h, acurv, acurvn, f, gf)
!
      return
      end subroutine uncstr

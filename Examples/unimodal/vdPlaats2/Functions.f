! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: objective and constraint functions
! SAOi:

      subroutine SAOi_funcs (n, m, ni, ne, x, f, c, iuser, luser, cuser,
     &                       ruser, eqn, lin, ictrl, lctrl, rctrl,
     &                       cctrl)
!----------------------------------------------------------------------!
!                                                                      !
!  Compute the objective function f and the inequality constraint      !
!  functions c(j), j=1,ni                                              !
!                                                                      !
!  Please see the users manual for type declarations and comments      !
!                                                                      !
!----------------------------------------------------------------------!
      implicit         none
      include          'ctrl.h'
      logical          eqn(*), lin(*)
      integer          i, j, n, m, ni, ne, nseg
      double precision f, x(n), c(ni)
      double precision y, yy, pload, young, xlen, xleni, aterm, bterm
      double precision ybar, a1, b1, xinert, xmomi, sigi, sbar
      double precision sigmabar
!
      nseg=n/2
      pload=50000.d0
      young=2.d7
      xlen=500.d0
      xleni=xlen/dble(nseg)
      sigmabar=14000.d0
      f=0.d0
      do i=1,nseg
        aterm=(xlen-dble(i)*xleni+2.d0/3.d0*xleni)
        bterm=(xlen+xleni/2.d0-dble(i)*xleni)
        xinert=(x(2*i-1)*x(2*i)**3)/12.d0
        xmomi=pload*(xlen+xleni-dble(i)*xleni)
        sigi=xmomi*x(2*i)/(2.d0*xinert)
        f=f+x(2*i-1)*x(2*i)*xleni
        c(2*i-1)=sigi/sigmabar-1.d0
        c(2*i  )=x(2*i)-20.d0*x(2*i-1)
      enddo
!
      return
      end subroutine SAOi_funcs

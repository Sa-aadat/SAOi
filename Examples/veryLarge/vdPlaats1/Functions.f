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
      double precision y, yy, pload, young, xlen, xleni
      double precision ybar, a1, b1, xinert, xmomi, sigi, sbar
      double precision aterm, bterm, scale1, sigmabar, ydisp, ydispp
!
      scale1   = 1.d3
      ruser(1) = scale1
!
      nseg     = n/2
      ydisp    = 0.d0
      ydispp   = 0.d0
      pload    = 50000.d0
      young    = 2.d7
      xlen     = 500.d0
      xleni    = xlen/dble(nseg)
      sigmabar = 14000.d0
      ybar     = 2.5d0
      f        = 0.d0
      do i=1,nseg
        aterm=(xlen-dble(i)*xleni+2.d0/3.d0*xleni)
        bterm=(xlen+xleni/2.d0-dble(i)*xleni)
        xinert=(x(2*i-1)*x(2*i)**3)/12.d0
!
        ydisp=ydisp+(pload*xleni**2)/(2.d0*young*xinert)*aterm
     &         +ydispp*xleni
        ydispp=ydispp+(pload*xleni)/(young*xinert)*bterm
!
        xmomi=pload*(xlen+xleni-dble(i)*xleni)
        sigi=xmomi*x(2*i)/(2.d0*xinert)
        f=f+x(2*i-1)*x(2*i)*xleni
        c(i)=sigi/sigmabar-1.d0
        c(nseg+i)=x(2*i)-20.d0*x(2*i-1)
      enddo
!
      c(ni)=ydisp/ybar*scale1-1.d0*scale1
!
      return
      end subroutine SAOi_funcs

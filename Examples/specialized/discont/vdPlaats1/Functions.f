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
      double precision f, x(n), c(ni), fact, tmp
      double precision y, yy, pload, young, xlen, slen(n/2), slenc
      double precision ybar, a1, b1, xinert, xmomi, sigi, sbar
!
      nseg  = n/2
      pload = 50000.d0
      young = 2.d7
      xlen  = 500.d0
      slenc = xlen/dble(nseg)
      sbar  = 14000.d0
      ybar  = 2.5d0
      fact  = 1.d-4
      ruser(1) = fact
!
      tmp = 0.d0
      do i=1,nseg-1,2
        slen(i) = slenc+slenc*rctrl(34)*fact ! extend segments by random number*fact
        tmp=tmp + slen(i)
      enddo
!
      do i=2,nseg,2
        slen(i) = slenc-slenc*rctrl(34)*fact ! shorten segments by random number*fact
        tmp=tmp + slen(i)
      enddo
!
      y  = 0.d0
      yy = 0.d0
      f  = 0.d0
      do i=1,nseg
        a1=(xlen-dble(i)*slen(i)+2.d0/3.d0*slen(i))
        b1=(xlen+slen(i)/2.d0-dble(i)*slen(i))
        xinert=(x(2*i-1)*x(2*i)**3)/12.d0
!
        y=y+(pload*slen(i)**2)/(2.d0*young*xinert)*a1+yy*slen(i)
        yy=yy+(pload*slen(i))/(young*xinert)*b1
!
        xmomi=pload*(xlen+slen(i)-dble(i)*slen(i))
        sigi=xmomi*x(2*i)/(2.d0*xinert)
        f=f+x(2*i-1)*x(2*i)*slen(i)
        c(2*i-1)=sigi/sbar-1.d0
        c(2*i  )=x(2*i)-20.d0*x(2*i-1)
      enddo
      c(ni)=y/ybar-1.d0
!
      return
      end subroutine SAOi_funcs

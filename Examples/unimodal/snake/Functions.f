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
      integer          i, j, n, m, ni, ne, l2030, ix, iy, iz
      integer          k, i1, i2, i3, i4
      double precision f, x(n), c(ni)
      double precision pi, zzk, pii2, zzi, alfaj, sumkva, xy1, zxy
!
      pi = 4.0d0*datan(1.0d0)
!
      l2030=10
      k=l2030
      zzk=dble(k)
      pii2=pi/2.d0
      f=0.d0
      sumkva=0.d0
      do i=1,k
        zzi=dble(i)
        alfaj=pii2*zzi/zzk - 2.d0*pii2/3.d0
        ix=i
        iy=k+i
        iz=2*k+i
        f=f+x(ix)*dcos(alfaj)+x(iy)*dsin(alfaj)-0.1d0*x(iz)
        sumkva=sumkva+x(ix)*x(ix)+x(iy)*x(iy)
      enddo
      do i=1,k
        zzi=dble(i)
        ix=i
        iy=k+i
        iz=2*k+i
        xy1 = x(ix)**2 + x(iy)**2 - 1.d0
        zxy = x(iz) - 2.d0*x(ix)*x(iy)
        i1=i
        i2=k+i
        i3=2*k+i
        i4=3*k+i
        c(i1)=10.d0*xy1 + (10.d0*xy1)**7
        c(i2)=10.d0*zxy + (10.d0*zxy)**7
        c(i3)=-c(i1)
        c(i4)=-c(i2)
      enddo
      c(ni)=sumkva
!
      do i=1,ni-1
        c(i)=c(i)-2.d0
      enddo
      c(ni)=c(ni)-dble(k)
!
      return
      end subroutine SAOi_funcs

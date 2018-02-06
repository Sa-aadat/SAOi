! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Specialized instructions, currently for local stress
! SAOi: constrained topology optimization
! SAOi:


      subroutine SAOi_special (n, ni, ne, x, f, c, h, iuser, luser, 
     &                         cuser, ruser, eqn, lin, ictrl, lctrl, 
     &                         rctrl, cctrl, x_lo, x_up, outerloop,
     &                         finished, kkt, violation, nactlo,
     &                         nacthi, nactc, delxnormc, delxnormi,
     &                         time)
!----------------------------------------------------------------------!
!                                                                      !
! Optional specialized instructions may be performed here. This        !
!     routine is called at the end of each OUTER iteration, but        !
!     also at iteration 0. If the algorithm has converged or           !
!     a stop is enforced, finished = .true.                            !
!                                                                      !
! Specialized.out is on unit 15                                        !
!                                                                      !
! Please see the users manual for type declarations and comments       !
!                                                                      !
!----------------------------------------------------------------------!
      implicit         none
      include          'ctrl.h'
      logical          eqn(*), lin(*), finished
      integer          i, j, n, ni, ne, outerloop
      integer          nactlo, nacthi, nactc
      double precision f, x(n), c(ni), h(ne), x_lo(n), x_up(n), time
      double precision kkt, violation, phiBW, delxnormc, delxnormi
      include          'ctrl_get.inc'

! In general, do nothing but return emptyhanded
      if (.not.special) return

! Else, do something clever here, like calculating the black-and-white
!     fraction in topology optimization, and implement stress
!     relaxation for locally constrained topology optimization,
!     etc.
!
      ! the code comes here !
!
      return
      end subroutine SAOi_special
!----------------------------------------------------------------------!
      subroutine filter1(nelx,nely,rmin,x,dc,dcold,n)
      implicit         none
      integer          nelx,nely,n,icol,irow,ii,jj,icount
      integer          imaxk,imaxl,imink,iminl,kk,ll
      double precision x(n),dc(n),dcold(n),rmin,sum,xfac
      double precision xmat(nely,nelx),xii,xjj,xkk,xll
      double precision dcmat(nely,nelx),dcoldmat(nely,nelx)
      
! pack x, dc and dcold in matrix form
      icount=0
      do icol=1,nelx
        do irow=1,nely
          icount=icount+1
          xmat(irow,icol)  = x(icount)
          dcoldmat(irow,icol) = dcold(icount)
          dcmat(irow,icol)    = 0.d0
        enddo
      enddo
 
! now filter sensitivities
      do ii=1,nelx
        do jj=1,nely
          sum=0.d0
          imink=nint(max(dble(ii)-rmin,1.d0))
          imaxk=nint(min(dble(ii)+rmin,dble(nelx)))
          iminl=nint(max(dble(jj)-rmin,1.d0))
          imaxl=nint(min(dble(jj)+rmin,dble(nely)))
          do kk=imink,imaxk
            do ll=iminl,imaxl
              xii = dble(ii)
              xjj = dble(jj)
              xkk = dble(kk)
              xll = dble(ll)
              xfac = rmin - dsqrt((xii-xkk)**2.d0 + (xjj-xll)**2.d0)
              sum = sum + dmax1(0.d0,xfac)
              dcmat(jj,ii) = dcmat(jj,ii) +
     &        (dmax1(0.d0,xfac)*xmat(ll,kk)*dcoldmat(ll,kk))
            enddo
          enddo
          dcmat(jj,ii) = dcmat(jj,ii)/(xmat(jj,ii)*sum)
        enddo
      enddo
 
! repack xmat, dcmat and dcoldmat into vector form
      icount=0
      do icol=1,nelx
        do irow=1,nely
            icount=icount+1
            dc(icount) = dcmat(irow,icol)
        enddo
      enddo
!
      return
      end

! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Specialized instructions, currently for local stress
! SAOi: constrained topology optimization
! SAOi:


      subroutine SAOi_special (n, ni, ne, x, f, c, h, Iuser, Luser, 
     &                         cuser, Ruser, eqn, lin, ictrl, lctrl, 
     &                         rctrl, cctrl, x_lo, x_up, outerloop,
     &                         finished, kkt, violation, nactlo,
     &                         nacthi, nactc)
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
      integer          nactlo, nacthi, nactc, ik,jk
      double precision f, x(n), c(ni), h(ne), x_lo(n), x_up(n)
      double precision kkt, violation, phiBW, eps
      include          'ctrl_get.inc'

! In general, do nothing but return emptyhanded
      if (.not.special) return

! Else, do something clever here, like calculating the black-and-white
!     fraction in topology optimization, and implement stress
!     relaxation for locally constrained topology optimization
!
! For EdSAP:
!
!       rho_min        = Ruser(1)
!       penal          = Ruser(2)
!       volfrac        = Ruser(3)
!       rmin           = Ruser(4)
!       p1             = Ruser(5)
!       p2             = Ruser(6)
!       penalq         = Ruser(7)
!       penal_bruns    = Ruser(8)
!       penal2         = Ruser(9)
!       sall           = Ruser(10)
!       conlim         = Ruser(11)
!       gradlim        = Ruser(12)
!       divval         = Ruser(13)
!       theta          = Ruser(14)
!       eps            = Ruser(15)
!       fscale         = Ruser(16)
!       nelx           = Iuser(1)
!       nely           = Iuser(2)
!       istrgrad       = Iuser(3)
!       iadjoint       = Iuser(4)
!       isdense        = Iuser(5)
!       jacobfilt      = Iuser(6)
!       jblockwrite    = Iuser(7)
!       mcompliance    = Iuser(8) 
!       loaded_knee    = Luser(1)
!
!     if (outerloop.gt.20.and.violation.lt.1.d-3) then
!       Ruser(2) = min(3.d0, Ruser(2)*1.1d0)             ! = penal
!     endif
!      
      open (18,file='DWW.out')
      rewind(18)
      read (18,*) ik,jk
      close(18)
      write(19,*) outerloop,ik,jk
!
      if (outerloop.ge.30.and.violation.lt.1.d-3) then
!       Ruser(2)  = min(3.d0, Ruser(2)*1.1d0)    
!       Ruser(15) = max(1.d-2, Ruser(15)/1.1d0)  
        Ruser(15) = min(1.d0, Ruser(15)*1.1d0)   
!       Ruser(1) = Ruser(15)**2                  
!       do i=1,n
!         x_lo(i)= Ruser(1)
!       end do
      endif
!
      write(15,1000) outerloop,Ruser(1:2),Ruser(4)
!
      if (finished) then
        write(15,*)' '
        phiBW = dble(nactlo+nacthi)/dble(n)
        write(15,*) f, violation, outerloop, kkt, phiBW, nactc
      endif
!
      include          'ctrl_set.inc'
!
      return
!
 1000 format (i9,100f12.4)     
!
      end subroutine SAOi_special
!----------------------------------------------------------------------!
      subroutine filter1(nelx,nely,rmin,x,dc,dcold,n)
      implicit         none
      integer          nelx,nely,n,icol,irow,ii,jj,icount
      integer          imaxk,imaxl,imink,iminl,kk,ll
      double precision x(n),dc(n),dcold(n),rmin,dsum,xfac
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
          dsum=0.d0
          imink=int(max(dble(ii)-rmin,1.d0))
          imaxk=int(min(dble(ii)+rmin,dble(nelx)))
          iminl=int(max(dble(jj)-rmin,1.d0))
          imaxl=int(min(dble(jj)+rmin,dble(nely)))
          do kk=imink,imaxk
            do ll=iminl,imaxl
              xii = dble(ii)
              xjj = dble(jj)
              xkk = dble(kk)
              xll = dble(ll)
              xfac = rmin - dsqrt((xii-xkk)**2.d0 + (xjj-xll)**2.d0)
              dsum = dsum + dmax1(0.d0,xfac)
              dcmat(jj,ii) = dcmat(jj,ii) +
     &        (dmax1(0.d0,xfac)*xmat(ll,kk)*dcoldmat(ll,kk))
            enddo
          enddo
          dcmat(jj,ii) = dcmat(jj,ii)/(xmat(jj,ii)*dsum)
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

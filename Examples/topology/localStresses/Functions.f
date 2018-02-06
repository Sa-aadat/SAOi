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
      integer          i, j, n, ni, ne, m
      double precision f, x(n), c(ni)
      double precision temp0, temp1, temp2, temp3
!
      logical          loaded_knee
      integer          nelx, nely, istrgrad, iadjoint, k, l
      integer          ind, icount,loop, ii, jj, mcompliance
      integer          isdense, jacobfilt, jblockwrite,ik,jk
      double precision penal, volfrac, rmin, p1, p2, penalq, ds
      double precision penal2, rho_min, sall, conlim, gcvol(n)
      double precision eps, eta, penal_bruns, clo, chi      
      double precision gradlim, divval, cstress(n), theta, fscale
!      
      include          'ctrl_get.inc'
!
      write(*,*) ' '
      write(*,*) '-----------------------------------------------'
      write(*,*) 'Entering F '
!
      rho_min        = ruser(1)
      penal          = ruser(2)
      volfrac        = ruser(3)
      rmin           = ruser(4)
      p1             = ruser(5)
      p2             = ruser(6)
      penalq         = ruser(7)
      penal_bruns    = ruser(8)
      penal2         = ruser(9)
      sall           = ruser(10)
      conlim         = ruser(11)
      gradlim        = ruser(12)
      divval         = ruser(13)
      theta          = ruser(14)
      eps            = ruser(15)
      fscale         = ruser(16)
      nelx           = iuser(1)
      nely           = iuser(2)
      istrgrad       = iuser(3)
      iadjoint       = iuser(4)
      isdense        = iuser(5)
      jacobfilt      = iuser(6)
      jblockwrite    = iuser(7)
      mcompliance    = iuser(8) 
      loaded_knee    = luser(1)
!
      if ((jacobfilt.eq.1).and.(gradlim.gt.1.0d10)) then
         jacobfilt = 2
         divval = 1000.0d0
      endif
!         
      open (22,file='DesVars',status='unknown',err=1005)
      rewind(22)
      do i=1,n
        write(22,1245) x(i)
      enddo
      eta=0.5d0
      write(22,1245) penal
      write(22,1245) volfrac
      write(22,1245) rmin
      write(22,1245) penal_bruns   ! volumetric penalization
      write(22,1245) eta
      write(22,1245) penalq        ! gray-scale filter in monotonic OC-methods
      write(22,1245) sall
      write(22,1245) eps
      write(22,1245) conlim
      write(22,1245) gradlim
      write(22,1245) divval
      write(22,1245) theta
      write(22,*) iadjoint
      write(22,*) istrgrad
      write(22,*) isdense
      write(22,*) jacobfilt
      write(22,*) jblockwrite
      close(22)
!
      call system('time ./sctopo8')
      write(*,*) '-----------------------------------------------'
!
      do j = 1,ni
         c(j) = 0.0d0
      enddo
!
      if (istrgrad.ne.0) then
        OPEN (27,FILE='Stresses.out',STATUS='OLD',
     &        FORM='UNFORMATTED',ERR=1005)
        read(27) cstress
        do ii=1,n
           c(ii) = cstress(ii)
        enddo
        close(27)
      endif
!
      f=0.d0
      if (mcompliance.eq.1) then
         open (23,file='TowerFuncs',status='old',err=1005)
         rewind(23)
         read(23,1245) f
         read(23,1245) c(ni)
         close(23)
      else
         do i=1,n
            f=f + fscale*x(i)**penal2
         enddo
      endif
!
!
      call system('rm -f TowerFuncs')
!
!
!       iactblo=0
!       iactbhi=1
!       do i=1,n
!         if (x(i).eq.rho_min) iactblo=iactblo+1
!         if (1.d0.eq.x(i))    iactbhi=iactbhi+1
!       end do
!       phiBW=dble(iactblo+iactbhi)/dble(n)
!       write(15,1111) f,c(1),phiBW,penal,penalq,penal_bruns
!  1111 format (2e15.6,4f10.5)
!
!       clo=  1.d16
!       chi= -1.d16
!       do j=1,ni
!         clo=min(clo,c(j))
!         chi=max(chi,c(j))
!       enddo
!
cc      write(*,*) ' Leaving F - lo and hi constraints are',clo,chi
!      
        !OPEN (18,FILE='DWW.out',STATUS='OLD',ERR=1005)
!
        return
!
 1005 stop ' fatal i/o error in Functions.f. period. '
 1245 format (d22.15)
 1255 format (i7,d24.15)
!
      end subroutine SAOi_funcs

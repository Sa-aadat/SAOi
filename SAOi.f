! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Main driver
! SAOi:


      program SAO_i
!----------------------------------------------------------------------!
!                                                                      !
!  The main driver for the SAOi algorithm                              !
!                                                                      !
!  Please see the users manual for type declarations and comments      !
!                                                                      !
!----------------------------------------------------------------------!
      implicit         none
      include          'ctrl.h'
      logical          eqn(mmax), lin(1)
      integer          i, n, m, ni, ne, loop, sloop, mxloop  
      integer          nnz, nnzh, ierr, iseed
      double precision x(nmax), x_l(nmax), x_u(nmax), shift(nmax)
      double precision lambda(mmax), z(nmax)
      double precision time0, timef, timeg, times, SAOi_seconds
      data             mxloop /1/

! open output files
      call SAOi_open (ictrl,lctrl,rctrl,cctrl,iuser,luser,cuser,ruser)

! execute the SAOi algorithm mxloop times (mxloop > 1 is useful for development and testing only)
      do loop=1,mxloop

! get the initial time
        time0 = SAOi_seconds ()
        
! initialize and pre-process
        call SAOi_defaults (n,m,ni,ne,ictrl,lctrl,rctrl,cctrl,
     &                      nnz,nnzh,eqn,lin,loop,iseed,shift,
     &                      iuser,luser,cuser,ruser)
     
! loop for sequential solution of partitioned problems        
        do sloop = 1,nmax
        
! reset the incremental time for mode extraction        
          if (sloop.gt.1) time0 = SAOi_seconds ()
 
! associated settings for mode extraction
          call SAOi_assocM (ictrl,lctrl,rctrl,cctrl,sloop,
     &                      iuser,luser,cuser,ruser)
          
! apply user settings
          call SAOi_init (n,ni,ne,x,x_l,x_u,ictrl,lctrl,rctrl,cctrl,
     &                    iuser,luser,cuser,ruser,nnz,nnzh,eqn,lin,
     &                    shift)                                     ! add lambda, z to Init later... !

! associated settings
          call SAOi_assoc (ictrl,lctrl,rctrl,cctrl,m,n,lambda,z,
     &                     iuser,luser,cuser,ruser)                  ! move lambda, z to Init later... !

! settings for the test bed (developers only)
          call SAOi_force (n,m,ni,ne,ictrl,lctrl,rctrl,cctrl,
     &                     iuser,luser,cuser,ruser)

! check the settings
          call SAOi_check (n,m,ni,ne,x,x_l,x_u,ictrl,lctrl,rctrl,cctrl,
     &                     nnz,nnzh,eqn,lin,ierr,iseed,shift,sloop,
     &                     lambda,z,iuser,luser,cuser,ruser)         ! move lambda, z to Init later... !
     
! conditially escape 
          if (ierr.ne.0) cycle ! not active - stop now temporarily enforced in SAOi_check

! call the secondary local/global driver
          call SAOi_split (n,m,ni,ne,x,x_l,x_u,timef,timeg,times,
     &                     ictrl,lctrl,rctrl,cctrl,iuser,luser,
     &                     cuser,ruser,nnz,nnzh,eqn,lin,loop,time0,
     &                     ierr,iseed,shift,lambda,sloop,z)

! make global settings available
          include          'ctrl_get.inc'
     
! conditially end loop for sequential solution of partitioned problems 
          if (mayReturn) exit
     
! end loop for sequential solution of partitioned problems
        enddo
     
! compute and output the run times
        call SAOi_time (time0,timef,timeg,times,ictrl,lctrl,rctrl,cctrl,
     &                  iuser,luser,cuser,ruser)

! end the mxloop loop(s) (normally only once)
      enddo

! close output files
      call SAOi_close ()

! finished; suppress IEEE_UNDERFLOW_FLAG IEEE_DENORMAL by not invoking stop
      end program SAO_i
!-----------------------------------------------------------------------

! SAOi: 
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Timer 
! SAOi: 


      double precision function SAOi_seconds()
      implicit         none
!
!     real             dummy(2), etime       ! single precision, i.e. real*4 !
      double precision time1
!
! replace the calls below by something suitable, depending on OS used and FORTRAN version
!
!     SAOi_seconds = dble(etime(dummy))
!
      call cpu_time (time1)
!
      SAOi_seconds = time1
!      
      return
      end function SAOi_seconds
!-----------------------------------------------------------------------

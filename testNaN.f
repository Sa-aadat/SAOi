! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: test for NaN; not all compilers support this. If not, simply comment
! SAOi: out the instructions below
! SAOi:


      subroutine testNaN (a, ierr)
      implicit         none
      integer          ierr
      double precision a
!
      if (isNaN(a)) then
        ! write(*,*) ' Terminal: the (sub)problem returned ',a
        ierr=-5
      endif
!
      return
      end subroutine testNaN
!----------------------------------------------------------------------c

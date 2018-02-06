! SAOi:
! SAOi: Wed Jul 14 14:26:25 SAST 2010, Albert Groenwold, Stellenbosch
! SAOi: Dummy driver for qpa (double) to solve sparse QP-like DQ subproblems,
! SAOi: based on a Galahad example by the original authors
! SAOi:


! THIS VERSION: GALAHAD 2.2 - 23/04/2008 AT 16:30 GMT.

! double precision version !

   subroutine GALAHAD_QPA (n, m, a_ne, f_a, c_l, c_u, x_l, x_u,  &
                           A_row, A_col, A_ptr, x, g, w, y, z,   &
                           A_val, string, ierr, ks, f)

!
   write(*,*) ' '
   write(*,*) ' STOP: This is a dummy driver for GALAHAD_QPA 64. ', &
   'Install Galahad, and run "make galahad" '
   write(*,*) ' '
!
   stop
!
   END SUBROUTINE GALAHAD_QPA

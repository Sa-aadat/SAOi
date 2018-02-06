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
      integer          i, j, n, m, ni, ne
      double precision f, x(n), c(m)
      double precision temp0, temp1, temp2, temp3
      double precision y(10), s, t
!
      y(1)=-6.089d0     
      y(2)=-17.164d0    
      y(3)=-34.054d0    
      y(4)=-5.914d0     
      y(5)=-24.721d0    
      y(6)=-14.986d0    
      y(7)=-24.1d0      
      y(8)=-10.708d0    
      y(9)=-26.662d0    
      y(10)=-22.179d0   
! 
      t=0.d0    
      do i=1,10      
        t=t+dexp(x(i))   
      enddo  
!
      s=0.d0    
      do i=1,10      
        s=s+dexp(x(i))*(y(i)+x(i)-dlog(t))
      enddo
      f=s   
!
      c(1)=dexp(x(1))+2.d0*dexp(x(2))+2.d0*dexp(x(3))  
     &    +dexp(x(6))+dexp(x(10))-2.d0  
      c(2)=dexp(x(4))+2.d0*dexp(x(5))+dexp(x(6))+dexp(x(7))-1.d0    
      c(3)=dexp(x(3))+dexp(x(7))+dexp(x(8))      
     &    +2.d0*dexp(x(9))+dexp(x(10))-1.d0  
!
      return
      end subroutine SAOi_funcs

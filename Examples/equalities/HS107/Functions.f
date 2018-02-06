! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: objective and constraint functions
! SAOi:

      subroutine SAOi_funcs (n, m, ni, ne, x, f, c, iuser, luser, cuser,
     &                       ruser, eqn, lin, ictrl, lctrl, rctrl,
     &                       cctrl)
!----------------------------------------------------------------------!
!                                                                      !
!  Compute the objective function f and the (in)equality constraint    !
!  functions c(j), j=1,m                                              !
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
      double precision y1, y2, y3, y4, y5, y6, b, d, v1
!
      f    = 3.d+3*x(1)+1.d+3*x(1)**3+2.d+3*x(2)+666.66666667d0*x(2)**3
!
      y1=sin(x(8))     
      y2=cos(x(8))     
      y3=sin(x(9))     
      y4=cos(x(9))     
      y5=sin(x(8)-x(9))
      y6=cos(x(8)-x(9))
      v1=48.4d0/50.176d0
      b=v1*sin(0.25d0) 
      d=v1*cos(0.25d0)
! 
      c(1)=0.4d0-x(1)+2.d0*b*x(5)**2+x(5)*x(6)
     &    *(-d*y1-b*y2)+x(5)*x(7)*(-d*y3-b*y4)    
      c(2)=0.4d0-x(2)+2.d0*b*x(6)**2+x(5)*x(6)
     &    *(d*y1-b*y2)+x(6)*x(7)*(d*y5-b*y6)     
      c(3)=0.8d0+2.d0*b*x(7)**2+x(5)*x(7)
     &    *(d*y3-b*y4)+x(6)*x(7)*(-d*y5-b*y6)   
      c(4)=0.2d0-x(3)+2.d0*d*x(5)**2-x(5)
     &    *x(6)*(-b*y1+d*y2)-x(5)*x(7)*(-b*y3+d*y4)     
      c(5)=0.2d0-x(4)+2.d0*d*x(6)**2
     &    -x(5)*x(6)*(b*y1+d*y2)-x(6)*x(7)*(b*y5+d*y6)       
      c(6)=-0.337d0+2.d0*d*x(7)**2-x(5)*x(7)     
     &    *(b*y3+d*y4)-x(6)*x(7)*(-b*y5+d*y6)     
!
      return
      end subroutine SAOi_funcs

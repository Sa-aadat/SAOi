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
      double precision y(10), s, t, dlogt
!
      c(1)=x(1)+2.d0*x(2)+2.d0*x(3)+x(6)+x(10)-2.0   
      c(2)=x(4)+2.d0*x(5)+x(6)+x(7)-1.0 
      c(3)=x(3)+x(7)+x(8)+2.d0*x(9)+x(10)-1.0  
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
      do 30 i=1,10      
30    t=t+x(i)    
!  
      if (t.lt.1.d-5) goto 34    
      dlogt=dlog(t)   
!  
      s=0.d0    
      do 31 i=1,10      
      if (x(i).lt.0.d0) goto 34  
   31 s=s+x(i)*(y(i)+dlog(x(i))-dlogt)  
      f=s
!
      return    
!
   34 s=0.d0    
!
      do 35 i=1,10      
   35 if (x(i).lt.0.d0) s=s+(x(i)-5.d0)**2       
      f=(s+1.d+3-47.8d0)
!
      return
      end subroutine SAOi_funcs

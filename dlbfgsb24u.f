! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Driver for l-bfgs-b 
! SAOi: This file was originally called driver1/2/3.f; it originated
! SAOi: with Ciyou Zhu
! SAOi:



c     --------------------------------------------------------------
c                SIMPLE DRIVER FOR L-BFGS-B (version 2.4)
c     --------------------------------------------------------------
c
c        L-BFGS-B is a code for solving large nonlinear optimization
c        problems with simple bounds on the variables.
c
c        The code can also be used for unconstrained problems and is
c        as efficient for these problems as the earlier limited memory
c        code L-BFGS.
c
c     References:
c
c        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c        memory algorithm for bound constrained optimization'',
c        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
c        Subroutines for Large Scale Bound Constrained Optimization''
c        Tech. Report, NAM-11,EECS Department, Northwestern University,
c        1994.
c
c        (Postscript files of these papers are available via anonymous
c        ftp to ece.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
c
c        NEOS,November 1994. (Latest revision April 1997.)
c        Optimization Technology Center.
c        Argonne National Laboratory and Northwestern University.
c
c        Written by
c                                Ciyou Zhu
c        in collaboration with R.H. Byrd,P. Lu-Chen and J. Nocedal.
c
c
c     NOTE: The user should adapt the subroutine 'timer' if 'etime' is
c           not available on the system.  An example for system
c           AIX Version 3.2 is available at the end of this driver.
c
c     **************
       subroutine drive_lbfgsb24u (n,x,ni,ne,x_h,
     &                            acurv,bcurv,ccurv,acurvn,
     &                            bcurvn,ccurvn,x_l,x_u,
     &                            f,c,h,gf,gc,gh,ksubiter,
     &                            xkkt,flam,cstage0,ostring,
     &                            ns,ictrl,lctrl,rctrl,cctrl,message,
     &                            yhi,outerloop,eqn,lin,acol,nnz,
     &                            do_diff,nfe,nge,
     &                            iuser,luser,cuser,ruser)

!  This driver calls the l-bfgs-b code for unconstrained (box constrained) PRIMAL problems

      implicit         none
      include          'ctrl.h'
      integer          umax
      parameter        (umax = 8)
      character*60     task,csave
      character*90     ostring
      logical          eqn(*),lin(*)
      logical          lsave(4),cstage0,warn,do_diff
      integer          ni,ne,n,ksubiter,m,iprints,acol(*),nnz
      integer          i,j,nbd(n),iwa(3*n),nfe,nge
      integer          isave(44),isparse,outerloop
      integer          ns,ncol(ns),nrow(ns),icnt,message,ifree
      double precision flambda,factr,pgtol,yhi
      double precision time1,time2
      double precision dsave(29),x(n),flam
      double precision wa(2*umax*n + 4*n + 11*umax*umax + 8*umax)
      double precision x_h(n),x_l(n),x_u(n) 
      double precision acurv,bcurv(ni),ccurv(ne)
      double precision acurvn(n),bcurvn(ni,n)
      double precision ccurvn(ne,n),xlamsml,xlambig
      double precision f,c(ni),h(ne),gf(n),gc(ni,n)
      double precision xkkt,gf_a(n),gc_a(ni,n)
      double precision gh(ne,n)
      double precision timef1,timef2,timef3
      include          'ctrl_get.inc'

! no output whatsoever
      iprints = -10

! specify the tolerances in the stopping criteria
      factr  =  1.0d1
      pgtol  =  1.0d-10

! specify the dimension n and the number m of limited memory corrections stored
      m      = 8
!
      if (m.gt.umax ) stop ' m.gt.umax  in dlbfgsb24u.f'

! set bounds on the variables
      do 10 i = 1,n
        nbd(i) = 2
   10 continue

! start the iteration by initializing task and timer
      call timer(time1)
      task = 'START'
      ostring=' subproblem terminated with standard settings'
      message=0

! entry point of the subproblem optimization loop
  111 continue

! call the l-bfgs-b code
      call setulb24(n,m,x,x_l,x_u,nbd,f,gf,factr,pgtol,wa,iwa,task,
     &              iprints,csave,lsave,isave,dsave)
!
      if (task(1:2) .eq. 'FG') then

! check the time limit
        call timer(time2)
!
        if (time2 - time1 .gt. tlimit) then
          task='STOP: CPU time exceeds limit.'
          ostring=' subproblem terminated on time limit'
          message=1
!
        else
!
          call SAOi_funcs (n,ni,ni-ne,ne,x,f,c,iuser,luser,cuser,ruser,
     &                     eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
          call formgrad (n,x,ni,ne,f,c,gf,gc,ictrl,lctrl,rctrl,cctrl,
     &                   nfe,nge,do_diff,iuser,luser,cuser,ruser,nnz,
     &                   Acol,ncol,eqn,lin)
!
          goto 111
        endif
!
      elseif (task(1:5) .eq. 'NEW_X') then
! the minimization routine has returned with a new iterate,and we have opted to (conditionally) continue the iterations

! terminate if the total number of f and g evaluations exceeds ksubmax
        if (isave(34) .ge. ksubmax) then
          task='STOP: allowed no. of subproblem evaluations exceeded'
          ostring=' subproblem terminated on number of samples of the du
     &al'
          message=2
        endif

! terminate if  |proj g|/(1 + |f|) < 1.0d-10.
        if (dsave(13) .le. 1.d-10*(1.0d0 + dabs(f))) then
          task='STOP: projected gradient is sufficiently small'
          ostring=' subproblem terminated with |proj g|/(1 + |f|) < 1.0d
     &-10'
          message=0
        endif
        goto 111

! terminate execution when task is neither FG nor NEW_X
      else
        ! do nothing

! end of the loop
      endif
!
      ksubiter = isave(34)
!
      return
      end subroutine drive_lbfgsb24u

c======================= The end of driver1 ============================
c     --------------------------------------------------------------
c             DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------
c
c     n is an INTEGER variable that must be set by the user to the
c       number of variables.  It is not altered by the routine.
c
c     m is an INTEGER variable that must be set by the user to the
c       number of corrections used in the limited memory matrix.
c       It is not altered by the routine.  Values of m < 3  are
c       not recommended,and large values of m can result in excessive
c       computing time. The range  3 <= m <= 20 is recommended.
c
c     x is a DOUBLE PRECISION array of length n.  On initial entry
c       it must be set by the user to the values of the initial
c       estimate of the solution vector.  Upon successful exit,it
c       contains the values of the variables at the best point
c       found (usually an approximate solution).
c
c     l is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the lower bounds on the variables. If
c       the i-th variable has no lower bound,l(i) need not be defined.
c
c     u is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the upper bounds on the variables. If
c       the i-th variable has no upper bound,u(i) need not be defined.
c
c     nbd is an INTEGER array of dimension n that must be set by the
c       user to the type of bounds imposed on the variables:
c       nbd(i)=0 if x(i) is unbounded,
c              1 if x(i) has only a lower bound,
c              2 if x(i) has both lower and upper bounds,
c              3 if x(i) has only an upper bound.
c
c     flambda is a DOUBLE PRECISION variable.  If the routine setulb24 returns
c       with task(1:2)= 'FG',then flambda must be set by the user to
c       contain the value of the function at the point x.
c
c     g is a DOUBLE PRECISION array of length n.  If the routine setulb24
c       returns with taskb(1:2)= 'FG',then g must be set by the user to
c       contain the components of the gradient at the point x.
c
c     factr is a DOUBLE PRECISION variable that must be set by the user.
c       It is a tolerance in the termination test for the algorithm.
c       The iteration will stop when
c
c        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
c
c       where epsmch is the machine precision which is automatically
c       generated by the code. Typical values for factr on a computer
c       with 15 digits of accuracy in double precision are:
c       factr=1.d+12 for low accuracy;
c             1.d+7  for moderate accuracy;
c             1.d+1  for extremely high accuracy.
c       The user can suppress this termination test by setting factr=0.
c
c     pgtol is a double precision variable.
c       On entry pgtol >= 0 is specified by the user.  The iteration
c         will stop when
c
c                 max{|proj g_i | i = 1,...,n} <= pgtol
c
c         where pg_i is the ith component of the projected gradient.
c       The user can suppress this termination test by setting pgtol=0.
c
c     wa is a DOUBLE PRECISION  array of length
c       (2umax + 4)nmaxb + 11umax^2 + 8umax used as workspace.
c       This array must not be altered by the user.
c
c     iwa is an INTEGER  array of length 3nmaxb used as
c       workspace. This array must not be altered by the user.
c
c     task is a CHARACTER string of length 60.
c       On first entry,it must be set to 'START'.
c       On a return with task(1:2)='FG',the user must evaluate the
c         function f and gradient g at the returned value of x.
c       On a return with task(1:5)='NEW_X',an iteration of the
c         algorithm has concluded,and f and g contain f(x) and g(x)
c         respectively.  The user can decide whether to continue or stop
c         the iteration.
c       When
c         task(1:4)='CONV',the termination test in L-BFGS-B has been
c           satisfied;
c         task(1:4)='ABNO',the routine has terminated abnormally
c           without being able to satisfy the termination conditions,
c           x contains the best approximation found,
c           f and g contain f(x) and g(x) respectively;
c         task(1:5)='ERROR',the routine has detected an error in the
c           input parameters;
c       On exit with task = 'CONV','ABNO' or 'ERROR',the variable task
c         contains additional information that the user can print.
c       This array should not be altered unless the user wants to
c          stop the run for some reason.  See driver2 or driver3
c          for a detailed explanation on how to stop the run
c          by assigning task(1:4)='STOP' in the driver.
c
c     iprints is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprints<0    no output is generated;
c        iprints=0    print only one line at the last iteration;
c        0<iprints<99 print also f and |proj g| every iprints iterations;
c        iprints=99   print details of every iteration except n-vectors;
c        iprints=100  print also the changes of active set and final x;
c        iprints>100  print details of every iteration including x and g;
c       When iprints > 0,the file iterate.dat will be created to
c                        summarize the iteration.
c
c     csave  is a CHARACTER working array of length 60.
c
c     lsave is a LOGICAL working array of dimension 4.
c       On exit with task = 'NEW_X',the following information is
c         available:
c       lsave(1) = .true.  the initial x did not satisfy the bounds;
c       lsave(2) = .true.  the problem contains bounds;
c       lsave(3) = .true.  each variable has upper and lower bounds.
c
c     isave is an INTEGER working array of dimension 44.
c       On exit with task = 'NEW_X',it contains information that
c       the user may want to access:
c         isave(30) = the current iteration number;
c         isave(34) = the total number of function and gradient
c                         evaluations;
c         isave(36) = the number of function value or gradient
c                                  evaluations in the current iteration;
c         isave(38) = the number of free variables in the current
c                         iteration;
c         isave(39) = the number of active constraints at the current
c                         iteration;
c
c       See the subroutine setulb24.f for a description of other
c       information contained in isave.lloopmax
c
c     dsave is a DOUBLE PRECISION working array of dimension 29.
c       On exit with task = 'NEW_X',it contains information that
c       the user may want to access:
c         dsave(2) = the value of f at the previous iteration;
c         dsave(5) = the machine precision epsmch generated by the code;
c         dsave(13) = the infinity norm of the projected gradient;
c
c       See the subroutine setulb24.f for a description of other
c       information contained in dsave.
c
c     --------------------------------------------------------------
c           END OF THE DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------
c
c     << An example of subroutine 'timer' for AIX Version 3.2 >>
c
c     subroutine timer(ttime)
c     double precision ttime
c     integer itemp,integer mclock
c
c     itemp = mclock()
c     ttime = dble(itemp)*1.0d-2
c     return
c     end
c=======================================================================
       subroutine drive_lbfgsb24u2 (n,x,f0,c,eqn,lin,ictrl,lctrl,
     &                              rctrl,cctrl,iuser,luser,cuser,ruser)

!  This driver calls the l-bfgs-b code for unconstrained (box constrained) DUAL problems
!    to estimate the initial dual variables at the starting point of an equality constrained
!    problem. Hence, n here in reality is m.

      implicit         none
      include          'ctrl.h'
      integer          umax
      parameter        (umax = 8)
      character*60     task,csave
      character*90     ostring
      logical          eqn(*),lin(*)
      logical          lsave(4) 
      integer          n,m,i,j,nbd(n),iwa(3*n),isave(44)
      integer          iprints,message 
      double precision factr,pgtol 
      double precision time1,time2
      double precision dsave(29),x(n) 
      double precision wa(2*umax*n + 4*n + 11*umax*umax + 8*umax)
      double precision x_l(n),x_u(n)
      double precision f,f0,c(n),gf(n)
      include          'ctrl_get.inc'

! no output whatsoever
      iprints = -10

! specify the tolerances in the stopping criteria
      factr  =  1.0d1
      pgtol  =  1.0d-10

! specify the dimension n and the number m of limited memory corrections stored
      m      = 8
!
      if (m.gt.umax ) stop ' m.gt.umax  in dlbfgsb24u.f'

! set bounds on the variables
      do 10 i = 1,n
        nbd(i) = 2
        x_l(i) = -biglam
        x_u(i) =  biglam
   10 continue
   
! get the initial function value (the Lagrangian)
      f = f0 
      do i = 1,n
        f = f + x(i)*c(i)
      enddo  
      f = f**2
      do i = 1,n
        gf(i) = 2.d0*dsqrt(f)*c(i)
      enddo  

! start the iteration by initializing task and timer
      call timer(time1)
      task = 'START'
      ostring=' subproblem terminated with standard settings'
      message=0

! entry point of the subproblem optimization loop
  111 continue

! call the l-bfgs-b code
      call setulb24(n,m,x,x_l,x_u,nbd,f,gf,factr,pgtol,wa,iwa,task,
     &            iprints,csave,lsave,isave,dsave)
!
      if (task(1:2) .eq. 'FG') then

! check the time limit
        call timer(time2)
!
        if (time2 - time1 .gt. tlimit) then
          task='STOP: CPU time exceeds limit.'
          ostring=' subproblem terminated on time limit'
          message=1
!
        else
!
          f = f0  ! form objective and derivatives
          do i = 1,n
            f = f + x(i)*c(i)
          enddo  
          f = f**2
          do i = 1,n
            gf(i) = 2.d0*dsqrt(f)*c(i)
          enddo  
!
          goto 111
!          
        endif
!
      elseif (task(1:5) .eq. 'NEW_X') then
! the minimization routine has returned with a new iterate,and we have opted to (conditionally) continue the iterations

! terminate if the total number of f and g evaluations exceeds ksubmax
        if (isave(34) .ge. ksubmax) then
          task='STOP: allowed no. of subproblem evaluations exceeded'
          ostring=' subproblem terminated on number of samples of the du
     &al'
          message=2
        endif

! terminate if  |proj g|/(1 + |f|) < 1.0d-10.
        if (dsave(13) .le. 1.d-10*(1.0d0 + dabs(f))) then
          task='STOP: projected gradient is sufficiently small'
          ostring=' subproblem terminated with |proj g|/(1 + |f|) < 1.0d
     &-10'
          message=0
        endif
        goto 111

! terminate execution when task is neither FG nor NEW_X
      else
        ! do nothing

! end of the loop
      endif
!
!     ksubiter = isave(34)
!
      return
      end subroutine drive_lbfgsb24u2
c=======================================================================

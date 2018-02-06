!
!-----------------------------------------------------------------------------!
!
!  IMPORTANT: This is a preimplementation of the Cplex solver in the SAOi
!             algorithm of Albert A. Groenwold
!
!             Only the dense and pseudo-sparse implementation is currently
!             implemented.
!
! This file contains subroutines which are used to call the cplex QP solver
! from SAOi_outer.f  
!
!-----------------------------------------------------------------------------!
!
       subroutine drive_cplex (nprimal,xprimal,ni,nq,f_a,c_a,
     &                         h_a,iact,nact,xlam,x_h,
     &                         acurv,bcurv,ccurv,acurvn,
     &                         bcurvn,ccurvn,x_l,x_u,
     &                         f,c,h,gf,gc,gh,ksubiter,
     &                         xkkt,flam,cstage0,ostring,xlam_h,
     &                         nsx,ictrl,lctrl,rctrl,cctrl,message,yhi,
     &                         outerloop,lloop,eqn,
     &                         iuser,luser,cuser,ruser)
!
!-----------------------------------------------------------------------------!
!
!  Here a description of this routine
!
! Input
!    nprimal  = Number of optimizations variables, used for allocation of
!               variables.
!    xprimal  =
!    ni       = Number of inequality constraints, used for allcoation of
!               variables.
!    nq       = (Not used) Number of equality constraints, not implemented yet.
!    f_a      = approximation of the objective function value based on the 
!               updated point xprimal.
!    c_a      = (Not used)
!    h_a      = (Not used)
!    iact     = (Not used)
!    nact     = (Not used)
!    xlam     =
!    x_h      = (Not used) history data of x befor xprimal
!    acurv    = Curvature data of f in the point xprimal
!    bcurv    = Curvature data of c in the point xprimal
!    ccurv    = (Not used)
!    acurvn   = (Not used)
!    bcurvn   = (Not used)
!    ccurvn   = (Not used)
!    x_l      = Lowerbounds on xprimal
!    x_u      = Uperbounds on xprimal
!    f        = Objective function value in xprimal
!    c        = Inequality constraint function values in xprimal
!    h        = (Not used) equality constraint function values in xprimal
!    gf       = Gradient of f in xprimal
!    gc       = Gradients of c in xprimal
!    gh       = (Not used) Gradients of h in xprimal
!    ksubiter = Number of subiterations at start this is 0
!    xkkt     = (Not used)
!    flam     = (Not used)
!    cstage0  = (Not used)
!    ostring  = (Not used)
!    xlam_h   = (Not used)
!    nsx      = (Not used)
!    ictrl    = Integer option values (see ctrl.h)
!    lctrl    = Logical option values (see ctrl.h)
!    rctrl    = Real option values (see ctrl.h)
!    message  = Output message which gives the status of the solving of the 
!               sub problem (0 if everything goes well)
!    yhi      = (Not used)
!    outerloop= Current number of outerloop (starting at 0)
!    lloop    = Current number of innerloop (starting at 0)
!    ictrl    = Integer global data array which contains the pointer to the
!               Matlab engine.
!    lctrl    = Logical global data array, not used here.
!    rctrl    = Real global data array, not used here.
!    cctrl    = Character global data array which contains the name of the
!               Matlab function file which is called in this subroutine.
!
!  Output
!    xprimal  = New calculated value for xprimal based on step in x.
!    f_a      = Approximated value of objective function in the updated point
!               xprimal
!    xlam     = Updated values for lagranges multipliers based on the found optimum
!
!  Local variables
!    i        = Index variable used for loops.
!    j        = Index variable used for loops.
!    status   = status variables used for error checking of calls to CPLEX
!    qpdata   = Q array of the subproblem
!    bl       = Lower bounds of the subproblem
!    bu       = Upper bounds of the subproblem
!    x        = step in xprimal which is determined by the CPLEX solver
!
!-----------------------------------------------------------------------------!
!
!.....Variable declarations
      implicit         none
      include          'ctrl.h'
      logical          cstage0, eqn(*)
      integer          n, ni,nq,nprimal,ksubiter,m,outerloop
      integer          nsx,nact,iact(nact),message
      integer          loops, lloop, i, j
      integer          status
      double precision xprimal(nprimal),flam, bu(nprimal), bl(nprimal)
      double precision x_h(nprimal),x_l(nprimal),x_u(nprimal)
      double precision qpdata(nprimal)
      double precision xlam(ni),gh(nq,nprimal), yhi
      double precision acurv,bcurv(ni),ccurv(nq),xlam_h(ni)
      double precision acurvn(nprimal),bcurvn(ni,nprimal)
      double precision ccurvn(nq,nprimal)
      double precision f,c(ni),h(nq),gf(nprimal),gc(ni,nprimal)
      double precision f_a,c_a(ni),h_a(nq)
      double precision xkkt, xlamj, bji, x(nprimal)
      character*90     ostring
!
!.....Get global option variables
      include          'ctrl_get.inc'
!
!.....determine if nprimal does not exceed nmax.
      if (nprimal.gt.nmax) stop ' nprimal.gt.nmax in drive_cplex.f'
!
      write(*,*) ' '
      stop ' This a dummy FORTRAN driver for CPLEX'
      write(*,*) ' '
! 
      return
 1000 format (100('='),/,'outerloop, innerloop: ',i5,',',i5)
 1001 format (/,'n: ',i7,/,'ni:',i7)
 1002 format (/,a30,/,10es18.4)
 1003 format (/,10es18.4)
 
      end
      
      subroutine cplex_stop (status)
      integer           status
!
      status = -1
!      
      write(*,*) ' '
      stop ' This a dummy FORTRAN driver for CPLEX'
      write(*,*) ' '
!
      end

      
       subroutine drive_cplexs_old (nprimal,xprimal,ni,nq,f_a,c_a,
     &                         iact,nact,xlam,x_h,
     &                         acurv,bcurv,amult,bmult,x_l,x_u,
     &                         f,c,gf,gc,ksubiter,
     &                         xkkt,flam,cstage0,ostring,xlam_h,
     &                         nsx,ictrl,lctrl,rctrl,cctrl,message,yhi,
     &                         outerloop,lloop,nnz,Acol,Aptr,
     &                         storind, delta_t,eqn,
     &                         iuser,luser,cuser,ruser)
!
!-----------------------------------------------------------------------------!
!
!  Here a description of this routine
!
! Input
!    nprimal  = Number of optimizations variables, used for allocation of
!               variables.
!    xprimal  =
!    ni       = Number of inequality constraints, used for allcoation of
!               variables.
!    nq       = (Not used) Number of equality constraints, not implemented yet.
!    f_a      = approximation of the objective function value based on the 
!               updated point xprimal.
!    c_a      = (Not used)
!    h_a      = (Not used)
!    iact     = (Not used)
!    nact     = (Not used)
!    xlam     =
!    x_h      = (Not used) history data of x befor xprimal
!    acurv    = Curvature data of f in the point xprimal
!    bcurv    = Curvature data of c in the point xprimal
!    ccurv    = (Not used)
!    acurvn   = (Not used)
!    bcurvn   = (Not used)
!    ccurvn   = (Not used)
!    x_l      = Lowerbounds on xprimal
!    x_u      = Uperbounds on xprimal
!    f        = Objective function value in xprimal
!    c        = Inequality constraint function values in xprimal
!    h        = (Not used) equality constraint function values in xprimal
!    gf       = Gradient of f in xprimal
!    gc       = Gradients of c in xprimal
!    gh       = (Not used) Gradients of h in xprimal
!    ksubiter = Number of subiterations at start this is 0
!    xkkt     = (Not used)
!    flam     = (Not used)
!    cstage0  = (Not used)
!    ostring  = (Not used)
!    xlam_h   = (Not used)
!    nsx      = (Not used)
!    ictrl    = Integer option values (see ctrl.h)
!    lctrl    = Logical option values (see ctrl.h)
!    rctrl    = Real option values (see ctrl.h)
!    message  = Output message which gives the status of the solving of the 
!               sub problem (0 if everything goes well)
!    yhi      = (Not used)
!    outerloop= Current number of outerloop (starting at 0)
!    lloop    = Current number of innerloop (starting at 0)
!    ictrl    = Integer global data array which contains the pointer to the
!               Matlab engine.
!    lctrl    = Logical global data array, not used here.
!    rctrl    = Real global data array, not used here.
!    cctrl    = Character global data array which contains the name of the
!               Matlab function file which is called in this subroutine.
!
!  Output
!    xprimal  = New calculated value for xprimal based on step in x.
!    f_a      = Approximated value of objective function in the updated point
!               xprimal
!    xlam     = Updated values for lagranges multipliers based on the found optimum
!
!  Local variables
!    i        = Index variable used for loops.
!    j        = Index variable used for loops.
!    status   = status variables used for error checking of calls to CPLEX
!    qpdata   = Q array of the subproblem
!    bl       = Lower bounds of the subproblem
!    bu       = Upper bounds of the subproblem
!    x        = step in xprimal which is determined by the CPLEX solver
!
!-----------------------------------------------------------------------------!
!
!.....Variable declarations
      implicit         none
      include          'ctrl.h'
      logical          cstage0, eqn(*)
      integer          ni,nq,nprimal,ksubiter,m,outerloop
      integer          nnz, Aptr(ni+1), Acol(nnz), storind(nnz)
      integer          nsx,nact,iact(nact),message
      integer          loops, lloop, i, j, ij, i1
      integer          status
      double precision xprimal(nprimal),flam, bu(nprimal), bl(nprimal)
      double precision x_h(nprimal),x_l(nprimal),x_u(nprimal)
      double precision qpdata(nprimal), delta_t, qpdatas
      double precision t_c1, t_c2, SAOi_seconds
      double precision xlam(ni),yhi
      double precision acurv,bcurv(ni),xlam_h(ni)
      double precision f,c(ni),gf(nprimal),gc(nnz)
      double precision f_a,c_a(ni)
      double precision xkkt, xlamj, bji, x(nprimal), gcij
      character*90     ostring
      double precision amult, bmult(ni)

!
!.....Get global option variables
      include          'ctrl_get.inc'
!
!.....determine if nprimal does not exceed nmax.
      if (nprimal.gt.nmax) stop ' nprimal.gt.nmax in drive_cplex.f'
!
      write(*,*) ' '
      stop ' This a dummy FORTRAN driver for CPLEX'
      write(*,*) ' '
! 
      return
 1000 format (100('='),/,'outerloop, innerloop: ',i5,',',i5)
 1001 format (/,'n: ',i7,/,'ni:',i7)
 1002 format (/,a30,/,10es18.4)
 1003 format (/,10es18.4)
  
      end
      
      subroutine drive_cplexSOCP(nprimal,xprimal,ni,nq,f_a,c_a,
     &                         h_a,iact,nact,xlam,x_h, f_h, c_h,
     &                         acurv,bcurv,ccurv,acurvn,
     &                         bcurvn,ccurvn,amult,bmult,x_l,x_u,
     &                         f,c,h,gf,gc,gh,ksubiter,
     &                         xkkt,flam,cstage0,ostring,xlam_h,
     &                         nsx,ictrl,lctrl,rctrl,cctrl,message,yhi,
     &                         outerloop,lloop,nnz,Acol,Aptr,
     &                         storind,delta_t,eqn,
     &                         iuser,luser,cuser,ruser)
!.....Variable declarations
      implicit         none
      include          'ctrl.h'
      logical          cstage0, eqn(ni)
      integer          ni,nq,nprimal,ksubiter,m,outerloop
      integer          nnz, Aptr(ni+1), Acol(nnz), storind(nnz)
      integer          nsx,nact,iact(nact),message
      integer          loops, lloop, i, j, ij, i1, approxc
      integer*8        status
      double precision xprimal(nprimal),flam, bu(nprimal), bl(nprimal)
      double precision x_h(nprimal),x_l(nprimal),x_u(nprimal)
      double precision qodata(nprimal), qcdata(nnz), delta_t, qpdatas
      double precision t_c1, t_c2, SAOi_seconds, qcdata_s(ni)
      double precision xlam(ni),gh(nq,nprimal), yhi, acurv, bcurv(ni)
      double precision ccurv(nq),xlam_h(ni), f_h, c_h(ni)
      double precision acurvn(nprimal),bcurvn(ni,nprimal)
      double precision ccurvn(nq,nprimal)
      double precision f,c(ni),h(nq),gf(nprimal),gc(nnz)
      double precision f_a,c_a(ni),h_a(nq) 
      double precision xkkt, xlamj, bji, x(nprimal), gcij
      character*90     ostring
      double precision amult, bmult(ni), delx(nprimal), gcji

!.....Get global option variables
      include          'ctrl_get.inc'
!
      if (nprimal.gt.nmax) stop ' nprimal.gt.nmax in drive_cplex.f'
!

! 
!        write(*,*) (t_c2 - t_c1)
      return
 1000 format (100('='),/,'outerloop, innerloop: ',i5,',',i5)
 1001 format (/,'n: ',i7,/,'ni:',i7)
 1002 format (/,a30,/,10es18.4)
 1003 format (/,10es18.4)
 
      end

        subroutine drive_SOCPd (nprimal,xprimal,ni,nq,f_a,c_a,
     &                         h_a,iact,nact,xlam,x_h,
     &                         acurv,bcurv,ccurv,acurvn,
     &                         bcurvn,ccurvn,x_l,x_u,
     &                         f,c,h,gf,gc,gh,ksubiter,
     &                         xkkt,flam,cstage0,ostring,xlam_h,
     &                         nsx,ictrl,lctrl,rctrl,cctrl,message,yhi,
     &                         outerloop,lloop,eqn,
     &                         iuser,luser,cuser,ruser)
!.....Variable declarations
      implicit         none
      include          'ctrl.h'
      logical          cstage0, eqn(ni)
      integer          n, ni,nq,nprimal,ksubiter,m,outerloop
      integer          nsx,nact,iact(nact),message
      integer          loops, lloop, i, j, k 
      integer*8        status
      double precision xprimal(nprimal),flam, bu(nprimal), bl(nprimal)
      double precision x_h(nprimal),x_l(nprimal),x_u(nprimal)
      double precision qodata(nprimal),qcdata(ni*nprimal)
      double precision qpdata(nprimal)
      double precision xlam(ni),gh(nq,nprimal), yhi
      double precision acurv,bcurv(ni),ccurv(nq),xlam_h(ni)
      double precision acurvn(nprimal),bcurvn(ni,nprimal)
      double precision ccurvn(nq,nprimal)
      double precision f,c(ni),h(nq),gf(nprimal),gc(ni,nprimal)
      double precision f_a,c_a(ni),h_a(nq) 
      double precision xkkt, xlamj, bij, x(nprimal)
      character*90     ostring
!      integer           Acol(*), Aptr(*), nnz
!
!.....Get global option variables
      include          'ctrl_get.inc'
!
!.....determine if nprimal does not exceed nmax.
      if (nprimal.gt.nmax) stop ' nprimal.gt.nmax in drive_cplex.f'
!
      return
 1000 format (100('='),/,'outerloop, innerloop: ',i5,',',i5)
 1001 format (/,'n: ',i7,/,'ni:',i7)
 1002 format (/,a30,/,10es18.4)
 1003 format (/,10es18.4)
 
      end

  

   

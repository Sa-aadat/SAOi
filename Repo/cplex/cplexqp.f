!                            #####CPLEX.f#####
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
!.....write input data to a file by first opening the file and then 
      if (debug) then
        if (outerloop.eq.1.and.lloop.eq.0) then
          open (20,file='cplex_loops.out')
        endif
        write(20,1000) outerloop, lloop
        write(20,1001) nprimal, ni
        write(20,1002) 'xstart: ', xprimal
        write(20,1002) 'x_low:', x_l
        write(20,1002) 'x_heigh:', x_u
        write(20,1002) 'fk: ', f  
        write(20,1002) 'ck: ', c
        write(20,1002) 'dFk:', gf
        write(20,1002) 'dCk:', gc(1,:)
        if (ni.ge.2) then
          do j=2,ni
            write(20,1003) gc(j,:)
          enddo
        endif
        write(20,1002) 'Qf:', acurvn
        write(20,1002) 'Qc:', bcurvn(1,:)
        if (ni.ge.2) then
          do j=2,ni
            write(20,1003) bcurvn(j,:)
          enddo
        endif
        write(20,1002) 'lambda:', xlam
      endif
!
!.....Initialize the cplex environment when this is the first outer and
!.....inner loop
      if (outerloop.eq.1.and.lloop.eq.0) then
        status = 0
        call cplex_init(status)

        if (status.eq.-10) then
          stop ' Cplex is not correctly installed.'
        endif
        
        if (status.ne.0.and.status.ne.-10) then
          call cplex_stop(status)
          stop ' initialization of CPLEX failed'
        endif
      endif
!      
!.....construct the primal variable bounds
      do i=1,nprimal
        bl(i)=x_l(i)-xprimal(i)
        bu(i)=x_u(i)-xprimal(i)
      enddo
!
!.....form qpdata which contains the diagonal terms of the second
!.....derivative to x of the lagrange function

      do i=1,nprimal
        qpdata(i) = acurvn(i)
       do j=1,ni
          bji   = bcurvn(j,i)
          xlamj = xlam(j)
          if (bji.gt.0.d0.and.xlamj.gt.0.d0)then
            qpdata(i) =  qpdata(i) + bji*xlamj
          endif
        enddo
      enddo



      do i = 1, ni
      do j = 1, nprimal
!           write (*,*) gc(i,j)
 
          enddo
          enddo   
!.....Solve the QP problem by putting it in the cplex workspace and then
!.....call the cplex QP solver\
      loops = outerloop + lloop
      status = 0
      call cplex_solve_sqp(status, nprimal, ni, c, gf, gc, bl, bu,
     &                     qpdata, x, xlam, loops, cplex_method,
     &                     cplex_feas, eqn)      
!
      if (status .ne. 0) then
        call cplex_stop(status)
      stop ' Solving of QP with CPLEX failed'
      endif
!
!.....update the primal variables
      do i=1,nprimal
        xprimal(i) = min(x_u(i),max(x_l(i),xprimal(i)+x(i)))
      enddo
!
!.....Get the approximate function value
      f_a = f
      do i=1,nprimal
        f_a = f_a + gf(i)*x(i) + 0.5*acurvn(i)*x(i)**2
      enddo
!
!.....Update the lagrange multipliers
      do j=1,ni
        xlam(j) = dabs(xlam(j))
      enddo

!.....write some iteration stuff for new solution
      if (debug) then
        write(20,1002) 'xstep:', x
        write(20,1002) 'xprimelk+1:',xprimal
        write(20,1002) 'approx. obj:', f_a
      endif
!
!.....message equals status of the cplex calls
      message = status
! 
      return
 1000 format (100('='),/,'outerloop, innerloop: ',i5,',',i5)
 1001 format (/,'n: ',i7,/,'ni:',i7)
 1002 format (/,a30,/,10es18.4)
 1003 format (/,10es18.4)
 
      end
      
       subroutine drive_cplexs_old (nprimal,xprimal,ni,nq,f_a,c_a,
     &                         iact,nact,xlam,x_h,
     &                         acurv,bcurv,amult,bmult,x_l,x_u,
     &                         f,c,gf,gc,ksubiter,
     &                         xkkt,flam,cstage0,ostring,xlam_h,
     &                         nsx,ictrl,lctrl,rctrl,cctrl,message,yhi,
     &                         outerloop,lloop,nnz,Acol,Aptr,
     &                         storind,delta_t,eqn,
     &                         iuser,luser,cuser,ruser,nnzh)!
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
      integer          ni,nq,nprimal,ksubiter,m,outerloop,nnzh
      integer          nnz, Aptr(ni+1), Acol(nnz), storind(nnz)
      integer          nsx,nact,iact(nact),message
      integer          loops, lloop, i, j, ij, i1
      integer          status
      double precision xprimal(nprimal),flam,bu(nprimal),bl(nprimal)
      double precision x_h(nprimal),x_l(nprimal),x_u(nprimal)
      double precision qpdata(nprimal), delta_t, qpdatas
      double precision t_c1, t_c2, SAOi_seconds
      double precision xlam(ni),yhi
      double precision acurv,bcurv(ni),xlam_h(ni),bcurvn(ni,nprimal)
      double precision f,c(ni),gf(nprimal),gc(nnz)
      double precision f_a,c_a(ni)
      double precision xkkt, xlamj, bji, x(nprimal), gcij
      character*90     ostring
      double precision amult, bmult(nprimal+ni),dl
      double precision Hval(nnzh)
      integer          Hptr(nprimal+1),Hcol(nnzh),Hne
!.....Get global option variables
      include          'ctrl_get.inc'
!
      if (nprimal.gt.nmax) stop ' nprimal.gt.nmax in drive_cplex.f'
!
!.....write input data to a file by first opening the file and then 
      if (debug) then
        if (outerloop.eq.1.and.lloop.eq.0) then
          open (20,file='cplex_loops.out')
        endif
        write(20,1000) outerloop, lloop
        write(20,1001) nprimal, ni
        write(20,1002) 'xstart: ', xprimal
        write(20,1002) 'x_low:', x_l
        write(20,1002) 'x_heigh:', x_u
        write(20,1002) 'fk: ', f  
        write(20,1002) 'ck: ', c
        write(20,1002) 'dFk:', gf
!        write(20,1002) 'dCk:', gc(1,:)
!         if (ni.ge.2) then
!           do j=2,ni
! !           write(20,1003) gc(j,:)
!           enddo
!         endif
!         write(20,1002) 'Qf:', acurvn
!         write(20,1002) 'Qc:', bcurvn(1,:)
!         if (ni.ge.2) then
!           do j=2,ni
!             write(20,1003) bcurvn(j,:)
!           enddo
!         endif
!         write(20,1002) 'lambda:', xlam
      endif
!
!.....Initialize the cplex environment when this is the first outer and
!.....inner loop
      if (outerloop.eq.1.and.lloop.eq.0) then
        status = 0
        call cplex_init(status)

        if (status.eq.-10) then
          stop ' Cplex is not correctly installed.'
        endif
        
        if (status.ne.0.and.status.ne.-10) then
          call cplex_stop(status)
          stop ' initialization of CPLEX failed'
        endif
      endif

      
!..updated Wed 15 Nov 2017 14:24:14 SAST
!..see equation 5 in Etman2012 or equation 21 and 22 in Munro2017a
!..construct the primal variable bounds

      dl=1.0d0       !init. move-limit
!      
!.....determine if nprimal does not exceed nmax.
!.....construct the primal variable bounds
      do i=1,nprimal
        bl(i)=max(-dl*(x_u(i)-x_l(i)),(x_l(i)-xprimal(i)))
        bu(i)=min(dl*(x_u(i)-x_l(i)),(x_u(i)-xprimal(i)))
            if(x_l(i).eq.x_u(i)) then
                bl(i)=0.0d0
                bu(i)=0.0d0
            endif
      enddo
      if (approx_f.eq.1.and.approx_c.eq.1) then
        qpdatas = max(atol1,acurv*amult)!atol1=1.0e-05
        do j = 1,ni
          if (bcurv(j).gt.0.d0.and.xlam(j).gt.0.d0) then!btol1=-1000
            qpdatas = qpdatas + max(btol1,bmult(j)*bcurv(j))*xlam(j)
          endif
        enddo
        do i = 1,nprimal
          qpdata(i) = qpdatas! all qpdata is the same!
        enddo
      elseif (approx_f.gt.1.or.approx_c.gt.1) then
        do i = 1,nprimal
           if (approx_f.eq.4) then
               qpdata(i) = max(atol1,2.d0*amult/x_h(i)*dabs(gf(i)))
           elseif (approx_f.eq.1) then
                qpdata(i) = max(atol1, acurv*amult)   
           endif
        enddo
        if (approx_f.eq.100.and.approx_c.eq.100) then
          call Sep_gg(xprimal,nprimal,ni,xlam,qpdata) !acurvn=0.d0!
!           qpdata = 2.d0*qpdata
        endif
        if (approx_f.lt.100.and.approx_c.lt.100) then
        ij = 0;
          do j = 1, ni
            xlamj = xlam(j)
            do i1 = 1, Aptr(j+1) - Aptr(j)
              ij = ij + 1
              i = Acol(ij)
              gcij = gc(ij)
              if(approx_c.eq.1) then
                bji = max(btol1, bmult(j)*bcurv(j))
              elseif(approx_c.eq.4) then
                bji = max(btol1,2.d0*bmult(j)/x_h(i)*dabs(gcij))
              endif
              if (bji.gt.0.d0.and.xlamj.gt.0.d0)then
                qpdata(i) = qpdata(i) + bji*xlam(j)
              endif
            enddo
          enddo
        endif
      endif
      
      if (approx_f.eq.101.and.approx_c.eq.101) then
        call Fem_gg(xprimal,nprimal,ni,xlam,nnzh,Hval,Hptr,Hcol,Hne)
        do i = 1, nprimal
          Hval(i) = 1.0d2*Hval(i)
        enddo
      endif

!.....Solve the QP problem by putting it in the cplex workspace and then
!.....call the cplex QP solver\
      loops = outerloop + lloop
      status = 0

      t_c2 = SAOi_seconds();
      if(approx_f.eq.101.and.approx_c.eq.101) then
        call cplex_solve_sqps(status,nprimal, ni, c, gf, gc, bl, bu,
     &                     x, xlam, loops, cplex_method,
     &                     cplex_feas,nnz,Acol,Aptr,nnzh,Hcol,Hptr,
     &                     Hval,storind,delta_t,eqn)!,qpdata
      else
        call cplex_solve_sqps_sep(status,nprimal, ni, c, gf, gc, bl, bu,
     &                     qpdata, x, xlam, loops, cplex_method,
     &                     cplex_feas, nnz, Acol, Aptr, storind,delta_t,
     &                     eqn)
      endif

      if (status .ne. 0) then
        call cplex_stop(status)
      stop ' Solving of QP with CPLEX failed'
      endif
!

!.....update the primal variables
      do i=1,nprimal
        xprimal(i) = min(x_u(i),max(x_l(i),xprimal(i)+x(i)))
      enddo
!
!.....Get the approximate function value
      f_a = f
      if (approx_f.eq.1) then
      do i=1,nprimal
        f_a = f_a + gf(i)*x(i) + 
     &         max(atol1,acurv*amult)/2.d0*x(i)**2
        enddo
      elseif (approx_f.eq.4) then
      do i = 1,nprimal
      f_a = f_a + gf(i)*x(i) +
     &      max(atol1,2.d0/x_h(i)*dabs(gf(i)))/2.d0*x(i)**2
      enddo
      elseif (approx_f.eq.100) then
      f_a = f
      do i=1,nprimal
        f_a = f_a + gf(i)*x(i)+0.5d0*qpdata(i)*x(i)**2.d0!
!      &         +max(atol1,acurv*amult)/2.d0*x(i)**2
      enddo
      endif

!
!.....Update the lagrange multipliers
      do j=1,ni
        xlam(j) = dabs(xlam(j))!-xlam(j)!!
      enddo
!
!.....write some iteration stuff for new solution
      if (debug) then
        write(20,1002) 'xstep:', x
        write(20,1002) 'xprimelk+1:',xprimal
        write(20,1002) 'approx. obj:', f_a
      endif
!
!.....message equals status of the cplex calls
      message = status

! 
!        write(*,*) (t_c2 - t_c1)
      return
 1000 format (100('='),/,'outerloop, innerloop: ',i5,',',i5)
 1001 format (/,'n: ',i7,/,'ni:',i7)
 1002 format (/,a30,/,10es18.4)
 1003 format (/,10es18.4)
 
      end
      

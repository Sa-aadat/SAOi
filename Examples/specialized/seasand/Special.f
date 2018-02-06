! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Specialized instructions, currently for local stress
! SAOi: constrained topology optimization
! SAOi:


      subroutine SAOi_special (n, ni, ne, x, f, c, h, iuser, luser, 
     &                         cuser, ruser, eqn, lin, ictrl, lctrl, 
     &                         rctrl, cctrl, x_lo, x_up, outerloop,
     &                         finished, kkt, violation, nactlo,
     &                         nacthi, nactc, delxnormc, delxnormi,
     &                         time)
!----------------------------------------------------------------------!
!                                                                      !
! Optional specialized instructions may be performed here. This        !
!     routine is called at the end of each OUTER iteration, but        !
!     also at iteration 0. If the algorithm has converged or           !
!     a stop is enforced, finished = .true.                            !
!                                                                      !
! Specialized.out is on unit 15                                        !
!                                                                      !
! Please see the users manual for type declarations and comments       !
!                                                                      !
!----------------------------------------------------------------------!
      implicit         none
      include          'ctrl.h'
      logical          eqn(*), lin(*), finished
      integer          i, j, n, ni, ne, outerloop
      integer          nactlo, nacthi, nactc
      double precision f, x(n), c(ni), h(ne), x_lo(n), x_up(n), time
      double precision kkt, violation, phiBW, delxnormc, delxnormi
      include          'ctrl_get.inc'

! In general, do nothing but return emptyhanded
      if (.not.special) return

! Else, do something clever here, like calculating the black-and-white
!     fraction in topology optimization, and implement stress
!     relaxation for locally constrained topology optimization,
!     etc.
!
      ! the code comes here !
!
      return
      end subroutine SAOi_special
!----------------------------------------------------------------------!
      subroutine filter1(nelx,nely,rmin,x,dc,dcold,n)
      implicit         none
      integer          nelx,nely,n,icol,irow,ii,jj,icount
      integer          imaxk,imaxl,imink,iminl,kk,ll
      double precision x(n),dc(n),dcold(n),rmin,sum,xfac
      double precision xmat(nely,nelx),xii,xjj,xkk,xll
      double precision dcmat(nely,nelx),dcoldmat(nely,nelx)
      
! pack x, dc and dcold in matrix form
      icount=0
      do icol=1,nelx
        do irow=1,nely
          icount=icount+1
          xmat(irow,icol)  = x(icount)
          dcoldmat(irow,icol) = dcold(icount)
          dcmat(irow,icol)    = 0.d0
        enddo
      enddo
 
! now filter sensitivities
      do ii=1,nelx
        do jj=1,nely
          sum=0.d0
          imink=nint(max(dble(ii)-rmin,1.d0))
          imaxk=nint(min(dble(ii)+rmin,dble(nelx)))
          iminl=nint(max(dble(jj)-rmin,1.d0))
          imaxl=nint(min(dble(jj)+rmin,dble(nely)))
          do kk=imink,imaxk
            do ll=iminl,imaxl
              xii = dble(ii)
              xjj = dble(jj)
              xkk = dble(kk)
              xll = dble(ll)
              xfac = rmin - dsqrt((xii-xkk)**2.d0 + (xjj-xll)**2.d0)
              sum = sum + dmax1(0.d0,xfac)
              dcmat(jj,ii) = dcmat(jj,ii) +
     &        (dmax1(0.d0,xfac)*xmat(ll,kk)*dcoldmat(ll,kk))
            enddo
          enddo
          dcmat(jj,ii) = dcmat(jj,ii)/(xmat(jj,ii)*sum)
        enddo
      enddo
 
! repack xmat, dcmat and dcoldmat into vector form
      icount=0
      do icol=1,nelx
        do irow=1,nely
            icount=icount+1
            dc(icount) = dcmat(irow,icol)
        enddo
      enddo
!
      return
      end
!
!
!
!
!
!
!
!
!
!
!
        subroutine clean(nprimal,ndual,nnz,xprimal,xprimalc,
     &  xdual,xdualc,nobj,con,jac,hes,aptr,aptrc,acol,xl,xu,
     &  eqn,nprimalc,ndualc,nnzc,bc)

        implicit none
        integer             i,j,k1,k2,k3,k4,k5,k6
        integer             nprimal,ndual,nnz
        integer             nprimalc,ndualc,nnzc
        integer             aptr(ndual+1),acol(nnz)
        integer             bc(nprimal)
        double precision    xprimal(nprimal),xdual(ndual)
        double precision    xprimalc(nprimal),xdualc(ndual)
        double precision    nobj(nprimal)
        double precision    con(nprimal),jac(nnz),hes(nprimal)
        double precision    xl(nprimal),xu(nprimal)
        real                eqn(ndual)
        integer             aptrc(ndual+1)

        do i=1,nprimal
        bc(i) = 0
        if(xl(i).eq.xu(i)) then
        bc(i)=1
        endif


!       if(xprimal(i).eq.xl(i).or.xprimal(i).eq.xu(i)) then
!       bc(i)=1
!       endif

        enddo
        
        ndualc=0
        k1=0
        k2=0
        k3=0
        k4=1
        k5=1
        k6=2
        nnzc=0
        aptrc(1)=1
        do i = 1,ndual
        if(bc(i).eq.0) then
        con(k4)=con(i)
        eqn(k4)=eqn(i)
        xdualc(k4)=xdual(i)
        k4=k4+1
        k5=k5+1
        aptrc(k5)=aptrc(k5-1)+aptr(k6)-aptr(k6-1)
        ndualc=ndualc+1
        endif
        k6=k6+1
        do j = 1,aptr(i+1)-aptr(i)
        k3=k3+1
        k2=0
        if(bc(i).eq.0) then
        nnzc=nnzc+1
        k2=k2+1
        k1=k1+1
        jac(k1)=jac(k3)
        acol(k1)=acol(k3)
        endif
        enddo
        enddo

        endsubroutine
        subroutine connectroutine(prob,mult,nsta,ndual,
     &             connect,rowcounter,nodes,ground)
    
        implicit            none
        integer             ncc,nex,ney
        integer             r,s,e,i,j,k
        integer             mult,nsta,ndual,nodec
        integer             connect(nsta,42)
        integer             rowcounter(ndual)
        integer             nodes(nsta),prob,nmat
        integer             ground
                
        nex=15*mult
        ney=5*mult
        if(ground.eq.3) then
        nex=25*mult
        ney=25*mult
        endif
        nmat = nex*ney

        ! Safety
        do i = 1,nsta
        nodes(i) = 0
        enddo
        do i = 1,nsta
        do j = 1,42
        connect(i,j) = 0
        enddo
        enddo
        do i = 1,nsta
        rowcounter(i) = 17
        enddo
        
        if(prob.eq.1) then
        do i = nsta+1,ndual
        rowcounter(i) = 9
        enddo
        else if(prob.eq.2) then
        do i = nsta+1,ndual
        rowcounter(i) = nex*ney
        enddo
        else if(prob.eq.3) then
        rowcounter(nsta+1) = nex*ney
        do i = nsta+2,ndual
        rowcounter(i) = 2
        enddo
        else if(prob.eq.4) then
        do i = nsta+1,nsta+nmat
        rowcounter(i) = 9
        enddo
        do i = nsta+nmat+1,ndual
        rowcounter(i) = 2
        enddo
        else if(prob.eq.99) then
        do i = nsta+1,ndual
        rowcounter(i) = nex*ney
        enddo
        elseif(prob.eq.100) then
        do i = nsta+1,ndual
        rowcounter(i) = 9
        enddo
        endif

        !Loop 0 
        nodec = 2*(2*nex+1)+3 
        s = nodec 
        e = 0 
        ncc = 1 
        do j = 1,ney
        do i =1,2*(nex-1) 
      
        do k = 1,10
        connect(nodec,k) = s - 2*(2*nex+1) + ncc*(2) + k - 5
        connect(nodec,16+k) = connect(nodec,k) + 2*(2*nex+1) + 2*(nex+1)
        enddo
        do k = 11,16
        connect(nodec,k) =  s - 2 + (k-11)
        enddo
        rowcounter(nodec) = 28
        nodec = nodec + 1
        e = e + 1
        
        if (e.eq.2) then
        s = s + 2
        ncc = ncc + 1
        e = 0
        endif
        enddo
        ncc = 1
        nodec = nodec + 1+ 2*(2*nex+1) + 3
        s = nodec
        enddo
                       
        !Loop X
        nodec = 2*(2*nex+1)+2*(nex+1) + 3
        s = nodec 
        e = 0 
        ncc = 1 
        do j = 1,(ney-1)
        do i =1,(2*nex) 
        do k = 1,6
        connect(nodec,k) = s - 2 - 2*(2*nex+1) - 2*(nex+1) + k - 1
        connect(nodec,20+k) = connect(nodec,k) + 4*(2*nex+1) + 4*(nex+1)
        enddo
        do k = 7,10
        connect(nodec,k) = s - 2*(nex+1) - 2*ncc + k - 7
        connect(nodec,10+k)=connect(nodec,k)+2*(2*nex+1)+2*(nex+1)
        enddo
        do k = 11,16
        connect(nodec,k) = s - 2 + k - 11
        enddo
        rowcounter(nodec) = 28
        nodec = nodec + 1
        e = e + 1
        r = r + 1
        if (e.eq.2) then
        s = s + 4
        ncc = ncc + 1
        e = 0
        nodec = nodec + 2
        endif
        enddo
        ncc = 1
        nodec = nodec + 1 + 2*(nex+1) + 1
        s = nodec
        enddo
        
        !Loop WL
        nodec = 2*(2*nex+1)+2*(nex+1) + 1
        do i = 1,(ney-1)
        do k = 1,26
        connect(nodec,k) = connect(nodec+2,k)
        connect(nodec+1,k) = connect(nodec+3,k)
        enddo

        rowcounter(nodec) = 28
        nodec = nodec + 1
        rowcounter(nodec) = 28
        nodec = nodec + 2*(2*nex+1) + 2*(nex+1) - 1
        enddo 
                  
        !Loop WR
        nodec = 4*(2*nex+1)+2*(nex+1) - 1 
        do i = 1,(ney-1)
        do k = 1,26
        connect(nodec,k) = connect(nodec-2,k)
        connect(nodec+1,k) = connect(nodec-1,k)
        enddo
            
        rowcounter(nodec) = 28
        nodec = nodec + 1
        rowcounter(nodec) = 28
        nodec = nodec + 2*(2*nex+1) + 2*(nex+1) - 1
        enddo  

        !Loop WU 	
        nodec = 5
        ncc = 1
        do i = 1,(nex-1)
        do k = 1,26
        connect(nodec,k) = connect(nodec + 2*(2*nex+1) - 2*(ncc),k)
        connect(nodec+1,k) = connect(nodec + 2*(2*nex+1) - 2*(ncc),k)
        enddo
            
        rowcounter(nodec) = 28
        nodec = nodec+1
        rowcounter(nodec) = 28
        nodec = nodec + 3
        ncc = ncc + 1
        enddo        
         
        !Loop WB
        ncc = 1
        nodec = nsta -2*(2*nex+1) + 5
        do i = 1,(nex-1)
        do k = 1,26
        connect(nodec,k)=connect(nodec-4*(ncc)-2*(nex+1) + 2*(ncc),k)
        connect(nodec+1,k)=connect(nodec-4*(ncc)-2*(nex+1) + 2*(ncc),k)
        enddo
        rowcounter(nodec) = 28
        nodec = nodec+1
        rowcounter(nodec) = 28
        nodec = nodec + 3
        ncc = ncc + 1
        enddo      
                     
        !Loop Y
        nodec = 2*(2*nex+1) + 2*(nex+1) + 5
        s = nodec
        ncc = 1
        e = 0
        do i = 1,(ney-1)
        do j = 1,(nex-1)
        do k = 1, 26
        connect(nodec,k) = connect(s - j*(2) -2*(nex+1), k)
        connect(nodec+1,k) = connect(s - j*(2) -2*(nex+1), k)
        enddo
        do k = 27, 42
        connect(nodec,k) = connect(nodec+2*(2*nex+1)-2-2*(j-1),k-16)
        connect(nodec+1,k) = connect(nodec+2*(2*nex+1)-2-2*(j-1),k-16)
        enddo
              
        rowcounter(nodec) = 46
        nodec = nodec + 1
        rowcounter(nodec) = 46
        nodec = nodec + 3
        s = s + 4
        enddo
        nodec = nodec + 6 + 2*(nex+1)
        s = nodec
        enddo 

        return
        end subroutine
!         subroutine drivecplex(outerloop,nprimal,ndual,nnz,
!      &  xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
!      &  xl,xu,dl,du,nsta,dx,feasible)
! 
!         implicit            none
!         integer             i,j,k
!         logical             eqn(ndual)
!         integer             nprimal,ndual,nnz,aptr(ndual+1)
!         double precision    xprimal(nprimal),xdual(ndual)
!         double precision    obj,nobj(nprimal),con(ndual)
!         double precision    jac(nnz),hes(nprimal),dnorm
!         double precision    xl(nprimal),xu(nprimal),biglam
!         double precision    qpbl(nprimal),qpbu(nprimal)
!         double precision    delta(nprimal),lambda(ndual)
!         double precision    dx(nprimal),dl,du,eps
!         integer             acol(nnz),outerloop,nsta
!         integer             status
!         logical             feasible
!         integer             bnprimal,bndual,bnnz
! 
! !       double precision,dimension(:),allocatable::bnobj
! !       double precision,dimension(:),allocatable::bjac
! !       double precision,dimension(:),allocatable::bqpbl
! !       double precision,dimension(:),allocatable::bqpbu
! !       double precision,dimension(:),allocatable::bdelta
! !       integer,dimension(:),allocatable::baptr
! !       integer,dimension(:),allocatable::bacol
! 
! !       bnprimal=nprimal+(ndual-nsta)
! !       bndual=ndual
! !       bnnz=nnz+(ndual-nsta)
! 
! !       allocate(bnobj(bnprimal))
! !       allocate(bdelta(bnprimal))
! !       allocate(bjac(bnnz))
! !       allocate(bqpbl(bnprimal))
! !       allocate(bqpbu(bnprimal))
! !       allocate(baptr(bndual+1))
! !       allocate(bacol(bnnz))
! 
! !       initialize the cplex environment
!         if (outerloop.eq.1) then
!         status = 0
!         if (status.eq.-10) then
!         stop ' Cplex is not correctly installed.'
!         endif
!         if (status.ne.0.and.status.ne.-10) then
!         stop ' initialization of CPLEX failed'
!         endif
!         endif
! 
! !       do i=1,nprimal
! !       bnobj(i)=0.0d0
! !       bdelta(i)=0.0d0
! !       bqpbl(i)=max(-dl*(xu(i)-xl(i)),(xl(i)-xprimal(i)))
! !       bqpbu(i)=min(dl*(xu(i)-xl(i)),(xu(i)-xprimal(i)))
! !       if(xl(i).eq.xu(i)) then
! !       bqpbl(i)=0.0d0
! !       bqpbu(i)=0.0d0
! !       endif
! !       enddo
! !       do i=nprimal+1,bnprimal
! !       bnobj(i)=1.0d0
! !       bdelta(i)=0.0d0
! !       bqpbl(i)=0.0d0
! !       bqpbu(i)=1d6
! !       enddo
! 
! !       do i=1,ndual+1
! !       baptr(i)
! !       enddo
! 
! !       k=1
! !       do i=1,ndual
! !       do j=1,aptr(i+1)-aptr(i)
! 
! !       k=k+1
! !       enddo
! !       enddo
! !       write(*,*) k, nnz
! 
! 
! !       call cplex_init(status)
! !       call cplex_solve_sqps(status,bnprimal,bndual,bnnz,eqn,
! !    &  bdelta,lambda,bnobj,con,bjac,bhes,baptr,bacol,
! !    &  bqpbl,bqpbu,1)
! !       call cplex_stop(status)
! 
! !       construct the primal variable bounds
!         do i=1,nprimal
!         delta(i)=0.0d0
!         dx(i)=0.0d0
! !       if(i.gt.nsta) then
!         qpbl(i)=max(-dl*(xu(i)-xl(i)),(xl(i)-xprimal(i)))
!         qpbu(i)=min(dl*(xu(i)-xl(i)),(xu(i)-xprimal(i)))
! !       else
! !       qpbl(i)=-du*dl
! !       qpbu(i)=du*dl
!         if(xl(i).eq.xu(i)) then
!         qpbl(i)=0.0d0
!         qpbu(i)=0.0d0
!         endif
! !       endif
!         enddo
!         do j=1,ndual
!         lambda(j) = 0.0d0
!         enddo
! 
!         call cplex_init(status)
!         call cplex_solve_sqps(status,nprimal,ndual,nnz,eqn,
!      &  delta,lambda,nobj,con,jac,hes,aptr,acol,qpbl,qpbu,2)
!         feasible=.true.
!         if (status.eq.1) then
!         stop ' Solving of QP with CPLEX failed, 1'
!         elseif (status.eq.2) then
!         feasible=.false.
!         endif
!         call cplex_stop(status)
!         
!         if(feasible) then
! !       update the primal variables
!         do i=1,nprimal
!         dx(i) = delta(i)
! !       xprimal(i)=min(xu(i),max(xl(i),xprimal(i)+delta(i)))
!         enddo
!         biglam = 1d16
! !       update the lagrange multipliers
!         do j=1,ndual
!         xdual(j)= max(-biglam,min(-lambda(j),biglam))
!         enddo
!         endif
! !       deallocate(bnobj)
! !       deallocate(bdelta)
! !       deallocate(bjac)
! !       deallocate(bqpbl)
! !       deallocate(bqpbu)
! !       deallocate(baptr)
! !       deallocate(bacol)
!         return
!         end subroutine
!         subroutine drivecplex_restoration(outerloop,nprimal,ndual,nnz, 
!      &  xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
!      &  xl,xu,dl,du,nsta,dx,feasible)
! !       call drivecplex_restoration(outerloop,nprimal,ndual,nnz, 
! !    &  xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
! !    &  xl,xu,dl,du,nsta,dx,feasible)
! 
!         implicit            none
!         integer             i,j
!         logical             eqn(ndual)
!         integer             nprimal,ndual,nnz,aptr(ndual+1)
!         double precision    xprimal(nprimal),xdual(ndual)
!         double precision    obj,nobj(nprimal),con(ndual)
!         double precision    jac(nnz),hes(nprimal),dnorm
!         double precision    xl(nprimal),xu(nprimal),biglam
!         double precision    qpbl(nprimal),qpbu(nprimal)
!         double precision    delta(nprimal),lambda(ndual)
!         double precision    dx(nprimal),dl,du
!         integer             acol(nnz),outerloop,nsta
!         integer             status
!         logical             feasible
! 
! !       initialize the cplex environment
!         if (outerloop.eq.1) then
!         status = 0
!         call cplex_init(status)
! 
!         if (status.eq.-10) then
!         stop ' Cplex is not correctly installed.'
!         endif
!         
!         if (status.ne.0.and.status.ne.-10) then
!         call cplex_stop(status)
!         stop ' initialization of CPLEX failed'
!         endif
!         endif
! 
! !       construct the primal variable bounds
!         do i=1,nprimal
!         delta(i)=0.0d0
!         dx(i)=0.0d0
! !       if(i.gt.nsta) then
!         qpbl(i)=max(-dl*(xu(i)-xl(i)),(xl(i)-xprimal(i)))
!         qpbu(i)=min(dl*(xu(i)-xl(i)),(xu(i)-xprimal(i)))
! !       else
! !       qpbl(i)=-du*dl
! !       qpbu(i)=du*dl
!         if(xl(i).eq.xu(i)) then
!         qpbl(i)=0.0d0
!         qpbu(i)=0.0d0
!         endif
! !       endif
!         enddo
!         do j=1,ndual
!         lambda(j) = 0.0d0
!         enddo
! 
! 
!         call cplex_solve_sqps(status,nprimal,ndual,nnz,eqn,
!      &  delta,lambda,nobj,con,jac,hes,aptr,acol,qpbl,qpbu,1)
! 
!         feasible=.true.
!         if (status.eq.1) then
!         call cplex_stop(status)
!         stop ' Solving of QP with CPLEX failed, 1'
!         elseif (status.eq.2) then
! !       call cplex_stop(status)
!         feasible=.false.
!         endif
!         
!         if(feasible) then
! !       update the primal variables
!         do i=1,nprimal
!         dx(i) = delta(i)
! !       xprimal(i)=min(xu(i),max(xl(i),xprimal(i)+delta(i)))
!         enddo
! 
!         biglam = 1d16
! !       update the lagrange multipliers
!         do j=1,ndual
!         xdual(j)= max(-biglam,min(-lambda(j),biglam))
!         enddo
!         endif
! 
!         return
!         end subroutine
        subroutine filter(nprimal,obj,viol,objh,violh,nobjh,hesh,
     &  dx,dl,dl0,innerloop,rem,filterlist,next,hmult,fstep,hstep,
     &  consv,dlm)

        implicit none
        logical             cond1,cond2,next
        logical             membrane,feasible
        integer             nprimal,innerloop,rem,i,j
        double precision    gama,beta,sigma,dl,dl0,kapa,zeta
        double precision    obj,viol,objh,violh,dx(nprimal)
        double precision    hesh(nprimal),nobjh(nprimal)
        double precision    filterlist(1000000,2),dlm
        double precision    delta_f, delta_q, hmult
        integer             fstep,hstep, consv 

        gama= 1.0d-7
        beta= 1.0d0-gama
        sigma = 1d-2
        kapa=1d-3
        zeta=2.0d0

        membrane=.false.
        cond1=.true.
        cond2=.true.

        do i = 1,rem
        if(viol.gt.beta*filterlist(i,2)) then
        cond1=.false.
        endif
        if(obj+gama*viol.gt.filterlist(i,1)) then
        cond2=.false.
        endif
        enddo


        if(cond1) then
        membrane=.true.
        endif
        if(cond2) then
        membrane=.true.
        endif
        if(obj.lt.0.0d0) then
        membrane=.false.
        endif

        if(membrane) then

        delta_q = 0.0d0
        delta_f = objh-obj
        do j = 1,nprimal
        delta_q=delta_q - dx(j)*nobjh(j)
     &  -dx(j)**2.0d0*hesh(j)/2.0d0
        enddo
           
        !if both feasible descent conditions are `false' 
        if(delta_f.lt.sigma*delta_q.and.
     &  delta_q.gt.kapa*violh**zeta) then
                if(consv.eq.1) then
                hmult=hmult*2.0d0
                else
                dl=dl/2.0d0
                endif
        next=.false.
        else !if one of them is true...

        !if 
        if(delta_q.gt.kapa*violh**zeta) then
         fstep=fstep+1
                if(consv.eq.1) then
                hmult=hmult/2.0d0
                hmult=max(hmult,1.0d0)
                else
                dl=dl*2d0
!               dl=min(dl0,dl)
                dl=min(dlm,dl)
                endif
         next=.true.
        else
         hstep=hstep+1
         rem=rem+1
                if(rem.gt.1000000) then
                write(*,*) 'WARNING: The filter is full'
                endif
         filterlist(rem,1) = objh
         filterlist(rem,2) = violh
                if(consv.eq.1) then
                hmult=hmult/2.0d0
                hmult=max(hmult,1.0d0)
                else
                dl=dl*2d0
!               dl=min(dl0,dl)
                dl=min(dlm,dl)
                endif
         next=.true.
        endif


        endif
        else
        if(consv.eq.1) then
        hmult=hmult*2.0d0
        else
        dl=dl/2.0d0
        endif
        next=.false.
        endif

        if(consv.eq.1) then
        do i =1,nprimal
        hesh(i)=hesh(i)*hmult
        enddo
        endif

        end subroutine
        subroutine forceroutine(mult,prob,ground,fvec,fxy)
        implicit            none
        integer             i,mult,ground,prob
        integer             fxy(2,2),nex,ney
        double precision    fvec(16),fx,l,fy
        
        nex = 15*mult
        ney = 5*mult
        l  = 3.0d0/15.0d0/mult

        if(ground.eq.1) then
        fx = 10.0d0
        fxy(1,1) =    nex/2 - 15*mult/10 + 1!2    !13!x-coor-1
        fxy(1,2) =    1!1    !y-coor-1
        fxy(2,1) =    nex/2 + 15*mult/10!2         !18!x-coor-2
        fxy(2,2) =    1!1    !y-coor-2

        do i =1,16
        fvec(i) = 0.0d0
        enddo
        fx = fx/2.0d0*l

        fvec(1) = 1.d0*fx/3.0d0
        fvec(3) = 1.d0*fx/3.0d0
        fvec(9) = 4.d0*fx/3.0d0

        elseif(ground.eq.2) then
        fx = 0.0d0
        fy = -1.0d0
        fxy(1,1) =    1!2    !13!x-coor-1
        fxy(1,2) =    ney!1    !y-coor-1
        fxy(2,1) =    1!2         !18!x-coor-2
        fxy(2,2) =    ney!1    !y-coor-2

        do i =1,16
        fvec(i) = 0.0d0
        enddo
        fx = fx/2.0d0*l
        fy = fy/2.0d0*l

!       fvec(7) = 1.d0*fx/3.0d0
        fvec(8) = -1.0d0!1.d0*fy/3.0d0
!       fvec(5) = 1.d0*fx/3.0d0
!       fvec(6) = 1.d0*fy/3.0d0
!       fvec(13) = 4.d0*fx/3.0d0
!       fvec(14) = 4.d0*fy/3.0d0
        
        elseif(ground.eq.3) then
        nex = 25*mult
        ney = 25*mult
        l  = 10.0d0/25.0d0/mult
        fx = 0.0d0
        fy = -2.0d0
        fxy(1,1) =    nex!2    !13!x-coor-1
        fxy(1,2) =    4*ney/10/2 - mult/2 + 1!1    !y-coor-1
        fxy(2,1) =    nex!2         !18!x-coor-2
        fxy(2,2) =    4*ney/10/2 + mult/2!8*mult!1    !y-coor-2

        do i =1,16
        fvec(i) = 0.0d0
        enddo
        fx = fx/2.0d0*l
        fy = fy/2.0d0*l
        fvec(4) = 1.d0*fy/3.0d0
        fvec(6) = 1.d0*fy/3.0d0
        fvec(12) = 4.d0*fy/3.0d0

        endif

        return
        
        end subroutine 
        
        
!         program global
! 
        subroutine globalD1 (idrive,nprimal,nmat,nsta,ndual,nnz,
     &      du,ground,prob)
        implicit none
        logical feasible
        integer mult,ground,approx,prob,solver,trust,idrive
        double precision filt,dl0,dlc,du,dlm
        integer exprm,nex,ney,nsta,nmat,bs,nprimal
        integer ndual,nnz,rs,glo,i,r,loops,consv
        integer warm, conti, scal, conv
        double precision  sorm, filt0, x0
        real t,tt,seconds,approxseconds
        character*12 :: filename
        character (len=100) :: temp

        double precision f_k,f_star,h_k,h_star,pr
        double precision,dimension(:),allocatable::x_k
        double precision,dimension(:),allocatable::x_star
 
        !!!
        include 'setsea.h'
        !!!

        filt=filt0

        !t=seconds()
        
        call init (mult,nex,ney,ground,prob,nsta,nmat,
     &  bs,nprimal,ndual,nnz,approx,filt,solver,
     &  dl0,dlc,du,trust,consv,rs)

! 
!         do glo = 1,exprm
! 
!         filt=filt0
! 
!         allocate(x_k(nprimal))
!         call seasand(mult,ground,approx,prob,filt,solver,
!      &  trust,consv,rs,dl0,dlc,f_k,h_k,feasible,x_k,
!      &  loops,glo,warm,approxseconds,conti,dlm,conv,x0)
! 
! !       write history
!         filename = 'global.dat'
!         if(glo.eq.1) then    
!         open(unit=9999, file=filename, action="write",
!      &  status="replace")
!         else
!         open(unit=9999, file=filename, action="write",
!      &  access="append")
!         endif
!         tt=seconds()
!         write(9999,*)glo,loops,f_k,h_k,feasible,tt-t
!         close(9999)
! 
!         if(glo.lt.10) then
!         WRITE(filename,'(a,i1.1,a)') "xsta",glo,".dat"  
!         elseif(glo.lt.100) then
!         WRITE(filename,'(a,i2.2,a)') "xsta",glo,".dat"
!         else
!         WRITE(filename,'(a,i3.3,a)') "xsta",glo,".dat"
!         endif
! 
! 
!         write(temp, 111) glo
!         CALL SYSTEM(temp)
! 
! 
!         open (unit = 666, file = filename)
!         write(666,*) 0,mult
!         write(666,*) 0,nprimal
!         do i = 1,nprimal
!         write(666,*) i,x_k(i)
!         x_k(i)=0.0d0
!         enddo
!         close(666)
!         deallocate(x_k)
!         enddo
! 
! 111     format ('cp stres.dat stres',I0,'.dat | grep -v Warning:')
        return
        end  

        
        
!         program global
! 
       subroutine globalD2 (idrive,x,nprimal,nmat,nsta,ndual,nnz,du,
     &      obj1,con1,eqn1,p) 

        implicit none
        logical feasible
        integer mult,ground,approx,prob,solver,trust,idrive
        double precision filt,dl0,dlc,du,dlm
        integer exprm,nex,ney,nsta,nmat,bs,nprimal
        integer ndual,nnz,rs,glo,i,r,loops,consv
        integer warm, conti, scal, conv
        double precision  sorm, filt0, x0,p
        real t,tt,seconds,approxseconds
        character*12 :: filename
        character (len=100) :: temp

        double precision f_k,f_star,h_k,h_star,pr
        double precision,dimension(:),allocatable::x_k
        double precision,dimension(:),allocatable::x_star
 
        double precision x(nprimal),objective,violation 
        double precision obj1,solution(nprimal)
        double precision con1(ndual)
        logical eqn1(ndual)
        double precision nobj1(nprimal)
        double precision jac1(nnz)
        integer aptr1(ndual+1)
        integer acol1(nnz)
        
        !!!
        include 'setsea.h'
        !!!

        filt=filt0

        t=seconds()
        call init (mult,nex,ney,ground,prob,nsta,nmat,
     &  bs,nprimal,ndual,nnz,approx,filt,solver,
     &  dl0,dlc,du,trust,consv,rs)


!        do glo = 1,exprm

        filt=filt0

        allocate(x_k(nprimal))
        
        x_k = x
        

        call seasand1(mult,ground,approx,prob,filt,solver,
     &  trust,consv,rs,dl0,dlc,objective,violation,feasible,solution,
     &  loops,glo,warm,approxseconds,conti,dlm,conv,x0,
     &  obj1,con1,nobj1,jac1,aptr1,acol1,eqn1,
     &  ndual,nprimal,nnz,p)     
     
     
     
     
     
! !       write history
!         filename = 'global.dat'
!         if(glo.eq.1) then    
!         open(unit=9999, file=filename, action="write",
!      &  status="replace")
!         else
!         open(unit=9999, file=filename, action="write",
!      &  access="append")
!         endif
!         tt=seconds()
!         write(9999,*)glo,loops,f_k,h_k,feasible,tt-t
!         close(9999)
! 
!         if(glo.lt.10) then
!         WRITE(filename,'(a,i1.1,a)') "xsta",glo,".dat"  
!         elseif(glo.lt.100) then
!         WRITE(filename,'(a,i2.2,a)') "xsta",glo,".dat"
!         else
!         WRITE(filename,'(a,i3.3,a)') "xsta",glo,".dat"
!         endif
! 
! 
!         write(temp, 111) glo
!         CALL SYSTEM(temp)
! 
! 
!         open (unit = 666, file = filename)
!         write(666,*) 0,mult
!         write(666,*) 0,nprimal
!         do i = 1,nprimal
!         write(666,*) i,x_k(i)
!         x_k(i)=0.0d0
!         enddo
!         close(666)
          deallocate(x_k)
!         enddo
! 
! 111     format ('cp stres.dat stres',I0,'.dat | grep -v Warning:')
        return
        end 

        
!         program global
! 
        subroutine globalD3 (idrive,x,nprimal,nmat,nsta,ndual,nnz,du,
     &      nobj1,jac1,aptr1,acol1,eqn1,p) 

        implicit none
        logical feasible
        integer mult,ground,approx,prob,solver,trust,idrive
        double precision filt,dl0,dlc,du,dlm
        integer exprm,nex,ney,nsta,nmat,bs,nprimal
        integer ndual,nnz,rs,glo,i,r,loops,consv
        integer warm, conti, scal, conv
        double precision  sorm, filt0, x0,p
        real t,tt,seconds,approxseconds
        character*12 :: filename
        character (len=100) :: temp

        double precision f_k,f_star,h_k,h_star,pr
        double precision,dimension(:),allocatable::x_k
        double precision,dimension(:),allocatable::x_star
 
        double precision x(nprimal),objective,violation
        double precision obj1,solution(nprimal)
        double precision con1(ndual) 
        double precision nobj1(nprimal) 
        double precision jac1(nnz) 
        integer aptr1(ndual+1) 
        integer acol1(nnz) 
        logical eqn1(ndual)
        
        !!!
        include 'setsea.h'
        !!!

        filt=filt0

        t=seconds()
        call init (mult,nex,ney,ground,prob,nsta,nmat,
     &  bs,nprimal,ndual,nnz,approx,filt,solver,
     &  dl0,dlc,du,trust,consv,rs)


!        do glo = 1,exprm

        filt=filt0

        allocate(x_k(nprimal))
        
        x_k = x
        
        call seasand1(mult,ground,approx,prob,filt,solver,
     &  trust,consv,rs,dl0,dlc,objective,violation,feasible,solution,
     &  loops,glo,warm,approxseconds,conti,dlm,conv,x0,
     &  obj1,con1,nobj1,jac1,aptr1,acol1,eqn1,
     &  ndual,nprimal,nnz,p)     
             
! !       write history
!         filename = 'global.dat'
!         if(glo.eq.1) then    
!         open(unit=9999, file=filename, action="write",
!      &  status="replace")
!         else
!         open(unit=9999, file=filename, action="write",
!      &  access="append")
!         endif
!         tt=seconds()
!         write(9999,*)glo,loops,f_k,h_k,feasible,tt-t
!         close(9999)
! 
!         if(glo.lt.10) then
!         WRITE(filename,'(a,i1.1,a)') "xsta",glo,".dat"  
!         elseif(glo.lt.100) then
!         WRITE(filename,'(a,i2.2,a)') "xsta",glo,".dat"
!         else
!         WRITE(filename,'(a,i3.3,a)') "xsta",glo,".dat"
!         endif
! 
! 
!         write(temp, 111) glo
!         CALL SYSTEM(temp)
! 
! 
!         open (unit = 666, file = filename)
!         write(666,*) 0,mult
!         write(666,*) 0,nprimal
!         do i = 1,nprimal
!         write(666,*) i,x_k(i)
!         x_k(i)=0.0d0
!         enddo
!         close(666)
          deallocate(x_k)
!         enddo
! 
! 111     format ('cp stres.dat stres',I0,'.dat | grep -v Warning:')
        return
        end 

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
! 
!         subroutine globalD1 (idrive,nprimal,nmat,nsta,ndual,nnz,
!      &      du,ground,prob) 
! 
!         implicit none
!         logical feasible
!         integer mult,ground,approx,prob,solver,trust
!         double precision filt,dl0,dlc,du
!         integer exprm,nex,ney,nsta,nmat,bs,nprimal
!         integer ndual,nnz,rs,glo,i,r,loops,consv
!         integer idrive
!         double precision conv, sorm, filt0
!         !external conv
!         real t,tt,seconds
!         character*12 :: filename
! 
!         double precision f_k,f_star,h_k,h_star,pr
! !         double precision,dimension(:),allocatable::x_k
! !         double precision,dimension(:),allocatable::x_star
!         
!         
!         
!     
!         !!!
!         include 'setsea.h'
!         !!!
! 
!         
!         
! !        t=seconds()
! 
!         call init (mult,nex,ney,ground,prob,nsta,nmat,
!      &  bs,nprimal,ndual,nnz,approx,filt,solver,
!      &  dl0,dlc,du,trust,consv,rs)
! 
! 
! !         do glo = 1,exprm
! ! 
! !         allocate(x_k(nprimal))
! !         filt=filt0
! !         call seasand(mult,ground,approx,prob,filt,solver,
! !      &  trust,consv,rs,dl0,dlc,f_k,h_k,feasible,x_k,
! !      &  loops,glo)
! ! 
! ! !       write history
! !         filename = 'global.dat'
! !         if(glo.eq.1) then    
! !         open(unit=9999, file=filename, action="write",
! !      &  status="replace")
! !         else
! !         open(unit=9999, file=filename, action="write",
! !      &  access="append")
! !         endif
! !         tt=seconds()
! !         write(9999,*)glo,loops,f_k,h_k,feasible,tt-t
! !         close(9999)
! ! 
! !         WRITE(filename,'(a,i4.4,a)') "xsta",glo,".dat"
! !         open (unit = 666, file = filename)
! !         write(666,*) 0,mult
! !         write(666,*) 0,nprimal
! !         do i = 1,nprimal
! !         write(666,*) i,x_k(i)
! !         x_k(i)=0.0d0
! !         enddo
! !         close(666)
! !         deallocate(x_k)
! !         enddo
!         
!         return 
!         end 
        
!         subroutine globalD2 (idrive,x,nprimal,nmat,nsta,ndual,nnz,du,
!      &      obj1,con1,eqn1,p) 
! 
!         implicit none
!         logical feasible
!         integer mult,ground,approx,prob,solver,trust
!         double precision filt,dl0,dlc,du
!         integer exprm,nex,ney,nsta,nmat,bs,nprimal
!         integer ndual,nnz,rs,glo,i,r,loops,consv
!         integer idrive
!         double precision conv, sorm, filt0,p
!         external conv
!         real t,tt,seconds
!         character*12 :: filename
! 
!         double precision f_k,f_star,h_k,h_star,pr
!         double precision,dimension(:),allocatable::x_k
!    
!         double precision x(nprimal) 
!         double precision obj1 
!         double precision con1(ndual)
!         logical eqn1(ndual)
!         double precision nobj1(nprimal)
!         double precision jac1(nnz)
!         integer aptr1(ndual+1)
!         integer acol1(nnz)
!     
!         !!!
!         include 'setsea.h'
!         !!!        
!         
! !       t=seconds()
! 
!         call init (mult,nex,ney,ground,prob,nsta,nmat,
!      &  bs,nprimal,ndual,nnz,approx,filt,solver,
!      &  dl0,dlc,du,trust,consv,rs)
! 
! ! 
! ! !         do glo = 1,exprm
! ! 
!           allocate(x_k(nprimal))
! !         
!           x_k = x
! !         
! !         filt=filt0
! 
!         call seasand1(mult,ground,approx,prob,filt,solver,
!      &  trust,consv,rs,dl0,dlc,f_k,h_k,feasible,x_k,
!      &  loops,glo,obj1,con1,nobj1,jac1,aptr1,acol1,eqn1,
!      &  ndual,nprimal,nnz,p)
!      
!    
! ! 
! ! ! !       write history
! ! !         filename = 'global.dat'
! ! !         if(glo.eq.1) then    
! ! !         open(unit=9999, file=filename, action="write",
! ! !      &  status="replace")
! ! !         else
! ! !         open(unit=9999, file=filename, action="write",
! ! !      &  access="append")
! ! !         endif
! ! !         tt=seconds()
! ! !         write(9999,*)glo,loops,f_k,h_k,feasible,tt-t
! ! !         close(9999)
! ! ! 
! ! !         WRITE(filename,'(a,i4.4,a)') "xsta",glo,".dat"
! ! !         open (unit = 666, file = filename)
! ! !         write(666,*) 0,mult
! ! !         write(666,*) 0,nprimal
! ! !         do i = 1,nprimal
! ! !         write(666,*) i,x_k(i)
! ! !         x_k(i)=0.0d0
! ! !         enddo
! ! !         close(666)
! !         deallocate(x_k)
! ! !         enddo
!         
!         
!         return 
!         end 
!         
        
!         subroutine globalD3 (idrive,x,nprimal,nmat,nsta,ndual,nnz,du,
!      &      nobj1,jac1,aptr1,acol1,eqn1,p) 
! 
!         implicit none
!         logical feasible
!         integer mult,ground,approx,prob,solver,trust
!         double precision filt,dl0,dlc,du
!         integer exprm,nex,ney,nsta,nmat,bs,nprimal
!         integer ndual,nnz,rs,glo,i,r,loops,consv
!         integer idrive
!         double precision conv, sorm, filt0,p
!         external conv
!         real t,tt,seconds
!         character*12 :: filename
! 
!         double precision f_k,f_star,h_k,h_star,pr
!         double precision,dimension(:),allocatable::x_k
!    
!         double precision x(nprimal) 
!         double precision obj1 
!         double precision con1(ndual) 
!         double precision nobj1(nprimal) 
!         double precision jac1(nnz) 
!         integer aptr1(ndual+1) 
!         integer acol1(nnz) 
!         logical eqn1(ndual)
!         
!         
!     
!         !!!
!         include 'setsea.h'
!         !!!
! 
!         
!         
! !        t=seconds()
! 
!         call init (mult,nex,ney,ground,prob,nsta,nmat,
!      &  bs,nprimal,ndual,nnz,approx,filt,solver,
!      &  dl0,dlc,du,trust,consv,rs)
! 
! 
! !         do glo = 1,exprm
! 
!         allocate(x_k(nprimal))
!         
!         x_k = x
!         
!         filt=filt0
!         call seasand1(mult,ground,approx,prob,filt,solver,
!      &  trust,consv,rs,dl0,dlc,f_k,h_k,feasible,x_k,
!      &  loops,glo,obj1,con1,nobj1,jac1,aptr1,acol1,eqn1,
!      &  ndual,nprimal,nnz,p)
! 
! ! !       write history
! !         filename = 'global.dat'
! !         if(glo.eq.1) then    
! !         open(unit=9999, file=filename, action="write",
! !      &  status="replace")
! !         else
! !         open(unit=9999, file=filename, action="write",
! !      &  access="append")
! !         endif
! !         tt=seconds()
! !         write(9999,*)glo,loops,f_k,h_k,feasible,tt-t
! !         close(9999)
! ! 
! !         WRITE(filename,'(a,i4.4,a)') "xsta",glo,".dat"
! !         open (unit = 666, file = filename)
! !         write(666,*) 0,mult
! !         write(666,*) 0,nprimal
! !         do i = 1,nprimal
! !         write(666,*) i,x_k(i)
! !         x_k(i)=0.0d0
! !         enddo
! !         close(666)
!         deallocate(x_k)
! !         enddo
!         
!         return 
!         end 
!         
        
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
        
        subroutine drivegurobi(outerloop,nprimal,ndual,nnz,
     &  xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
     &  xl,xu,dl,du,nsta,dx,feasible)

        implicit            none
        integer             i,j,k
        logical             eqn(ndual)
        integer             nprimal,ndual,nnz,aptr(ndual+1)
        double precision    xprimal(nprimal),xdual(ndual)
        double precision    obj,nobj(nprimal),con(ndual)
        double precision    jac(nnz),hes(nprimal),dnorm
        double precision    xl(nprimal),xu(nprimal),biglam
        double precision    qpbl(nprimal),qpbu(nprimal)
        double precision    delta(nprimal),lambda(ndual)
        double precision    dx(nprimal),dl,du
        integer             acol(nnz),outerloop,nsta
        integer             status
        logical             feasible
        integer             bnprimal,bndual,bnnz

!       double precision,dimension(:),allocatable::bnobj

!       construct the primal variable bounds
        do i=1,nprimal
        delta(i)=0.0d0
        dx(i)=0.0d0
!       if(i.gt.nsta) then
        qpbl(i)=max(-dl*(xu(i)-xl(i)),(xl(i)-xprimal(i)))
        qpbu(i)=min(dl*(xu(i)-xl(i)),(xu(i)-xprimal(i)))
!       else
!       qpbl(i)=-du
!       qpbu(i)=du
        if(xl(i).eq.xu(i)) then
        qpbl(i)=0.0d0
        qpbu(i)=0.0d0
        endif
!       endif
        enddo
        do j=1,ndual
        lambda(j) = 0.0d0
        enddo

!       call gurobi_solve_sqps(status,nprimal,ndual,nnz,eqn,
!    &  delta,lambda,nobj,con,jac,hes,aptr,acol,qpbl,qpbu)

        feasible=.true.
        if (status.eq.2) then
        feasible=.false.
        endif

!       update the primal variables
        do i=1,nprimal
        dx(i) = delta(i)
!       xprimal(i)=min(xu(i),max(xl(i),xprimal(i)+delta(i)))
        enddo

        biglam = 1d6
!       update the lagrange multipliers
        do j=1,ndual
        xdual(j)= max(-biglam,min(-lambda(j),biglam))
        enddo

        return
        end subroutine
        
        subroutine init (mult,nex,ney,ground,prob,nsta,
     &  nmat,bs,nprimal,ndual,nnz,approx,filt,solver,
     &  dl0,dlc,du,trust,consv,rs)
     
     
     
        implicit none
        integer             mult,ground,prob,nsta,nmat
        integer             nprimal,ndual,nnz,approx
        integer             nex,ney,tmp1,tmp2,tmp3,tmp4,bs
        double precision    filt, dl0, dlc, du
        integer             solver,trust,consv,rs

        du=1d6

        nex=15*mult
        ney=5*mult
        
        if(ground.eq.3) then
        nex=25*mult
        ney=25*mult
        endif

        bs = (nex-1)*ney+(ney-1)*nex
        nsta=2*(2*nex+1)*(ney+1)+2*ney*(nex+1)
        nmat=nex*ney
        tmp1=ney*nex*256-36*(ney*(nex-1)+nex*(ney-1))+4*(ney-1)*(nex-1)
        tmp2=24+2*2*(nex-2)+2*2*(ney-2)+8*(nex-1)+8*(ney-1)
        tmp3=8*(nex-1)*(ney-1)+4*(nex-1)*ney+4*(ney-1)*nex

        if(prob.eq.1) then
        nprimal=nsta+nmat
        ndual=nsta+nmat
        tmp4=9*nex*ney
        elseif(prob.eq.2) then
        nprimal=nsta+nmat
        ndual=nsta+1
        tmp4=nmat
        elseif(prob.eq.3) then
        nprimal=nsta+nmat
        ndual=nsta+1+2*bs
        tmp4=nmat+2*2*bs
        elseif(prob.eq.4) then
        nprimal=nsta+nmat
        ndual=nsta+nmat+2*bs
        tmp4=9*nex*ney+2*2*bs
        elseif(prob.eq.99) then
        nprimal=nsta+nmat
        ndual=nsta+1
        tmp4=nmat
        elseif(prob.eq.100) then
        nprimal=nsta+nmat
        ndual=nsta+nmat
        tmp4=9*nex*ney
        endif
        nnz=tmp1+tmp2+tmp3+tmp4
        write(*,*) ndual, nmat, nsta
        end subroutine init
        subroutine maproutine(ex,ey,nex,ney,map)
        implicit none
        integer         ex,ey,nex,ney
        integer         map(16)

!       Local to global mapping
        !-->u1
        map(1) = 2*(ey-1)*(2*nex+1) + 2*(ey-1)*(nex+1) + 2*2*(ex-1) +1
        !-->v1
        map(2) = map(1) + 1
        !-->u2
        map(3) = 2*(ey-1)*(2*nex+1) + 2*(ey-1)*(nex+1) + 2*2*(ex-1) +5
        !-->v2
        map(4) = map(3) + 1
        !-->u3
        map(5) = 2*ey*(2*nex+1) + 2*ey*(nex+1) + 2*2*(ex-1) +5
        !-->v3
        map(6) = map(5) + 1
        !-->u4
        map(7) = 2*ey*(2*nex+1) + 2*ey*(nex+1) + 2*2*(ex-1) +1
        !-->v4
        map(8) = map(7) + 1
        !-->u5
        map(9) = 2*(ey-1)*(2*nex+1) + 2*(ey-1)*(nex+1) + 2*2*(ex-1) +3
        !-->v5
        map(10) = map(9) + 1
        !-->u6
        map(11) = 2*ey*(2*nex+1) + 2*(ey-1)*(nex+1) + 2*ex +1
        !-->v6
        map(12) = map(11) + 1
        !-->u7
        map(13) = 2*ey*(2*nex+1) + 2*ey*(nex+1) + 2*2*(ex-1) +3
        !-->v7
        map(14) = map(13) + 1
        !-->u8
        map(15) = 2*ey*(2*nex+1) + 2*(ey-1)*(nex+1) + 2*ex -1
        !-->v8
        map(16) = map(15) + 1

        return
        end subroutine
!         subroutine init_random_seed()
! 
!         INTEGER :: i, n, clock
!         INTEGER, DIMENSION(:), ALLOCATABLE :: seed
! 
!         CALL RANDOM_SEED(size = n)
!         ALLOCATE(seed(n))
! 
!         CALL SYSTEM_CLOCK(COUNT=clock)
! 
!         seed = clock + 37 * (/ (i - 1, i = 1, n) /)
!         CALL RANDOM_SEED(PUT = seed)
! 
!         DEALLOCATE(seed)
!         end
! seasand:
! seasand: mon august 10 17:12 sast 2015, dirk munro, bellville
! seasand: main driver
! seasand:

        subroutine seasand1(mult,ground,approx,prob,filt,solver,
     &  trust,consv,rs,dl0,dlc,objective,violation,feasible,solution,
     &  loops,glo,warm,approxseconds,conti,dlm,conv,x0,
     &  obj1,con1,nobj1,jac1,aptr1,acol1,eqn1,
     &  ndual,nprimal,nnz,p)
     

        implicit none
        integer i,j,k,bs,nex,ney,glo
        integer nprimal,ndual,nnz,nmat
        integer nprimalc,ndualc,nnzc
        integer mult,ground,prob,nsta
        integer outerloop,innerloop,ff
        integer actst,actstv,actsl,conti
        integer loops,approx,rem,rs,warm
        integer n_lo,n_up,solver,trust,conv
        integer fstep,hstep,consv,xpc,io
        logical feasible,next
        double precision filterlist(1000000,2), dmax
        double precision obj,objh,dnorm,violh,xp,cdx
        double precision objt,violt,viol,p,filt,dlc,x0
        double precision dl,dl0,du,temp,hmult,curv0,dlm
        double precision objective,violation,solution(*)
        double precision xmax,xmin,mstres,filt0,curvp
        real cputimes, cputimef, time(3),seconds
        real approxs,approxf,approxseconds
        character*12 :: filename
        logical,dimension(:),allocatable::eqn
        double precision,dimension(:),allocatable::xprimal
        double precision,dimension(:),allocatable::xprimalc
        double precision,dimension(:),allocatable::xprimalh
        double precision,dimension(:),allocatable::xprimalt
        double precision,dimension(:),allocatable::xdual
        double precision,dimension(:),allocatable::xdualc
        double precision,dimension(:),allocatable::xdualh
        double precision,dimension(:),allocatable::xdualt
        double precision,dimension(:),allocatable::nobj
        double precision,dimension(:),allocatable::nobjh
        double precision,dimension(:),allocatable::nobjt
        double precision,dimension(:),allocatable::con
        double precision,dimension(:),allocatable::conh
        double precision,dimension(:),allocatable::cont
        double precision,dimension(:),allocatable::jac
        double precision,dimension(:),allocatable::jach
        double precision,dimension(:),allocatable::jact
        double precision,dimension(:),allocatable::hes
        double precision,dimension(:),allocatable::hesh
        double precision,dimension(:),allocatable::hesk
        double precision,dimension(:),allocatable::hest
        double precision,dimension(:),allocatable::xl
        double precision,dimension(:),allocatable::xu
        double precision,dimension(:),allocatable::stres
        double precision,dimension(:),allocatable::dx
        double precision,dimension(:),allocatable::turvj
        integer,dimension(:),allocatable::aptr
        integer,dimension(:),allocatable::aptrc
        integer,dimension(:),allocatable::acol
        integer,dimension(:),allocatable::bc
        integer,dimension(:),allocatable::void
        
        double precision obj1 
        double precision con1(ndual)
        double precision nobj1(nprimal)
        double precision jac1(nnz)
        integer aptr1(ndual+1)
        integer acol1(nnz)
        logical eqn1(ndual)

        
        
        
    
        ff=0
        rem=0
        time(1) = 0
        time(2) = 0
        time(3) = seconds()
        cputimes = 0
        cputimef = 0
        approxf=0
        approxs=0

!       initialize the groundstructure / problem
        call init (mult,nex,ney,ground,prob,nsta,nmat,
     &  bs,nprimal,ndual,nnz,approx,filt,solver,
     &  dl0,dlc,du,trust,consv,rs)

!       initialize storage vectors
        allocate (xprimal(nprimal))
        allocate (xprimalc(nprimal))
        allocate (xprimalh(nprimal))
        allocate (stres(nmat))
        allocate (xprimalt(nprimal))
        allocate (eqn(ndual))
        allocate (xdual(ndual))
        allocate (xdualc(ndual))
        allocate (xdualh(ndual))
        allocate (xdualt(ndual))
        allocate (nobj(nprimal))
        allocate (nobjh(nprimal))
        allocate (nobjt(nprimal))
        allocate (con(ndual))
        allocate (conh(ndual))
        allocate (cont(ndual))
        allocate (jac(nnz))
        allocate (jach(nnz))
        allocate (jact(nnz))
        allocate (hes(nprimal))
        allocate (hesh(nprimal))
        allocate (hesk(nprimal))
        allocate (turvj(ndual))
        allocate (hest(nprimal))
        allocate (aptr(ndual+1))
        allocate (aptrc(ndual+1))
        allocate (acol(nnz))
        allocate (xl(nprimal))
        allocate (xu(nprimal))
        allocate (dx(nprimal))
        allocate (bc(nprimal))
        allocate (void(nprimal))

        !warmstart
        k=0
        xpc=0
        if(warm.eq.1) then
        open(unit=99999,file='xwarm.dat')
        do
        read(99999,*,iostat=io) xpc,xp
        if(io/=0) exit
        xprimal(xpc) = xp
        enddo
        close(99999)
        if(xpc.ne.nprimal) then
        write(*,*) 'ERROR', xpc, nprimal
        stop
        endif
        endif

        !call init_random_seed
        feasible=.true.

!        p=3.0d0
        obj=0.0d0
        objh=0.0d0
        objt=0.0d0
        do i=1,nprimal
        dx(i) = 0.0d0
        nobj(i)=0.0d0
        nobjh(i)=0.0d0
        nobjt(i)=0.0d0
        hes(i)=0.0d0
        hesh(i)=0.0d0
        hest(i)=0.0d0
        hesk(i)=0.0d0
        if(i.le.nsta) then
        if(warm.eq.0) then
        xprimal(i)= 0.0d0
        endif
        xl(i)=-du
        xu(i)=du
        endif
        if(i.gt.nsta) then
        if(rs.eq.1) then
        call random_number (temp)
        xprimal(i)=temp*0.01d0+0.99d0
        elseif (rs.eq.0) then
        xprimal(i)=x0
        endif
        xprimalh(i)=xprimal(i)
        xprimalt(i)=xprimal(i)
        xl(i)=0d0
        xu(i)=1d0
        if(void(i).eq.1) then
        xu(i)=0.0d0
        endif
        endif
        enddo
        do i = 1,nmat
        stres(i) = 0.0d0
        enddo

        do i=1,ndual
        xdual(i)=0.0d0
        con(i)=0.0d0
        if(i.le.nsta) then
        eqn(i)=.true.
        endif
        if(i.gt.nsta) then
        eqn(i)=.false.
        endif
        enddo
        do i=1,nnz
        jac(i)=0.0d0
        enddo
        dnorm = 0.0d0

        rem=1                     ! AG2017
        dl=dl0
        fstep=0
        hstep=0
        hmult=1.0d0
        !do i =1,1000000          ! AG2017
        filterlist(1,1) = 1.0d16
        filterlist(1,2) = 1.0d1   ! AG2017
        !enddo

!       call simulation at point x_0
        cputimes=seconds()
        if(prob.eq.1) then
        call simu_1(mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
     &           xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
     &           xl,xu,p,viol,stres)
        elseif(prob.eq.2) then
        call simu_2(mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
     &           xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
     &           xl,xu,p,viol)
        elseif(prob.eq.3) then
        call simu_3(mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
     &          xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
     &          xl,xu,p,bs,filt,viol)
        elseif(prob.eq.4) then
        call simu_4(mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
     &          xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
     &          xl,xu,p,bs,filt,viol,stres,xprimalh,conh,
     &          turvj)
        elseif(prob.eq.5) then
!         call simu_5(mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
!      &          xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
!      &          xl,xu,p,bs,filt,viol,stres)
        endif 
        cputimef=seconds()
        time(1) = time(1) + cputimef - cputimes

                obj1   =  obj  
        con1   =  con  
        nobj1  =  nobj  
        jac1   =  jac  
        aptr1  =  aptr  
        acol1  =  acol  
        eqn1   =  eqn  

        return
        
        
        loops = 0
!       outer loop
        do outerloop = 1,999

        approxs=seconds()
!       curvature approximation and convexity
        cdx=0.0d0
        curv0=0.0d0
        if(approx.eq.2.or.approx.eq.3) then
        if(outerloop.gt.1) then
        curv0=objh-obj
        do i=1,nprimal
        curv0=curv0-nobj(i)*(xprimalh(i)-xprimal(i))
        cdx=cdx+(xprimalh(i)-xprimal(i))**2.0d0
        enddo
        cdx=max(cdx,1d-6)
        curv0=2.0d0*curv0/cdx
        do j=1,nsta+nmat
        turvj(j)=2.0d0*turvj(j)/cdx
        curv0=curv0+turvj(j)
        enddo
        else
        curv0=1d0
        endif
        endif

        curvp=0.0d0
        if(approx.eq.0) then
        do i = 1,nprimal
        hes(i) = max(1d-6,hes(i))
        curvp=max(curvp,hes(i))
        enddo
        elseif(approx.eq.1) then
        do i = 1,nprimal
        hes(i) = max(abs(hes(i)),1d-6)
        curvp=max(curvp,hes(i))
        enddo
        elseif(approx.eq.2) then
        do i = 1,nprimal
        hes(i)=max(1d-6,curv0)
        curvp=max(curvp,hes(i))
        enddo
        elseif(approx.eq.3) then
        do i = 1,nprimal
        hes(i)=max(1d-6,abs(curv0))
        curvp=max(curvp,hes(i))
        enddo
        endif

        
        curvp=0.0d0
        do i = 1,nprimal
        hes(i)=min(1d4,hes(i))
        curvp=max(curvp,hes(i))
        enddo
        
        
        approxf=seconds()
        approxseconds=approxseconds+approxf-approxs
!       set psub_h
        call storesub(nprimal,ndual,nnz,xprimal,xdual,
     &  obj,nobj,hes,jac,con,viol,violh,xprimalh,xdualh,objh,
     &  nobjh,hesh,conh,jach)

        fstep=0
        hstep=0
        do innerloop=1,1

!       solve psub
        feasible = .false.
        cputimes=seconds()
        loops=loops+1

!         if(solver.eq.1) then
!         call drivecplex(outerloop,nprimal,ndual,nnz, 
!      &  xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
!      &  xl,xu,dl,du,nsta,dx,feasible)
!         elseif(solver.eq.2) then
!         call drivegurobi(outerloop,nprimal,ndual,nnz,
!      &  xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
!      &  xl,xu,dl,du,nsta,dx,feasible)
!         endif
        cputimef=seconds()
        time(2) = time(2) + cputimef - cputimes

        if(feasible) then
          do i = 1,nprimal
            xprimal(i)=min(xu(i),max(xl(i),xprimal(i)+dx(i)))
          enddo
        else
          rem=rem+1
          if(rem.gt.1000000) then
            stop 'STOP: The filter is full'
          endif
          filterlist(rem,1) = objh
          filterlist(rem,2) = violh
!         simple backtracking
          do i = 1,nprimal
            dx(i) = xprimalt(i) - xprimalh(i)
          enddo
          call fetchsub(nprimal,ndual,nnz,xprimal,xdual,
     &    obj,nobj,hes,jac,con,viol,violt,xprimalt,xdualt,
     &    objt,nobjt,hest,cont,jact)
          dl=1.0d0
!         rem=0
          exit
        endif

!       call simulation at current point x_k
        cputimes=seconds()
        if(prob.eq.1) then
        call simu_1(mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
     &           xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
     &           xl,xu,p,viol,stres)
        elseif(prob.eq.2) then
        call simu_2(mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
     &           xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
     &           xl,xu,p,viol)
        elseif(prob.eq.3) then
        call simu_3(mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
     &         xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
     &         xl,xu,p,bs,filt,viol)
        elseif(prob.eq.4) then
        call simu_4(mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
     &          xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
     &          xl,xu,p,bs,filt,viol,stres,xprimalh,conh,
     &          turvj)
        elseif(prob.eq.5) then
!         call simu_5(mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
!      &          xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
!      &          xl,xu,p,bs,filt,viol,stres)
        endif
        cputimef=seconds()
        time(1) = time(1) + cputimef - cputimes

        
        !filter
!       if(feasible) then
        if(trust.eq.1) then
        next=.false.
        call filter(nprimal,obj,viol,objh,violh,nobjh,hesh,
     &  dx,dl,dl0,innerloop,rem,filterlist,next,hmult,
     &  fstep,hstep,consv,dlm)
        if(next) then
        exit
        else
        call fetchsub(nprimal,ndual,nnz,xprimal,xdual,
     &  obj,nobj,hes,jac,con,viol,violh,xprimalh,xdualh,objh,
     &  nobjh,hesh,conh,jach)
        endif
        else
        exit
        endif
!       endif

        enddo

        if(feasible) then
        call storesub(nprimal,ndual,nnz,xprimalh,xdualh,
     &  objh,nobjh,hesh,jach,conh,violh,violt,xprimalt,
     &  xdualt,objt,nobjt,hest,cont,jact)
        endif

        dnorm = 0.0d0
        dmax = 0.0d0
!       active bounds
        n_lo = 0
        n_up = 0
        do i = 1,nprimal
        dmax = max(dmax, abs(dx(i)))
        dnorm=dnorm+dx(i)**2.0d0
        if(i.gt.nsta) then
        if (xprimal(i).gt.1.0d0-5d-2) then
        n_up = n_up + 1
        endif
        endif
        if(i.gt.nsta) then
        if (xprimal(i).lt.1d-3) then
        n_lo = n_lo + 1
        endif
        endif
        enddo
        dnorm=dnorm**(1.0d0/2.0d0)
        if (conv.eq.1) then
        dnorm=dmax
        endif
        
        if(glo.lt.10) then
        WRITE(filename,'(a,i1.1,a)') "hist",glo,".dat"  
        elseif(glo.lt.100) then
        WRITE(filename,'(a,i2.2,a)') "hist",glo,".dat"
        else
        WRITE(filename,'(a,i3.3,a)') "hist",glo,".dat"
        endif
        

        mstres=0.0d0
        do i =1,nmat
        if(xprimal(nsta+i).gt.1d-6) then
        mstres=max(mstres,stres(i))
        endif
        enddo

        actst=0
        actstv=0
        actsl=0
        do i = 1,ndual
!       if(eqn(i)) then
!       else
        if(xdual(i).ne.0d0) then
        if(i.gt.nmat+nsta) then
        actsl=actsl+1
        else
        actst=actst+1
        if(xprimal(i).gt.0d0) then
        actstv=actstv+1
        endif
        endif
!       endif
        endif
        enddo


        if(outerloop.eq.1) then    
        open(unit=9999, file=filename, action="write",
     &  status="replace")
        else
        open(unit=9999, file=filename, action="write",
     &  access="append")
        endif
        write(9999,1100)outerloop,innerloop,dl,obj,viol,
     &  n_lo,n_up,dnorm, feasible, curvp, fstep, hstep,mstres,
     &  actst,actstv,actsl
        write(*,1100)outerloop,innerloop,dl,obj,viol,
     &  n_lo,n_up,dnorm, feasible, curvp, fstep, hstep,mstres,
     &  actst,actstv,actsl
        close(9999)
        if(ground.eq.3) then
!       call lopoprint_k(mult,xprimal,nprimal,nex,ney)
        else
!       call stresprint_k(mult,stres,nmat)
        call topoprint_k(mult,xprimal,nprimal,nex,ney)
        endif

        write(*,*) dlc, dnorm

!       termination?
        if(feasible) then
         if(conti.eq.0) then
          if(dnorm.lt.dlc) then
          exit
          endif
         endif
         if(conti.eq.1) then
          if(dnorm.lt.dlc.and.filt.gt.dble(nex)) then
          exit
          endif
          if(dnorm.lt.dlc.and.filt.le.dble(nex)) then
          if(ff.eq.0) then
          WRITE(filename,'(a)') "con1_bc.dat"
          open (unit = 7777, file = filename)
          write(7777,*) 0,mult
          write(7777,*) 0,ndual
          do i = 1,ndual
          write(7777,*) i,con(i),xdual(i)
          enddo
          close(7777)

          WRITE(filename,'(a)') "xsta1_bc.dat"
          open (unit = 777, file = filename)
          write(777,*) 0,mult
          write(777,*) 0,nprimal
          do i = 1,nprimal
          write(777,*) i,xprimal(i)
          enddo
          close(777)

          CALL SYSTEM('cp stres.dat stres1_bc.dat 
     &    | grep -v Warning:')
          CALL SYSTEM('cp hist1.dat hist1_bc.dat 
     &    | grep -v Warning:')
          CALL SYSTEM('cp topo_k.png topo_bc.png
     &    | grep -v Warning:')
          ff=1
          endif
          filt=filt*1.1d0
          rem=0
          endif
         endif
        endif

         if(conti.eq.0) then
          if(dnorm.eq.0.0d0.or.obj.eq.0.0d0) then
!         rem=0
          exit
          endif
         endif

        if(conti.eq.2) then
         if(dnorm.eq.0.0d0.or.obj.eq.0.0d0) then
          if(ff.eq.0) then

          WRITE(filename,'(a)') "con1_bc.dat"
          open (unit = 7777, file = filename)
          write(7777,*) 0,mult
          write(7777,*) 0,ndual
          do i = 1,ndual
          write(7777,*) i,con(i),xdual(i)
          enddo
          close(7777)

          WRITE(filename,'(a)') "xsta1_bc.dat"
          open (unit = 777, file = filename)
          write(777,*) 0,mult
          write(777,*) 0,nprimal
          do i = 1,nprimal
          write(777,*) i,xprimal(i)
          enddo
          close(777)

          CALL SYSTEM('cp stres.dat stres1_bc.dat 
     &    | grep -v Warning:')
          CALL SYSTEM('cp hist1.dat hist1_bc.dat 
     &    | grep -v Warning:')
          CALL SYSTEM('cp topo_k.png topo_bc.png
     &    | grep -v Warning:')
          ff=1
          endif
          if(filt.gt.dble(nex)) then
          exit
          endif
          if(outerloop.le.10) then
          exit
          endif
          filt=filt*1.1d0
          rem=0
         endif
        endif

        enddo

        cputimef=seconds()
        time(3) = cputimef- time(3)

        if (ground.eq.3) then
        objective=obj/(20**2.0d0*mult**2.0d0)
        else
        objective=obj/nex/ney
        endif
        violation=viol
        do i=1,nprimal
        solution(i)=xprimal(i)
        enddo
        if(glo.lt.10) then
        WRITE(filename,'(a,i1.1,a)') "con",glo,".dat"  
        elseif(glo.lt.100) then
        WRITE(filename,'(a,i2.2,a)') "con",glo,".dat"
        else
        WRITE(filename,'(a,i3.3,a)') "con",glo,".dat"
        endif
        open (unit = 6666, file = filename)
        write(6666,*) 0,mult
        write(6666,*) 0,ndual
        do i = 1,ndual
        write(6666,*) i,con(i),xdual(i)
        enddo
        close(6666)

!       kill storage vectors
        deallocate (xprimal)
        deallocate (xprimalc)
        deallocate (xprimalh)
        deallocate (xprimalt)
        deallocate (eqn)
        deallocate (xdual)
        deallocate (xdualc)
        deallocate (xdualh)
        deallocate (xdualt)
        deallocate (nobj)
        deallocate (nobjh)
        deallocate (nobjt)
        deallocate (con)
        deallocate (conh)
        deallocate (cont)
        deallocate (jac)
        deallocate (jach)
        deallocate (jact)
        deallocate (hes)
        deallocate (hesh)
        deallocate (hest)
        deallocate (hesk)
        deallocate (turvj)
        deallocate (aptr)
        deallocate (aptrc)
        deallocate (acol)
        deallocate (stres)
        deallocate (xl)
        deallocate (xu)
        deallocate (dx)
        deallocate (bc)
        deallocate (void)
        return

1000    format ('  k  l',4x,' dl',9x,'f',15x,'c',15x,
     &  'n_lo',9x,'n_up',9x,'d',9x,'eps')
1100    format (i3,1x,i2,4x,E7.1,6x,ES10.3,7x,ES9.3,7x,
     &  i6,6x,i6,9x,ES9.3,9x,L1,8x,ES10.3,9x,i2,9x,i2,
     &   9x,ES9.3,3x,i12,3x,i12,3x,i12)    

        end subroutine

        subroutine storesub(nprimal,ndual,nnz,xprimal,xdual,
     &  obj,nobj,hes,jac,con,viol,viols,xprimals,xduals,objs,
     &  nobjs,hess,cons,jacs)
        integer nprimal,ndual,nnz
        double precision obj,objs,viol,viols
        double precision xprimal(nprimal),xdual(ndual)
        double precision xprimals(nprimal),xduals(ndual)
        double precision nobj(nprimal),nobjs(nprimal)
        double precision hes(nprimal),hess(nprimal)
        double precision jac(nnz),jacs(nnz),s
        double precision con(ndual),cons(ndual)
        objs=obj
        viols=viol
        do i =1,nprimal
        xprimals(i)=xprimal(i)
        nobjs(i)=nobj(i)
        hess(i)=hes(i)
        enddo
        do i =1,ndual
        xduals(i)=xdual(i)
        cons(i)=con(i)
        enddo
        do i =1,nnz
        jacs(i)=jac(i)
        enddo
        end subroutine

        subroutine fetchsub(nprimal,ndual,nnz,xprimal,xdual,
     &  obj,nobj,hes,jac,con,viol,viols,xprimals,xduals,objs,
     &  nobjs,hess,cons,jacs)
        integer nprimal,ndual,nnz
        double precision obj,objs,viol,viols
        double precision xprimal(nprimal),xdual(ndual)
        double precision xprimals(nprimal),xduals(ndual)
        double precision nobj(nprimal),nobjs(nprimal)
        double precision hes(nprimal),hess(nprimal)
        double precision jac(nnz),jacs(nnz), s
        double precision con(ndual),cons(ndual)
        obj=objs
        viol=viols
        do i =1,nprimal
        xprimal(i)=xprimals(i)
        nobj(i)=nobjs(i)
        hes(i)=hess(i)
        enddo
        do i =1,ndual
        xdual(i)=xduals(i)
        con(i)=cons(i)
        enddo
        do i =1,nnz
        jac(i)=jacs(i)
        enddo
        end subroutine
        subroutine simu_1 (mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
     &               xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
     &               xl,xu,p,viol,stres) !stres,p,Es,scal)
!       call simu_1(mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
!    &           xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
!    &           xl,xu,p,viol)
        implicit none
        integer          i, j, k
        integer          mult,ground,prob,nsta,nmat,nprimal,ndual,nnz
        integer          nex,ney,ep,ex,ey,map(16),fxy(2,2),nsf
        double precision xprimal(nprimal),nobj(nprimal)
        double precision obj,con(ndual),jac(nnz),hes(nprimal),p,ke,sy
        double precision u5,u6,u7,u8,v5,v6,v7,v8,sig(3),E,nu,l,fvec(16)
        integer          connect(nsta,42),nodes(nsta),acol(nnz)
        integer          rowcounter(ndual),elmd(16),aptr(ndual+1)
        integer          contmp(nsta), bmap(6)
        double precision x_elm, xdual(ndual), xl(nprimal), xu(nprimal)
        double precision stres(nmat)
        double precision viol
        logical          eqn(ndual)
!
        nu=0.3d0
        l=3.0d0/15.0d0/mult
        sy=20.0d0
        do i =1,nmat
        if(xprimal(nsta+i)**p.le.1d-3) then
        xprimal(nsta+i)=0.0d0
        endif
        enddo

        elmd(1) = 1
        elmd(2) = 2
        elmd(3) = 9
        elmd(4) = 10
        elmd(5) = 3
        elmd(6) = 4
        elmd(7) = 15
        elmd(8) = 16
        elmd(9) = 11
        elmd(10) = 12
        elmd(11) = 7
        elmd(12) = 8
        elmd(13) = 13
        elmd(14) = 14
        elmd(15) = 5
        elmd(16) = 6

        ! safety
        obj=0.0d0
        do i=1,nprimal
        nobj(i)=0.0d0
        hes(i)=0.0d0
        enddo
        do i=1,ndual
        con(i)=0.0d0
        aptr(i)=0
        enddo
        aptr(ndual+1)=0
        do i=1,nnz
        jac(i)=0.0d0
        enddo

        ! geometry
        nex = 15*mult
        ney = 5*mult

        ! connect
        call connectroutine(prob,mult,nsta,ndual,connect,
     &  rowcounter,nodes,ground)
        aptr(1) = 1   
        do i=1,ndual
        aptr(i+1)=aptr(i)+rowcounter(i)
        enddo
        if(aptr(ndual+1).ne.nnz+1) then
        write(*,*) 'ERROR: aptr',aptr(ndual+1),nnz+1
        endif
    
        ! force
        call forceroutine(mult,prob,ground,fvec,fxy)

        ! FE loop
        do ey = 1, ney
        do ex = 1, nex
        ep = (ey-1)*nex+ex

        ! get mapping vector
        call maproutine(ex,ey,nex,ney,map)

        if(ground.eq.1) then
        if(ey.eq.mult*5) then
        xdual(map(7)) = 0.d0
        xdual(map(5)) = 0.0d0
        xdual(map(6)) = 0.0d0
        xdual(map(8)) = 0.d0
        xdual(map(13)) = 0.0d0
        xdual(map(14)) = 0.0d0
        endif
        elseif(ground.eq.2) then
        if(ex.eq.1) then
        xdual(map(7)) = 0.d0
        xdual(map(15)) = 0.d0
        xdual(map(1)) = 0.d0
        endif
        if(ey.le.mult/2.and.ex.eq.nex) then
        xdual(map(4)) = 0.d0
        xdual(map(12)) = 0.d0
        endif
        endif


        ! objective
        obj = obj + xprimal(nsta + ep)
        nobj(nsta+ep) = 1.0d0

        ! forces
        if(ex.ge.fxy(1,1).and.ex.le.fxy(2,1)) then
        if(ey.ge.fxy(1,2).and.ey.le.fxy(2,2)) then
        do i = 1,16
        con(map(i)) = con(map(i)) -  fvec(i)
        enddo
        endif
        endif

        ! violations
        do i = 1,16
        do j = 1,16
        call stiffroutine(i, j, ke)                  ! AG2017 , E)

        con(map(i))=con(map(i))+
     &  ke*xprimal(map(j))*(xprimal(nsta+ep))**p

        hes(nsta+ep) = hes(nsta+ep)
     &  +ke*xprimal(map(j))*p*(p-1.0d0)
     &  *xprimal(nsta+ep)**(p-2.0d0)
     &  *xdual(map(i))
        enddo
!       enddo

        !jacobian
        x_elm = xprimal(nsta+ep)
!       do i = 1, 16                 
        nsf = Aptr(map(elmd(i))) - 1
        if(rowcounter(map(elmd(i))).eq.17) then
            do j = 1, 16
            call stiffroutine(elmd(i), elmd(j), ke)  ! AG2017 , E)
            jac(nsf+j)     = jac(nsf+j)+
     &      ke*x_elm**p
            acol(nsf+j)   = map(elmd(j))
            jac(nsf+17)   = jac(nsf+17)
     &      +ke*xprimal(map(elmd(j)))
     &         *p*x_elm**(p-1.d0)
            enddo 
            acol(nsf+17)  = nsta + ep
        elseif(rowcounter(map(elmd(i))).eq.28) then
            do j = 1,26
            contmp(connect(map(elmd(i)), j)) = j
            enddo
            do j = 1,16
            call stiffroutine(elmd(i), elmd(j), ke)  ! AG2017 , E)
      jac(nsf+contmp(map(elmd(j))))=jac(nsf+contmp(map(elmd(j))))
     &      +ke*x_elm**p
            acol(nsf+contmp(map(elmd(j))))  = map(elmd(j))
            jac(nsf+27+nodes(map(elmd(i))))  =  
     &      jac(nsf+27+nodes(map(elmd(i))))  
     &      + ke*xprimal(map(elmd(j)))
     &         *p*x_elm**(p-1.d0)
            enddo
            acol(nsf+27+nodes(map(elmd(i)))) = nsta + ep
            nodes(map(elmd(i))) = nodes(map(elmd(i))) + 1
        elseif(rowcounter(map(elmd(i))).eq.46) then
         do j = 1,42
            contmp(connect(map(elmd(i)), j)) = j
            enddo
            do j = 1,16
            call stiffroutine(elmd(i), elmd(j), ke)  ! AG2017 , E)
      jac(nsf+contmp(map(elmd(j))))=jac(nsf+contmp(map(elmd(j))))
     &       +ke*x_elm**p
            acol(nsf+contmp(map(elmd(j))))  = map(elmd(j))
            jac(nsf+43+nodes(map(elmd(i))))  =  
     &      jac(nsf+43+nodes(map(elmd(i))))  
     &      + ke*xprimal(map(elmd(j)))
     &       *p*x_elm**(p-1.d0)
            enddo
            acol(nsf+43+nodes(map(elmd(i)))) = nsta + ep
            nodes(map(elmd(i))) = nodes(map(elmd(i))) + 1
        endif
        enddo

        ! BC
        if(ground.eq.1) then
        if(ey.eq.mult*5) then
        xl(map(7)) = 0.d0
        xu(map(7)) = 0.d0
        xl(map(5)) = 0.0d0
        xu(map(5)) = 0.0d0
        xl(map(6)) = 0.0d0
        xu(map(6)) = 0.0d0
        xl(map(8)) = 0.d0
        xu(map(8)) = 0.d0
        xl(map(13)) = 0.0d0
        xl(map(14)) = 0.0d0
        xu(map(13)) = 0.0d0
        xu(map(14)) = 0.0d0
        endif
        if(ey.eq.mult*5) then
        con(map(7)) = 0.d0
        con(map(5)) = 0.0d0
        con(map(6)) = 0.0d0
        con(map(8)) = 0.d0
        con(map(13)) = 0.0d0
        con(map(14)) = 0.0d0
        endif
        bmap(1) = 7
        bmap(2) = 8
        bmap(3) = 5
        bmap(4) = 6
        bmap(5) = 13
        bmap(6) = 14
        if(ey.eq.mult*5) then
        do i = 1,6
        do j = 1,aptr(map(bmap(i))+1)-aptr(map(bmap(i))) 
        jac(aptr(map(bmap(i)))-1+j) = 0.d0
        enddo
        enddo
        endif
        elseif(ground.eq.2) then
        if(ex.eq.1) then
        xl(map(7)) = 0.d0
        xl(map(15)) = 0.d0
        xl(map(1)) = 0.d0
        xu(map(7)) = 0.d0
        xu(map(15)) = 0.d0
        xu(map(1)) = 0.d0
        con(map(7)) = 0.d0
        con(map(15)) = 0.d0
        con(map(1)) = 0.d0
        endif
        if(ey.le.mult/2.and.ex.eq.nex) then
        xl(map(4)) = 0.d0
        xu(map(4)) = 0.d0
        con(map(4)) = 0.d0
        xl(map(12)) = 0.d0
        xu(map(12)) = 0.d0
        con(map(12)) = 0.d0
        endif
        if(ex.eq.1) then
        do j = 1,rowcounter(map(1))
        jac(aptr(map(1))-1+j) = 0.d0
        enddo
        do j = 1,rowcounter(map(7))  !7
        jac(aptr(map(7))-1+j) = 0.d0
        enddo
        do j = 1,rowcounter(map(15)) !15
        jac(aptr(map(15))-1+j) = 0.d0
        enddo
        endif
        if(ey.le.mult/2.and.ex.eq.nex) then
!       do j = 1,rowcounter(map(5))
!       jac(aptr(map(5))-1+j) = 0.d0
!       enddo
!       do j = 1,rowcounter(map(6))
!       jac(aptr(map(6))-1+j) = 0.d0
!       enddo
!       do j = 1,rowcounter(map(3))
!       jac(aptr(map(3))-1+j) = 0.d0
!       enddo
        do j = 1,rowcounter(map(4))
        jac(aptr(map(4))-1+j) = 0.d0
        enddo
        do j = 1,rowcounter(map(12))
        jac(aptr(map(12))-1+j) = 0.d0
        enddo
        endif
        endif
    
        call stressroutine_1(xprimal,xdual,
     &  map,nu,sy,l,nsta,nmat,nnz,ep,
     &  con,jac,hes,acol,stres)
        enddo
        enddo
        ! FE loop

        viol = 0.0d0
        do i = 1,ndual
        if(eqn(i)) then
        viol = max(viol,abs(con(i)))
        else
        viol = max(viol,max(con(i),0.0d0))
        endif
        enddo

        return
        end subroutine
        subroutine simu_2 (mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
     &               xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
     &               xl,xu,p,viol)
        implicit none
        integer          i, j, k
        integer          mult,ground,prob,nsta,nmat,nprimal,ndual,nnz
        integer          nex,ney,ep,ex,ey,map(16),fxy(2,2),nsf
        double precision xprimal(nprimal),nobj(nprimal)
        double precision obj,con(ndual),jac(nnz),hes(nprimal),p,ke,sy
        double precision u5,u6,u7,u8,v5,v6,v7,v8,sig(3),E,nu,l,fvec(16)
        integer          connect(nsta,42),nodes(nsta),acol(nnz)
        integer          rowcounter(ndual),elmd(16),aptr(ndual+1)
        integer          contmp(nsta), bmap(6)
        double precision x_elm, xdual(ndual), xl(nprimal), xu(nprimal)
        logical          eqn(ndual)
        double precision viol

!       E=1.0d0
        nu=0.3d0
        l=3.0d0/15.0d0/mult
        sy=20.0d0

        elmd(1) = 1
        elmd(2) = 2
        elmd(3) = 9
        elmd(4) = 10
        elmd(5) = 3
        elmd(6) = 4
        elmd(7) = 15
        elmd(8) = 16
        elmd(9) = 11
        elmd(10) = 12
        elmd(11) = 7
        elmd(12) = 8
        elmd(13) = 13
        elmd(14) = 14
        elmd(15) = 5
        elmd(16) = 6

        ! safety
        obj=0.0d0
        do i=1,nprimal
        nobj(i)=0.0d0
        hes(i)=0.0d0
        enddo
        do i=1,ndual
        con(i)=0.0d0
        aptr(i)=0
        enddo
        aptr(ndual+1)=0
        do i=1,nnz
        jac(i)=0.0d0
        enddo

        ! geometry
        nex = 15*mult
        ney = 5*mult
        ! lbracket specific
        if(ground.eq.3) then
        nex = 35*mult
        ney = 35*mult
        do ey = 1, ney
        do ex = 1, nex
        ep = (ey-1)*nex+ex
        if(ex.gt.14*mult.and.ey.gt.14*mult) then
        xl(nsta+ep) = 0.0d0
        xu(nsta+ep) = 0.0d0
        xprimal(nsta+ep) = 0.0d0
        endif
        enddo
        enddo
        endif

        ! connect
        call connectroutine(prob,mult,nsta,ndual,connect,
     &  rowcounter,nodes,ground)
        aptr(1) = 1   
        do i=1,ndual
        aptr(i+1)=aptr(i)+rowcounter(i)
        enddo
        if(aptr(ndual+1).ne.nnz+1) then
        write(*,*) 'ERROR: aptr',aptr(ndual+1),nnz+1
        endif
    
        ! force
        call forceroutine(mult,prob,ground,fvec,fxy)

        ! FE loop
        do ey = 1, ney
        do ex = 1, nex
        ep = (ey-1)*nex+ex

        ! get mapping vector
        call maproutine(ex,ey,nex,ney,map)

        ! forces and objective
        if(ex.ge.fxy(1,1).and.ex.le.fxy(2,1)) then
        if(ey.ge.fxy(1,2).and.ey.le.fxy(2,2)) then
        do i = 1,16
        obj = obj + xprimal(map(i))*fvec(i)
        nobj(map(i)) = fvec(i)
        con(map(i)) = con(map(i)) -  fvec(i)
        enddo
        endif
        endif

        ! violations
        do i = 1,16
        do j = 1,16
        call stiffroutine(i, j, ke)

        con(map(i))=con(map(i))+
     &  ke*xprimal(map(j))*(xprimal(nsta+ep))**p

        hes(nsta+ep) = hes(nsta+ep)
     &  +ke*xprimal(map(j))*p*(p-1.0d0)
     &  *xprimal(nsta+ep)**(p-2.0d0)
     &  *xdual(map(i))
        enddo
!       enddo

        !jacobian
        x_elm = xprimal(nsta+ep)
!       do i = 1, 16                 
        nsf = Aptr(map(elmd(i))) - 1
        if(rowcounter(map(elmd(i))).eq.17) then
            do j = 1, 16
            call stiffroutine(elmd(i), elmd(j), ke)
            jac(nsf+j)     = jac(nsf+j)+
     &      ke*x_elm**p
            acol(nsf+j)   = map(elmd(j))
            jac(nsf+17)   = jac(nsf+17)
     &      +ke*xprimal(map(elmd(j)))
     &         *p*x_elm**(p-1.d0)
            enddo 
            acol(nsf+17)  = nsta + ep
        elseif(rowcounter(map(elmd(i))).eq.28) then
            do j = 1,26
            contmp(connect(map(elmd(i)), j)) = j
            enddo
            do j = 1,16
            call stiffroutine(elmd(i), elmd(j), ke)
      jac(nsf+contmp(map(elmd(j))))=jac(nsf+contmp(map(elmd(j))))
     &      +ke*x_elm**p
            acol(nsf+contmp(map(elmd(j))))  = map(elmd(j))
            jac(nsf+27+nodes(map(elmd(i))))  =  
     &      jac(nsf+27+nodes(map(elmd(i))))  
     &      + ke*xprimal(map(elmd(j)))
     &         *p*x_elm**(p-1.d0)
            enddo
            acol(nsf+27+nodes(map(elmd(i)))) = nsta + ep
            nodes(map(elmd(i))) = nodes(map(elmd(i))) + 1
        elseif(rowcounter(map(elmd(i))).eq.46) then
         do j = 1,42
            contmp(connect(map(elmd(i)), j)) = j
            enddo
            do j = 1,16
            call stiffroutine(elmd(i), elmd(j), ke)
      jac(nsf+contmp(map(elmd(j))))=jac(nsf+contmp(map(elmd(j))))
     &       +ke*x_elm**p
            acol(nsf+contmp(map(elmd(j))))  = map(elmd(j))
            jac(nsf+43+nodes(map(elmd(i))))  =  
     &      jac(nsf+43+nodes(map(elmd(i))))  
     &      + ke*xprimal(map(elmd(j)))
     &       *p*x_elm**(p-1.d0)
            enddo
            acol(nsf+43+nodes(map(elmd(i)))) = nsta + ep
            nodes(map(elmd(i))) = nodes(map(elmd(i))) + 1
        endif
        enddo

        con(nsta+1) = con(nsta+1) + xprimal(nsta+ep)
        jac(nnz - nmat + ep) = 1.d0
        acol(nnz - nmat + ep) = nsta+ep

        ! BC
        if(ground.eq.1) then

            if(ey.eq.mult*5) then
            xl(map(7)) = 0.d0
            xu(map(7)) = 0.d0
            xl(map(5)) = 0.0d0
            xu(map(5)) = 0.0d0
            xl(map(6)) = 0.0d0
            xu(map(6)) = 0.0d0
            xl(map(8)) = 0.d0
            xu(map(8)) = 0.d0
            xl(map(13)) = 0.0d0
            xl(map(14)) = 0.0d0
            xu(map(13)) = 0.0d0
            xu(map(14)) = 0.0d0
            con(map(7)) = 0.d0
            con(map(5)) = 0.0d0
            con(map(6)) = 0.0d0
            con(map(8)) = 0.d0
            con(map(13)) = 0.0d0
            con(map(14)) = 0.0d0
            endif
        
            bmap(1) = 7
            bmap(2) = 8
            bmap(3) = 5
            bmap(4) = 6
            bmap(5) = 13
            bmap(6) = 14
            if(ey.eq.mult*5) then
            do i = 1,6
            do j = 1,aptr(map(bmap(i))+1)-aptr(map(bmap(i))) 
            jac(aptr(map(bmap(i)))-1+j) = 0.d0
            enddo
            enddo
            endif
    
        elseif(ground.eq.2) then

            if(ex.eq.1) then
            xl(map(7)) = 0.d0
            xu(map(7)) = 0.d0
            xl(map(15)) = 0.d0
            xu(map(15)) = 0.d0
            xl(map(1)) = 0.d0
            xu(map(1)) = 0.d0
            con(map(7)) = 0.0d0
            con(map(15)) = 0.0d0
            con(map(1)) = 0.0d0
            endif
            if(ey.eq.1.and.ex.eq.nex) then
            xl(map(4)) = 0.d0
            xu(map(4)) = 0.d0
            con(map(4)) = 0.0d0
            endif

            if(ex.eq.1) then
            do j = 1,rowcounter(map(1))
            jac(aptr(map(1))-1+j) = 0.d0
            enddo
            do j = 1,rowcounter(map(7))  !7
            jac(aptr(map(7))-1+j) = 0.d0
            enddo
            do j = 1,rowcounter(map(15)) !15
            jac(aptr(map(15))-1+j) = 0.d0
            enddo
            endif
            if(ey.eq.1.and.ex.eq.nex) then
            do j = 1,rowcounter(map(4))
            jac(aptr(map(4))-1+j) = 0.d0
            enddo
            endif
        
        elseif(ground.eq.3) then

        if(ey.eq.35*mult) then
        xl(map(7)) = 0.d0
        xu(map(7)) = 0.d0
        xl(map(5)) = 0.0d0
        xu(map(5)) = 0.0d0
        xl(map(6)) = 0.0d0
        xu(map(6)) = 0.0d0
        xl(map(8)) = 0.d0
        xu(map(8)) = 0.d0
        xl(map(13)) = 0.0d0
        xl(map(14)) = 0.0d0
        xu(map(13)) = 0.0d0
        xu(map(14)) = 0.0d0
        endif
        if(ey.eq.35*mult) then
        con(map(7)) = 0.d0
        con(map(5)) = 0.0d0
        con(map(6)) = 0.0d0
        con(map(8)) = 0.d0
        con(map(13)) = 0.0d0
        con(map(14)) = 0.0d0
        endif
        bmap(1) = 7
        bmap(2) = 8
        bmap(3) = 5
        bmap(4) = 6
        bmap(5) = 13
        bmap(6) = 14
        if(ey.eq.35*mult) then
        do i = 1,6
        do j = 1,aptr(map(bmap(i))+1)-aptr(map(bmap(i))) 
        jac(aptr(map(bmap(i)))-1+j) = 0.d0
        enddo
        enddo
        endif

        endif

        enddo
        enddo
        ! FE loop
        con(nsta+1) = con(nsta+1) - 0.5d0*nmat

        viol = 0.0d0
        do i = 1,ndual
        if(eqn(i)) then
        viol = max(viol,abs(con(i)))
        else
        viol = max(viol,max(con(i),0.0d0))
        endif
        enddo

        return
        end subroutine
        subroutine simu_3(mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
     &          xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
     &          xl,xu,p,bs,filt,viol)
        implicit none
        integer          i, j, k
        integer          mult,ground,prob,nsta,nmat,nprimal,ndual,nnz
        integer          nex,ney,ep,ex,ey,map(16),fxy(2,2),nsf
        double precision xprimal(nprimal),nobj(nprimal)
        double precision obj,con(ndual),jac(nnz),hes(nprimal),p,ke,sy
        double precision u5,u6,u7,u8,v5,v6,v7,v8,sig(3),E,nu,l,fvec(16)
        integer          connect(nsta,42),nodes(nsta),acol(nnz), soc1
        integer          rowcounter(ndual),elmd(16),aptr(ndual+1), soc2
        integer          contmp(nsta),bmap(6),bs,bsc,bscc,bscon(bs,2)
        double precision x_elm, xdual(ndual), xl(nprimal), xu(nprimal)
        double precision filt, viol
        logical          eqn(ndual)

!       E=1.0d0
        nu=0.3d0
        l=3.0d0/15.0d0/mult
        sy=20.0d0

        elmd(1) = 1
        elmd(2) = 2
        elmd(3) = 9
        elmd(4) = 10
        elmd(5) = 3
        elmd(6) = 4
        elmd(7) = 15
        elmd(8) = 16
        elmd(9) = 11
        elmd(10) = 12
        elmd(11) = 7
        elmd(12) = 8
        elmd(13) = 13
        elmd(14) = 14
        elmd(15) = 5
        elmd(16) = 6

        ! safety
        obj=0.0d0
        do i=1,nprimal
        nobj(i)=0.0d0
        hes(i)=0.0d0
        enddo
        do i=1,ndual
        con(i)=0.0d0
        aptr(i)=0
        enddo
        aptr(ndual+1)=0
        do i=1,nnz
        jac(i)=0.0d0
        enddo

        ! geometry
        nex = 15*mult
        ney = 5*mult
        ! lbracket specific
        if(ground.eq.3) then
        nex = 35*mult
        ney = 35*mult
        do ey = 1, ney
        do ex = 1, nex
        ep = (ey-1)*nex+ex
        if(ex.gt.14*mult.and.ey.gt.14*mult) then
        xl(nsta+ep) = 0.0d0
        xu(nsta+ep) = 0.0d0
        xprimal(nsta+ep) = 0.0d0
        endif
        enddo
        enddo
        endif

        ! connect
        call connectroutine(prob,mult,nsta,ndual,connect,
     &  rowcounter,nodes,ground)
        aptr(1) = 1   
        do i=1,ndual
        aptr(i+1)=aptr(i)+rowcounter(i)
        enddo
        if(aptr(ndual+1).ne.nnz+1) then
        write(*,*) 'ERROR: aptr',aptr(ndual+1),nnz+1
        endif
    
        ! force
        call forceroutine(mult,prob,ground,fvec,fxy)

        ! FE loop
        do ey = 1, ney
        do ex = 1, nex
        ep = (ey-1)*nex+ex

        ! get mapping vector
        call maproutine(ex,ey,nex,ney,map)

        ! forces and objective
        if(ex.ge.fxy(1,1).and.ex.le.fxy(2,1)) then
        if(ey.ge.fxy(1,2).and.ey.le.fxy(2,2)) then
        do i = 1,16
        obj = obj + xprimal(map(i))*fvec(i)
        nobj(map(i)) = fvec(i)
        con(map(i)) = con(map(i)) -  fvec(i)
        enddo
        endif
        endif

        ! violations
        do i = 1,16
        do j = 1,16
        call stiffroutine(i, j, ke)

        con(map(i))=con(map(i))+
     &  ke*xprimal(map(j))*(xprimal(nsta+ep))**p

        hes(nsta+ep) = hes(nsta+ep)
     &  +ke*xprimal(map(j))*p*(p-1.0d0)
     &  *xprimal(nsta+ep)**(p-2.0d0)
     &  *xdual(map(i))
        enddo
!       enddo

        !jacobian
        x_elm = xprimal(nsta+ep)
!       do i = 1, 16                 
        nsf = Aptr(map(elmd(i))) - 1
        if(rowcounter(map(elmd(i))).eq.17) then
            do j = 1, 16
            call stiffroutine(elmd(i), elmd(j), ke)
            jac(nsf+j)     = jac(nsf+j)+
     &      ke*x_elm**p
            acol(nsf+j)   = map(elmd(j))
            jac(nsf+17)   = jac(nsf+17)
     &      +ke*xprimal(map(elmd(j)))
     &         *p*x_elm**(p-1.d0)
            enddo 
            acol(nsf+17)  = nsta + ep
        elseif(rowcounter(map(elmd(i))).eq.28) then
            do j = 1,26
            contmp(connect(map(elmd(i)), j)) = j
            enddo
            do j = 1,16
            call stiffroutine(elmd(i), elmd(j), ke)
      jac(nsf+contmp(map(elmd(j))))=jac(nsf+contmp(map(elmd(j))))
     &      +ke*x_elm**p
            acol(nsf+contmp(map(elmd(j))))  = map(elmd(j))
            jac(nsf+27+nodes(map(elmd(i))))  =  
     &      jac(nsf+27+nodes(map(elmd(i))))  
     &      + ke*xprimal(map(elmd(j)))
     &         *p*x_elm**(p-1.d0)
            enddo
            acol(nsf+27+nodes(map(elmd(i)))) = nsta + ep
            nodes(map(elmd(i))) = nodes(map(elmd(i))) + 1
        elseif(rowcounter(map(elmd(i))).eq.46) then
         do j = 1,42
            contmp(connect(map(elmd(i)), j)) = j
            enddo
            do j = 1,16
            call stiffroutine(elmd(i), elmd(j), ke)
      jac(nsf+contmp(map(elmd(j))))=jac(nsf+contmp(map(elmd(j))))
     &       +ke*x_elm**p
            acol(nsf+contmp(map(elmd(j))))  = map(elmd(j))
            jac(nsf+43+nodes(map(elmd(i))))  =  
     &      jac(nsf+43+nodes(map(elmd(i))))  
     &      + ke*xprimal(map(elmd(j)))
     &       *p*x_elm**(p-1.d0)
            enddo
            acol(nsf+43+nodes(map(elmd(i)))) = nsta + ep
            nodes(map(elmd(i))) = nodes(map(elmd(i))) + 1
        endif
        enddo

        con(nsta+1) = con(nsta+1) + xprimal(nsta+ep)
        jac(nnz - 4*bs - nmat + ep) = 1.d0
        acol(nnz - 4*bs - nmat + ep) = nsta+ep

        ! BC
        if(ground.eq.1) then

            if(ey.eq.mult*5) then
            xl(map(7)) = 0.d0
            xu(map(7)) = 0.d0
            xl(map(5)) = 0.0d0
            xu(map(5)) = 0.0d0
            xl(map(6)) = 0.0d0
            xu(map(6)) = 0.0d0
            xl(map(8)) = 0.d0
            xu(map(8)) = 0.d0
            xl(map(13)) = 0.0d0
            xl(map(14)) = 0.0d0
            xu(map(13)) = 0.0d0
            xu(map(14)) = 0.0d0
            con(map(7)) = 0.d0
            con(map(5)) = 0.0d0
            con(map(6)) = 0.0d0
            con(map(8)) = 0.d0
            con(map(13)) = 0.0d0
            con(map(14)) = 0.0d0
            endif
        
            bmap(1) = 7
            bmap(2) = 8
            bmap(3) = 5
            bmap(4) = 6
            bmap(5) = 13
            bmap(6) = 14
            if(ey.eq.mult*5) then
            do i = 1,6
            do j = 1,aptr(map(bmap(i))+1)-aptr(map(bmap(i))) 
            jac(aptr(map(bmap(i)))-1+j) = 0.d0
            enddo
            enddo
            endif
    
        elseif(ground.eq.2) then

            if(ex.eq.1) then
            xl(map(7)) = 0.d0
            xu(map(7)) = 0.d0
            xl(map(15)) = 0.d0
            xu(map(15)) = 0.d0
            xl(map(1)) = 0.d0
            xu(map(1)) = 0.d0
            con(map(7)) = 0.0d0
            con(map(15)) = 0.0d0
            con(map(1)) = 0.0d0
            endif
            if(ey.eq.1.and.ex.eq.nex) then
            xl(map(4)) = 0.d0
            xu(map(4)) = 0.d0
            con(map(4)) = 0.0d0
            endif

            if(ex.eq.1) then
            do j = 1,rowcounter(map(1))
            jac(aptr(map(1))-1+j) = 0.d0
            enddo
            do j = 1,rowcounter(map(7))  !7
            jac(aptr(map(7))-1+j) = 0.d0
            enddo
            do j = 1,rowcounter(map(15)) !15
            jac(aptr(map(15))-1+j) = 0.d0
            enddo
            endif
            if(ey.eq.1.and.ex.eq.nex) then
            do j = 1,rowcounter(map(4))
            jac(aptr(map(4))-1+j) = 0.d0
            enddo
            endif
        elseif(ground.eq.3) then

        if(ey.eq.35*mult) then
        xl(map(7)) = 0.d0
        xu(map(7)) = 0.d0
        xl(map(5)) = 0.0d0
        xu(map(5)) = 0.0d0
        xl(map(6)) = 0.0d0
        xu(map(6)) = 0.0d0
        xl(map(8)) = 0.d0
        xu(map(8)) = 0.d0
        xl(map(13)) = 0.0d0
        xl(map(14)) = 0.0d0
        xu(map(13)) = 0.0d0
        xu(map(14)) = 0.0d0
        endif
        if(ey.eq.35*mult) then
        con(map(7)) = 0.d0
        con(map(5)) = 0.0d0
        con(map(6)) = 0.0d0
        con(map(8)) = 0.d0
        con(map(13)) = 0.0d0
        con(map(14)) = 0.0d0
        endif
        bmap(1) = 7
        bmap(2) = 8
        bmap(3) = 5
        bmap(4) = 6
        bmap(5) = 13
        bmap(6) = 14
        if(ey.eq.35*mult) then
        do i = 1,6
        do j = 1,aptr(map(bmap(i))+1)-aptr(map(bmap(i))) 
        jac(aptr(map(bmap(i)))-1+j) = 0.d0
        enddo
        enddo
        endif
        endif

        enddo
        enddo
        ! FE loop
        con(nsta+1) = con(nsta+1) - 0.5d0*nmat

        ! Slope constraints
        bsc = 1
        bscc = 1
        do j = 1, ney
        do i = 1, nex-1
        bscon(bsc,1) = bscc
        bscon(bsc,2) = bscc+1
        bsc = bsc + 1
        bscc = bscc + 1
        enddo
        bscc = bscc + 1
        enddo
        bsc = ney*(nex-1) + 1
        bscc = 1
        do j = 1, ney-1
        do i = 1, nex
        bscon(bsc,1) = bscc
        bscon(bsc,2) = bscc + nex
        bsc = bsc + 1
        bscc = bscc + 1
        enddo
        enddo

        do i = 1, bs
        if(xl(nsta+bscon(i,1)).eq.xu(nsta+bscon(i,1)).and.
     &  ground.eq.3) then     
        con(nsta+1+i) =                                         
     &        -filt/DBLE(nex)                                   
        elseif(xl(nsta+bscon(i,2)).eq.xu(nsta+bscon(i,2)).and.
     &  ground.eq.3) then 
        con(nsta+1+i) =                                         
     &        -filt/DBLE(nex)                                   
        else                                                    
        con(nsta+1+i) =
     &        xprimal(nsta+bscon(i,1))-xprimal(nsta+bscon(i,2))
     &        -filt/DBLE(nex)
        endif                                                   
        enddo
        do i = 1, bs

        if(xl(nsta+bscon(i,1)).eq.xu(nsta+bscon(i,1)).and.
     &  ground.eq.3) then     
        con(nsta+1+bs+i) =                                      
     &        -filt/DBLE(nex)                                   
        elseif(xl(nsta+bscon(i,2)).eq.xu(nsta+bscon(i,2)).and.
     &  ground.eq.3) then 
        con(nsta+1+bs+i) =                                      
     &        -filt/DBLE(nex)                                   
        else                                                    
        con(nsta+1+bs+i) =
     &        xprimal(nsta+bscon(i,2))-xprimal(nsta+bscon(i,1))
     &        -filt/DBLE(nex)
        endif                                                   
        enddo

        soc1 = 1
        soc2 = 2
        do i = 1, bs
        jac(nnz-4*bs+soc1) = 1.0d0
        jac(nnz-4*bs+soc2) = -1.0d0
        acol(nnz-4*bs+soc1) = nsta+bscon(i,1)
        acol(nnz-4*bs+soc2) = nsta+bscon(i,2)
        soc1 = soc1 + 2
        soc2 = soc2 + 2
        enddo
!
        soc1 = 1
        soc2 = 2
        do i = 1, bs
        jac(nnz-2*bs+soc1) = -1.0d0
        jac(nnz-2*bs+soc2) = 1.0d0
        acol(nnz-2*bs+soc1) = nsta+bscon(i,1)
        acol(nnz-2*bs+soc2) = nsta+bscon(i,2)
        soc1 = soc1 + 2
        soc2 = soc2 + 2
        enddo

        viol = 0.0d0
        do i = 1,ndual
        if(eqn(i)) then
        viol = max(viol,abs(con(i)))
        else
        viol = max(viol,max(con(i),0.0d0))
        endif
        enddo

        return
        end subroutine
        subroutine simu_4(mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
     &          xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
     &          xl,xu,p,bs,filt,viol,stres,xprimalh,conh,
     &          turvj)
        implicit none
        integer          i, j, k
        integer          mult,ground,prob,nsta,nmat,nprimal,ndual,nnz
        integer          nex,ney,ep,ex,ey,map(16),fxy(2,2),nsf
        double precision xprimal(nprimal),nobj(nprimal)
        double precision xprimalh(nprimal), conh(ndual)
        double precision obj,con(ndual),jac(nnz),hes(nprimal),p,ke,sy
        double precision u5,u6,u7,u8,v5,v6,v7,v8,sig(3),E,nu,l,fvec(16)
        integer          connect(nsta,42),nodes(nsta),acol(nnz), soc1
        integer          rowcounter(ndual),elmd(16),aptr(ndual+1), soc2
        integer          contmp(nsta),bmap(6),bs,bsc,bscc,bscon(bs,2)
        double precision x_elm, xdual(ndual), xl(nprimal), xu(nprimal)
        double precision filt,viol,stres(nmat),turvj(ndual)
        logical          eqn(ndual)

!       E=1.0d0
        nu=0.3d0
        l=3.0d0/15.0d0/mult
        sy=20.0d0

        do i =1,nmat
        if(xprimal(nsta+i).le.1d-2) then
        xprimal(nsta+i)=0.0d0
        endif
        enddo

        elmd(1) = 1
        elmd(2) = 2
        elmd(3) = 9
        elmd(4) = 10
        elmd(5) = 3
        elmd(6) = 4
        elmd(7) = 15
        elmd(8) = 16
        elmd(9) = 11
        elmd(10) = 12
        elmd(11) = 7
        elmd(12) = 8
        elmd(13) = 13
        elmd(14) = 14
        elmd(15) = 5
        elmd(16) = 6

        ! safety
        obj=0.0d0
        do i=1,nprimal
        nobj(i)=0.0d0
        hes(i)=0.0d0
        enddo
        do i=1,ndual
        turvj(i)=conh(i)*xdual(i)
        con(i)=0.0d0
        aptr(i)=0
        enddo
        aptr(ndual+1)=0
        do i=1,nnz
        jac(i)=0.0d0
        enddo

        ! geometry
        nex = 15*mult
        ney = 5*mult
!       ! mbb specific
        if (ground.eq.2) then
        do ey = 1, mult/2
        do ex = 1, nex
        ep = (ey-1)*nex+ex
        if (ex.gt.(nex-mult/2)) then
        xl(nsta+ep) = 1.0d0
        xu(nsta+ep) = 1.0d0
        xprimal(nsta+ep) = 1.0d0
        endif
        enddo
        enddo
        endif
        ! lbracket specific
        if(ground.eq.3) then
        nex = 25*mult
        ney = 25*mult
        l=10.0d0/25.0d0/mult
        sy=4.0d0
        do ey = 1, ney
        do ex = 1, nex
        ep = (ey-1)*nex+ex
        if(ex.gt.10*mult.and.ey.gt.10*mult) then
        xl(nsta+ep) = 0.0d0
        xu(nsta+ep) = 0.0d0
        xprimal(nsta+ep) = 0.0d0
        endif
        if(ex.gt.nex-mult/2) then
        if(ey.ge.4*ney/10/2-mult/2+1.and.ey.le.4*ney/10/2+mult/2) then!1    !y-coor-1
        xl(nsta+ep) = 1.0d0
        xu(nsta+ep) = 1.0d0
        xprimal(nsta+ep) = 1.0d0
        endif
        endif
        enddo
        enddo
        endif

        ! connect
        call connectroutine(prob,mult,nsta,ndual,connect,
     &  rowcounter,nodes,ground)
        aptr(1) = 1   
        do i=1,ndual
        aptr(i+1)=aptr(i)+rowcounter(i)
        enddo
        if(aptr(ndual+1).ne.nnz+1) then
        write(*,*) 'ERROR: aptr',aptr(ndual+1),nnz+1
        endif
    
        ! force
        call forceroutine(mult,prob,ground,fvec,fxy)

        ! FE loop
        do ey = 1, ney
        do ex = 1, nex
        ep = (ey-1)*nex+ex

        ! get mapping vector
        call maproutine(ex,ey,nex,ney,map)

        ! objective
        obj = obj + xprimal(nsta + ep)
        nobj(nsta+ep) = 1.0d0
        ! forces
        if(ex.ge.fxy(1,1).and.ex.le.fxy(2,1)) then
        if(ey.ge.fxy(1,2).and.ey.le.fxy(2,2)) then
        do i = 1,16
        con(map(i)) = con(map(i)) -  fvec(i)
        turvj(map(i))=turvj(map(i))+fvec(i)*xdual(map(i))
        enddo
        endif
        endif

        ! violations
        do i = 1,16
        do j = 1,16
        call stiffroutine(i, j, ke)

        con(map(i))=con(map(i))+
     &  ke*xprimal(map(j))*(xprimal(nsta+ep))**p

        turvj(map(i))=turvj(map(i))-
     &  ke*xprimal(map(j))*(xprimal(nsta+ep))**p
     &  *xdual(map(i))

        hes(nsta+ep) = hes(nsta+ep)
     &  +ke*xprimal(map(j))*p*(p-1.0d0)
     &  *xprimal(nsta+ep)**(p-2.0d0)
     &  *xdual(map(i))
        enddo
!       enddo

        !jacobian
        x_elm = xprimal(nsta+ep)
!       do i = 1, 16                 
        nsf = Aptr(map(elmd(i))) - 1
        if(rowcounter(map(elmd(i))).eq.17) then
            do j = 1, 16
            call stiffroutine(elmd(i), elmd(j), ke)
            jac(nsf+j)     = jac(nsf+j)+
     &      ke*x_elm**p
            acol(nsf+j)   = map(elmd(j))
            jac(nsf+17)   = jac(nsf+17)
     &      +ke*xprimal(map(elmd(j)))
     &         *p*x_elm**(p-1.d0)

            turvj(map(elmd(i)))=turvj(map(elmd(i)))
     &      -ke*x_elm**p
     &      *(xprimalh(map(elmd(j)))-xprimal(map(elmd(j))))
     &      *xdual(map(elmd(i)))
!           turvj(map(elmd(i)))=turvj(map(elmd(i)))
     &      -ke*xprimal(map(elmd(j)))
     &         *p*x_elm**(p-1.d0)
     &      *(xprimalh(nsta+ep)-xprimal(nsta+ep))
     &      *xdual(map(elmd(i)))

            enddo 
            acol(nsf+17)  = nsta + ep
        elseif(rowcounter(map(elmd(i))).eq.28) then
            do j = 1,26
            contmp(connect(map(elmd(i)), j)) = j
            enddo
            do j = 1,16
            call stiffroutine(elmd(i), elmd(j), ke)
      jac(nsf+contmp(map(elmd(j))))=jac(nsf+contmp(map(elmd(j))))
     &      +ke*x_elm**p
            acol(nsf+contmp(map(elmd(j))))  = map(elmd(j))
            jac(nsf+27+nodes(map(elmd(i))))  =  
     &      jac(nsf+27+nodes(map(elmd(i))))  
     &      + ke*xprimal(map(elmd(j)))
     &         *p*x_elm**(p-1.d0)

            turvj(map(elmd(i)))=turvj(map(elmd(i)))
     &      -ke*x_elm**p
     &      *(xprimalh(map(elmd(j)))-xprimal(map(elmd(j))))
     &      *xdual(map(elmd(i)))
!           turvj(map(elmd(i)))=turvj(map(elmd(i)))
     &      -ke*xprimal(map(elmd(j)))
     &         *p*x_elm**(p-1.d0)
     &      *(xprimalh(nsta+ep)-xprimal(nsta+ep))
     &      *xdual(map(elmd(i)))

            enddo
            acol(nsf+27+nodes(map(elmd(i)))) = nsta + ep
            nodes(map(elmd(i))) = nodes(map(elmd(i))) + 1
        elseif(rowcounter(map(elmd(i))).eq.46) then
         do j = 1,42
            contmp(connect(map(elmd(i)), j)) = j
            enddo
            do j = 1,16
            call stiffroutine(elmd(i), elmd(j), ke)
      jac(nsf+contmp(map(elmd(j))))=jac(nsf+contmp(map(elmd(j))))
     &       +ke*x_elm**p
            acol(nsf+contmp(map(elmd(j))))  = map(elmd(j))
            jac(nsf+43+nodes(map(elmd(i))))  =  
     &      jac(nsf+43+nodes(map(elmd(i))))  
     &      + ke*xprimal(map(elmd(j)))
     &       *p*x_elm**(p-1.d0)

            turvj(map(elmd(i)))=turvj(map(elmd(i)))
     &      -ke*x_elm**p
     &      *(xprimalh(map(elmd(j)))-xprimal(map(elmd(j))))
     &      *xdual(map(elmd(i)))
!           turvj(map(elmd(i)))=turvj(map(elmd(i)))
     &      -ke*xprimal(map(elmd(j)))
     &         *p*x_elm**(p-1.d0)
     &      *(xprimalh(nsta+ep)-xprimal(nsta+ep))
     &      *xdual(map(elmd(i)))

            enddo
            acol(nsf+43+nodes(map(elmd(i)))) = nsta + ep
            nodes(map(elmd(i))) = nodes(map(elmd(i))) + 1
        endif
        enddo

        ! BC
        if(ground.eq.1) then
         if(ey.eq.mult*5) then
         xl(map(7)) = 0.d0
         xu(map(7)) = 0.d0
         xl(map(5)) = 0.0d0
         xu(map(5)) = 0.0d0
         xl(map(6)) = 0.0d0
         xu(map(6)) = 0.0d0
         xl(map(8)) = 0.d0
         xu(map(8)) = 0.d0
         xl(map(13)) = 0.0d0
         xl(map(14)) = 0.0d0
         xu(map(13)) = 0.0d0
         xu(map(14)) = 0.0d0
         con(map(7)) = 0.d0
         con(map(5)) = 0.0d0
         con(map(6)) = 0.0d0
         con(map(8)) = 0.d0
         con(map(13)) = 0.0d0
         con(map(14)) = 0.0d0
         turvj(map(7)) = 0.d0
         turvj(map(5)) = 0.0d0
         turvj(map(6)) = 0.0d0
         turvj(map(8)) = 0.d0
         turvj(map(13)) = 0.0d0
         turvj(map(14)) = 0.0d0
         endif
       
         bmap(1) = 7
         bmap(2) = 8
         bmap(3) = 5
         bmap(4) = 6
         bmap(5) = 13
         bmap(6) = 14
         if(ey.eq.mult*5) then
         do i = 1,6
         do j = 1,aptr(map(bmap(i))+1)-aptr(map(bmap(i))) 
         jac(aptr(map(bmap(i)))-1+j) = 0.d0
         enddo
         enddo
         endif
        elseif(ground.eq.2) then
         if(ex.eq.1) then
         xl(map(7)) = 0.d0
         xu(map(7)) = 0.d0
         xl(map(15)) = 0.d0
         xu(map(15)) = 0.d0
         xl(map(1)) = 0.d0
         xu(map(1)) = 0.d0
         con(map(7)) = 0.0d0
         con(map(15)) = 0.0d0
         con(map(1)) = 0.0d0
         turvj(map(7)) = 0.0d0
         turvj(map(15)) = 0.0d0
         turvj(map(1)) = 0.0d0
         endif
         if(ey.le.mult/2.and.ex.eq.nex) then
         xl(map(4)) = 0.d0
         xu(map(4)) = 0.d0
         con(map(4)) = 0.d0
         turvj(map(4)) = 0.0d0
         xl(map(12)) = 0.d0
         xu(map(12)) = 0.d0
         con(map(12)) = 0.d0
         turvj(map(12)) = 0.d0
         endif
         if(ex.eq.1) then
         do j = 1,rowcounter(map(1))
         jac(aptr(map(1))-1+j) = 0.d0
         enddo
         do j = 1,rowcounter(map(7))  !7
         jac(aptr(map(7))-1+j) = 0.d0
         enddo
         do j = 1,rowcounter(map(15)) !15
         jac(aptr(map(15))-1+j) = 0.d0
         enddo
         endif
         if(ey.le.mult/2.and.ex.eq.nex) then
         do j = 1,rowcounter(map(4))
         jac(aptr(map(4))-1+j) = 0.d0
         enddo
         do j = 1,rowcounter(map(12))
         jac(aptr(map(12))-1+j) = 0.d0
         enddo
         endif
        elseif(ground.eq.3) then
         if(ey.eq.25*mult) then
         xl(map(7)) = 0.d0
         xu(map(7)) = 0.d0
         xl(map(5)) = 0.0d0
         xu(map(5)) = 0.0d0
         xl(map(6)) = 0.0d0
         xu(map(6)) = 0.0d0
         xl(map(8)) = 0.d0
         xu(map(8)) = 0.d0
         xl(map(13)) = 0.0d0
         xl(map(14)) = 0.0d0
         xu(map(13)) = 0.0d0
         xu(map(14)) = 0.0d0
         turvj(map(7)) = 0.d0
         turvj(map(5)) = 0.0d0
         turvj(map(6)) = 0.0d0
         turvj(map(8)) = 0.d0
         turvj(map(13)) = 0.0d0
         turvj(map(14)) = 0.0d0
         endif
         if(ey.eq.25*mult) then
         con(map(7)) = 0.d0
         con(map(5)) = 0.0d0
         con(map(6)) = 0.0d0
         con(map(8)) = 0.d0
         con(map(13)) = 0.0d0
         con(map(14)) = 0.0d0
         endif
         bmap(1) = 7
         bmap(2) = 8
         bmap(3) = 5
         bmap(4) = 6
         bmap(5) = 13
         bmap(6) = 14
         if(ey.eq.25*mult) then
         do i = 1,6
         do j = 1,aptr(map(bmap(i))+1)-aptr(map(bmap(i))) 
         jac(aptr(map(bmap(i)))-1+j) = 0.d0
         enddo
         enddo
         endif
        endif

        call stressroutine_11(xprimal,xdual,
     &  map,nu,sy,l,nsta,nmat,nnz,ep,
     &  con,jac,hes,acol,stres,bs)

        turvj(nsta+ep)=turvj(nsta+ep)-con(nsta+ep)
     &  *xdual(nsta+ep)
!       turvj(nsta+ep)=turvj(nsta+ep)
     &  -jac(nnz-4*bs-9*nmat+(ep-1)*9+1)
     &  *(xprimalh(map(9))-xprimal(map(9)))
     &  *xdual(nsta+ep)
!       turvj(nsta+ep)=turvj(nsta+ep)
     &  -jac(nnz-4*bs-9*nmat+(ep-1)*9+2)
     &  *(xprimalh(map(10))-xprimal(map(10)))
     &  *xdual(nsta+ep)
!       turvj(nsta+ep)=turvj(nsta+ep)
     &  -jac(nnz-4*bs-9*nmat+(ep-1)*9+3)
     &  *(xprimalh(map(15))-xprimal(map(15)))
     &  *xdual(nsta+ep)
!       turvj(nsta+ep)=turvj(nsta+ep)
     &  -jac(nnz-4*bs-9*nmat+(ep-1)*9+4)
     &  *(xprimalh(map(16))-xprimal(map(16)))
     &  *xdual(nsta+ep)
!       turvj(nsta+ep)=turvj(nsta+ep)
     &  -jac(nnz-4*bs-9*nmat+(ep-1)*9+5)
     &  *(xprimalh(map(11))-xprimal(map(11)))
     &  *xdual(nsta+ep)
!       turvj(nsta+ep)=turvj(nsta+ep)
     &  -jac(nnz-4*bs-9*nmat+(ep-1)*9+6)
     &  *(xprimalh(map(12))-xprimal(map(12)))
     &  *xdual(nsta+ep)
!       turvj(nsta+ep)=turvj(nsta+ep)
     &  -jac(nnz-4*bs-9*nmat+(ep-1)*9+7)
     &  *(xprimalh(map(13))-xprimal(map(13)))
     &  *xdual(nsta+ep)
!       turvj(nsta+ep)=turvj(nsta+ep)
     &  -jac(nnz-4*bs-9*nmat+(ep-1)*9+8)
     &  *(xprimalh(map(14))-xprimal(map(14)))
     &  *xdual(nsta+ep)
!       turvj(nsta+ep)=turvj(nsta+ep)
     &  -jac(nnz-4*bs-9*nmat+(ep-1)*9+9)
     &  *(xprimalh(nsta+ep)-xprimal(nsta+ep))
     &  *xdual(nsta+ep)

        if(xu(nsta+ep).eq.0.0d0.and.xl(nsta+ep).eq.0.0d0) then
        turvj(nsta+ep)=0.0d0
        con(nsta+ep)=0.0d0
        xdual(nsta+ep)=0.0d0
        stres(ep)=0.0d0
        endif

        enddo
        enddo

        ! Slope constraints
        bsc = 1
        bscc = 1
        do j = 1, ney
        do i = 1, nex-1
        bscon(bsc,1) = bscc
        bscon(bsc,2) = bscc+1
        bsc = bsc + 1
        bscc = bscc + 1
        enddo
        bscc = bscc + 1
        enddo
        bsc = ney*(nex-1) + 1
        bscc = 1
        do j = 1, ney-1
        do i = 1, nex
        bscon(bsc,1) = bscc
        bscon(bsc,2) = bscc + nex
        bsc = bsc + 1
        bscc = bscc + 1
        enddo
        enddo

        do i = 1, bs
        if(xl(nsta+bscon(i,1)).eq.xu(nsta+bscon(i,1)).and.
     &  ground.eq.3) then
        con(nsta+nmat+i) = -1d6
!    &        -filt/DBLE(nex)
        elseif(xl(nsta+bscon(i,2)).eq.xu(nsta+bscon(i,2)).and.
     &  ground.eq.3) then
        con(nsta+nmat+i) = -1d6
!    &        -filt/DBLE(nex)
        else
        con(nsta+nmat+i) =
     &        xprimal(nsta+bscon(i,1))-xprimal(nsta+bscon(i,2))
     &        -filt/DBLE(nex)
        endif
        enddo
        do i = 1, bs
        if(xl(nsta+bscon(i,1)).eq.xu(nsta+bscon(i,1)).and.
     &  ground.eq.3) then
        con(nsta+nmat+bs+i) = -1d6
!    &        -filt/DBLE(nex)
        elseif(xl(nsta+bscon(i,2)).eq.xu(nsta+bscon(i,2)).and.
     &  ground.eq.3) then
        con(nsta+nmat+bs+i) = -1d6
!    &        -filt/DBLE(nex)
        else
        con(nsta+nmat+bs+i) =
     &        xprimal(nsta+bscon(i,2))-xprimal(nsta+bscon(i,1))
     &        -filt/DBLE(nex)
        endif
        enddo

        soc1 = 1
        soc2 = 2
        do i = 1, bs
        jac(nnz-4*bs+soc1) = 1.0d0
        jac(nnz-4*bs+soc2) = -1.0d0
        acol(nnz-4*bs+soc1) = nsta+bscon(i,1)
        acol(nnz-4*bs+soc2) = nsta+bscon(i,2)
        soc1 = soc1 + 2
        soc2 = soc2 + 2
        enddo
!
        soc1 = 1
        soc2 = 2
        do i = 1, bs
        jac(nnz-2*bs+soc1) = -1.0d0
        jac(nnz-2*bs+soc2) = 1.0d0
        acol(nnz-2*bs+soc1) = nsta+bscon(i,1)
        acol(nnz-2*bs+soc2) = nsta+bscon(i,2)
        soc1 = soc1 + 2
        soc2 = soc2 + 2
        enddo

        viol = 0.0d0
        do i = 1,ndual
        if(eqn(i)) then
        viol = max(viol,abs(con(i)))
        else
        viol = max(viol,max(con(i),0.0d0))
        endif
        enddo

        return
        end subroutine
!         subroutine simu_5(mult,ground,prob,nsta,nmat,nprimal,ndual,nnz,
!      &          xprimal,xdual,eqn,obj,nobj,con,jac,hes,aptr,acol,
!      &          xl,xu,p,bs,filt,viol,stres,eps)
!         implicit none
!         integer          i, j, k
!         integer          mult,ground,prob,nsta,nmat,nprimal,ndual,nnz
!         integer          nex,ney,ep,ex,ey,map(16),fxy(2,2),nsf
!         double precision xprimal(nprimal),nobj(nprimal)
!         double precision obj,con(ndual),jac(nnz),hes(nprimal),p,ke,sy
!         double precision u5,u6,u7,u8,v5,v6,v7,v8,sig(3),E,nu,l,fvec(16)
!         integer          connect(nsta,42),nodes(nsta),acol(nnz), soc1
!         integer          rowcounter(ndual),elmd(16),aptr(ndual+1), soc2
!         integer          contmp(nsta),bmap(6),bs,bsc,bscc,bscon(bs,2)
!         double precision x_elm, xdual(ndual), xl(nprimal), xu(nprimal)
!         double precision filt, viol, eps, stres(nmat)
!         logical          eqn(ndual)
! 
! !       E=1.0d0
!         nu=0.3d0
!         l=3.0d0/15.0d0/mult
!         sy=20.0d0
! 
!         do i =1,nmat
!         if(xprimal(nsta+i).le.1d-6) then
!         xprimal(nsta+i)=0.0d0
!         endif
!         enddo
! 
!         elmd(1) = 1
!         elmd(2) = 2
!         elmd(3) = 9
!         elmd(4) = 10
!         elmd(5) = 3
!         elmd(6) = 4
!         elmd(7) = 15
!         elmd(8) = 16
!         elmd(9) = 11
!         elmd(10) = 12
!         elmd(11) = 7
!         elmd(12) = 8
!         elmd(13) = 13
!         elmd(14) = 14
!         elmd(15) = 5
!         elmd(16) = 6
! 
!         ! safety
!         obj=0.0d0
!         do i=1,nprimal
!         nobj(i)=0.0d0
!         hes(i)=0.0d0
!         enddo
!         do i=1,ndual
!         con(i)=0.0d0
!         aptr(i)=0
!         enddo
!         aptr(ndual+1)=0
!         do i=1,nnz
!         jac(i)=0.0d0
!         enddo
! 
!         ! geometry
!         nex = 15*mult
!         ney = 5*mult
!         ! lbracket specific
!         if(ground.eq.3) then
!         nex = 35*mult
!         ney = 35*mult
!         do ey = 1, ney
!         do ex = 1, nex
!         ep = (ey-1)*nex+ex
!         if(ex.gt.14*mult.and.ey.gt.14*mult) then
!         xl(nsta+ep) = 0.0d0
!         xu(nsta+ep) = 0.0d0
!         xprimal(nsta+ep) = 0.0d0
!         endif
!         enddo
!         enddo
!         endif
! 
!         ! connect
!         call connectroutine(prob,mult,nsta,ndual,connect,
!      &  rowcounter,nodes,ground)
!         aptr(1) = 1   
!         do i=1,ndual
!         aptr(i+1)=aptr(i)+rowcounter(i)
!         enddo
!         if(aptr(ndual+1).ne.nnz+1) then
!         write(*,*) 'ERROR: aptr',aptr(ndual+1),nnz+1
!         endif
!     
!         ! force
!         call forceroutine(mult,prob,ground,fvec,fxy)
! 
!         ! FE loop
!         do ey = 1, ney
!         do ex = 1, nex
!         ep = (ey-1)*nex+ex
! 
!         ! get mapping vector
!         call maproutine(ex,ey,nex,ney,map)
! 
!         ! objective
! !       obj = obj + xprimal(nsta + ep)
! !       nobj(nsta+ep) = 1.0d0
!         ! forces
!         if(ex.ge.fxy(1,1).and.ex.le.fxy(2,1)) then
!         if(ey.ge.fxy(1,2).and.ey.le.fxy(2,2)) then
!         do i = 1,16
!         con(map(i)) = con(map(i)) -  fvec(i)
!         obj = obj + xprimal(map(i))*fvec(i)
!         nobj(map(i)) = fvec(i)
!         enddo
!         endif
!         endif
! 
!         ! violations
!         do i = 1,16
!         do j = 1,16
!         call stiffroutine(i, j, ke)
! 
!         con(map(i))=con(map(i))+
!      &  ke*xprimal(map(j))*(xprimal(nsta+ep))**p
! 
!         hes(nsta+ep) = hes(nsta+ep)
!      &  +ke*xprimal(map(j))*p*(p-1.0d0)
!      &  *xprimal(nsta+ep)**(p-2.0d0)
!      &  *xdual(map(i))
!         enddo
! !       enddo
! 
!         !jacobian
!         x_elm = xprimal(nsta+ep)
! !       do i = 1, 16                 
!         nsf = Aptr(map(elmd(i))) - 1
!         if(rowcounter(map(elmd(i))).eq.17) then
!             do j = 1, 16
!             call stiffroutine(elmd(i), elmd(j), ke)
!             jac(nsf+j)     = jac(nsf+j)+
!      &      ke*x_elm**p
!             acol(nsf+j)   = map(elmd(j))
!             jac(nsf+17)   = jac(nsf+17)
!      &      +ke*xprimal(map(elmd(j)))
!      &         *p*x_elm**(p-1.d0)
!             enddo 
!             acol(nsf+17)  = nsta + ep
!         elseif(rowcounter(map(elmd(i))).eq.28) then
!             do j = 1,26
!             contmp(connect(map(elmd(i)), j)) = j
!             enddo
!             do j = 1,16
!             call stiffroutine(elmd(i), elmd(j), ke)
!       jac(nsf+contmp(map(elmd(j))))=jac(nsf+contmp(map(elmd(j))))
!      &      +ke*x_elm**p
!             acol(nsf+contmp(map(elmd(j))))  = map(elmd(j))
!             jac(nsf+27+nodes(map(elmd(i))))  =  
!      &      jac(nsf+27+nodes(map(elmd(i))))  
!      &      + ke*xprimal(map(elmd(j)))
!      &         *p*x_elm**(p-1.d0)
!             enddo
!             acol(nsf+27+nodes(map(elmd(i)))) = nsta + ep
!             nodes(map(elmd(i))) = nodes(map(elmd(i))) + 1
!         elseif(rowcounter(map(elmd(i))).eq.46) then
!          do j = 1,42
!             contmp(connect(map(elmd(i)), j)) = j
!             enddo
!             do j = 1,16
!             call stiffroutine(elmd(i), elmd(j), ke)
!       jac(nsf+contmp(map(elmd(j))))=jac(nsf+contmp(map(elmd(j))))
!      &       +ke*x_elm**p
!             acol(nsf+contmp(map(elmd(j))))  = map(elmd(j))
!             jac(nsf+43+nodes(map(elmd(i))))  =  
!      &      jac(nsf+43+nodes(map(elmd(i))))  
!      &      + ke*xprimal(map(elmd(j)))
!      &       *p*x_elm**(p-1.d0)
!             enddo
!             acol(nsf+43+nodes(map(elmd(i)))) = nsta + ep
!             nodes(map(elmd(i))) = nodes(map(elmd(i))) + 1
!         endif
!         enddo
! 
! !       con(nsta+1) = con(nsta+1) + xprimal(nsta+ep)
! !       jac(nnz - 4*bs - nmat + ep) = 1.d0
! !       acol(nnz - 4*bs - nmat + ep) = nsta+ep
! 
!         ! BC
!         if(ground.eq.1) then
! 
!             if(ey.eq.mult*5) then
!             xl(map(7)) = 0.d0
!             xu(map(7)) = 0.d0
!             xl(map(5)) = 0.0d0
!             xu(map(5)) = 0.0d0
!             xl(map(6)) = 0.0d0
!             xu(map(6)) = 0.0d0
!             xl(map(8)) = 0.d0
!             xu(map(8)) = 0.d0
!             xl(map(13)) = 0.0d0
!             xl(map(14)) = 0.0d0
!             xu(map(13)) = 0.0d0
!             xu(map(14)) = 0.0d0
!             con(map(7)) = 0.d0
!             con(map(5)) = 0.0d0
!             con(map(6)) = 0.0d0
!             con(map(8)) = 0.d0
!             con(map(13)) = 0.0d0
!             con(map(14)) = 0.0d0
!             endif
!         
!             bmap(1) = 7
!             bmap(2) = 8
!             bmap(3) = 5
!             bmap(4) = 6
!             bmap(5) = 13
!             bmap(6) = 14
!             if(ey.eq.mult*5) then
!             do i = 1,6
!             do j = 1,aptr(map(bmap(i))+1)-aptr(map(bmap(i))) 
!             jac(aptr(map(bmap(i)))-1+j) = 0.d0
!             enddo
!             enddo
!             endif
!     
!         elseif(ground.eq.2) then
! 
!             if(ex.eq.1) then
!             xl(map(7)) = 0.d0
!             xu(map(7)) = 0.d0
!             xl(map(15)) = 0.d0
!             xu(map(15)) = 0.d0
!             xl(map(1)) = 0.d0
!             xu(map(1)) = 0.d0
!             con(map(7)) = 0.0d0
!             con(map(15)) = 0.0d0
!             con(map(1)) = 0.0d0
!             endif
!         if(ey.le.mult/2.and.ex.eq.nex) then
!         xl(map(4)) = 0.d0
!         xu(map(4)) = 0.d0
!         con(map(4)) = 0.d0
!         xl(map(12)) = 0.d0
!         xu(map(12)) = 0.d0
!         con(map(12)) = 0.d0
!         endif
!             if(ex.eq.1) then
!             do j = 1,rowcounter(map(1))
!             jac(aptr(map(1))-1+j) = 0.d0
!             enddo
!             do j = 1,rowcounter(map(7))  !7
!             jac(aptr(map(7))-1+j) = 0.d0
!             enddo
!             do j = 1,rowcounter(map(15)) !15
!             jac(aptr(map(15))-1+j) = 0.d0
!             enddo
!             endif
!             if(ey.le.mult/2.and.ex.eq.nex) then
!             do j = 1,rowcounter(map(4))
!             jac(aptr(map(4))-1+j) = 0.d0
!             enddo
!             do j = 1,rowcounter(map(12))
!             jac(aptr(map(12))-1+j) = 0.d0
!             enddo
!             endif
!         elseif(ground.eq.3) then
! 
!         if(ey.eq.35*mult) then
!         xl(map(7)) = 0.d0
!         xu(map(7)) = 0.d0
!         xl(map(5)) = 0.0d0
!         xu(map(5)) = 0.0d0
!         xl(map(6)) = 0.0d0
!         xu(map(6)) = 0.0d0
!         xl(map(8)) = 0.d0
!         xu(map(8)) = 0.d0
!         xl(map(13)) = 0.0d0
!         xl(map(14)) = 0.0d0
!         xu(map(13)) = 0.0d0
!         xu(map(14)) = 0.0d0
!         endif
!         if(ey.eq.35*mult) then
!         con(map(7)) = 0.d0
!         con(map(5)) = 0.0d0
!         con(map(6)) = 0.0d0
!         con(map(8)) = 0.d0
!         con(map(13)) = 0.0d0
!         con(map(14)) = 0.0d0
!         endif
!         bmap(1) = 7
!         bmap(2) = 8
!         bmap(3) = 5
!         bmap(4) = 6
!         bmap(5) = 13
!         bmap(6) = 14
!         if(ey.eq.35*mult) then
!         do i = 1,6
!         do j = 1,aptr(map(bmap(i))+1)-aptr(map(bmap(i))) 
!         jac(aptr(map(bmap(i)))-1+j) = 0.d0
!         enddo
!         enddo
!         endif
!         endif
! 
!         call stressroutine_11(xprimal,xdual,
!      &  map,nu,sy,l,nsta,nmat,nnz,ep,
!      &  con,jac,hes,acol,stres,eps,bs)
!         enddo
!         enddo
! 
!         ! Slope constraints
!         bsc = 1
!         bscc = 1
!         do j = 1, ney
!         do i = 1, nex-1
!         bscon(bsc,1) = bscc
!         bscon(bsc,2) = bscc+1
!         bsc = bsc + 1
!         bscc = bscc + 1
!         enddo
!         bscc = bscc + 1
!         enddo
!         bsc = ney*(nex-1) + 1
!         bscc = 1
!         do j = 1, ney-1
!         do i = 1, nex
!         bscon(bsc,1) = bscc
!         bscon(bsc,2) = bscc + nex
!         bsc = bsc + 1
!         bscc = bscc + 1
!         enddo
!         enddo
! 
!         do i = 1, bs
!         con(nsta+nmat+i) =
!      &        xprimal(nsta+bscon(i,1))-xprimal(nsta+bscon(i,2))
!      &        -filt/DBLE(nex)
!         enddo
!         do i = 1, bs
!         con(nsta+nmat+bs+i) =
!      &        xprimal(nsta+bscon(i,2))-xprimal(nsta+bscon(i,1))
!      &        -filt/DBLE(nex)
!         enddo
! 
!         soc1 = 1
!         soc2 = 2
!         do i = 1, bs
!         jac(nnz-4*bs+soc1) = 1.0d0
!         jac(nnz-4*bs+soc2) = -1.0d0
!         acol(nnz-4*bs+soc1) = nsta+bscon(i,1)
!         acol(nnz-4*bs+soc2) = nsta+bscon(i,2)
!         soc1 = soc1 + 2
!         soc2 = soc2 + 2
!         enddo
! !
!         soc1 = 1
!         soc2 = 2
!         do i = 1, bs
!         jac(nnz-2*bs+soc1) = -1.0d0
!         jac(nnz-2*bs+soc2) = 1.0d0
!         acol(nnz-2*bs+soc1) = nsta+bscon(i,1)
!         acol(nnz-2*bs+soc2) = nsta+bscon(i,2)
!         soc1 = soc1 + 2
!         soc2 = soc2 + 2
!         enddo
! 
!         viol = 0.0d0
!         do i = 1,ndual
!         if(eqn(i)) then
!         viol = max(viol,abs(con(i)))
!         else
!         viol = max(viol,max(con(i),0.0d0))
!         endif
!         enddo
! 
!         return
!         end subroutine
        subroutine soil(nprimal,ndual,nprimalc,ndualc,xprimal,xprimalc,
     &  xdual,xdualc,bc,dx)
        
        implicit none
        integer             i,k
        integer             nprimal,ndual
        integer             nprimalc,ndualc
        double precision    xprimal(nprimal),xprimalc(nprimal)
        double precision    xdual(ndual),xdualc(ndual)
        integer             bc(nprimal),dx(nprimal)

        k=0
        do i = 1,ndual
        if(bc(i).eq.0) then
        k=k+1
        xdual(i)=xdualc(k)
        else
        xdual(i)=0.0d0
        endif
        enddo

        end subroutine
        subroutine stiffroutine(i,j, k)
        integer i, j
        double precision h,k,E,nu,ke(16,16)
        
        h = 1.0d0
        E = 1.00d2
        nu = 0.3d0
        
        k = h*E/(1.0d0-nu**2.0d0)
        ke(1,1) = -13.d0/45.d0*nu+13.d0/15.d0
        ke(2,1) = 17.d0/72.d0*nu + 17.d0/72.d0
        ke(3,1) = -17.d0/180.d0*nu+73.d0/180.d0
        ke(4,1) = -1.d0/8.d0*nu+1.d0/24.d0
        ke(5,1) = -23.d0/180.d0*nu+23.d0/60.d0
        ke(6,1) = 7.d0/72.d0*nu+7.d0/72.d0
        ke(7,1) = -7.d0/45.d0*nu+31.d0/90.d0
        ke(8,1) = 1.d0/8.d0*nu-1.d0/24.d0
        ke(9,1) = -1.d0/30.d0*nu-77.d0/90.d0
        ke(10,1) = 7.d0/18.d0*nu-5.d0/18.d0
        ke(11,1) = 2.d0/9.d0*nu-13.d0/45.d0
        ke(12,1) = -1.d0/18.d0*nu-1.d0/18.d0
        ke(13,1) = 1.d0/30.d0*nu-43.d0/90.d0
        ke(14,1) = -1.d0/18.d0*nu-1.d0/18.d0
        ke(15,1) = 4.d0/9.d0*nu-17.d0/45.d0
        ke(16,1) = -11.d0/18.d0*nu+1.d0/18.d0
    
        ke(1,2) =  17.d0/72.d0*nu+17.d0/72.d0
        ke(2,2) = -13.d0/45.d0*nu+13.d0/15.d0
        ke(3,2) = 1.d0/8.d0*nu-1.d0/24.d0
        ke(4,2) = -7.d0/45.d0*nu+31.d0/90.d0
        ke(5,2) = 7.d0/72.d0*nu+7.d0/72.d0
        ke(6,2) = -23.d0/180.d0*nu+23.d0/60.d0
        ke(7,2) = -1.d0/8.d0*nu+1.d0/24.d0
        ke(8,2) = -17.d0/180.d0*nu+73.d0/180.d0
        ke(9,2) = -11.d0/18.d0*nu+1.d0/18.d0
        ke(10,2) =4.d0/9.d0*nu-17.d0/45.d0
        ke(11,2) =-1.d0/18.d0*nu-1.d0/18.d0
        ke(12,2) =1.d0/30.d0*nu-43.d0/90.d0
        ke(13,2) =-1.d0/18.d0*nu-1.d0/18.d0
        ke(14,2) =2.d0/9.d0*nu-13.d0/45.d0
        ke(15,2) =7.d0/18.d0*nu-5.d0/18.d0
        ke(16,2) =-1.d0/30.d0*nu-77.d0/90.d0
    
        ke(1,3) =-17.d0/180.d0*nu+73.d0/180.d0 
        ke(2,3) =1.d0/8.d0*nu-1.d0/24.d0
        ke(3,3) =-13.d0/45.d0*nu+13.d0/15.d0
        ke(4,3) =-17.d0/72.d0*nu-17.d0/72.d0
        ke(5,3) =-7.d0/45.d0*nu+31.d0/90.d0
        ke(6,3) =-1.d0/8.d0*nu+1.d0/24.d0
        ke(7,3) = -23.d0/180.d0*nu+23.d0/60.d0 
        ke(8,3) = -7.d0/72.d0*nu-7.d0/72.d0
        ke(9,3) = -1.d0/30.d0*nu-77.d0/90.d0
        ke(10,3) = -7.d0/18.d0*nu+5.d0/18.d0
        ke(11,3) = 4.d0/9.d0*nu-17.d0/45.d0
        ke(12,3) = 11.d0/18.d0*nu-1.d0/18.d0
        ke(13,3) = 1.d0/30.d0*nu-43.d0/90.d0
        ke(14,3) = 1.d0/18.d0*nu+1.d0/18.d0
        ke(15,3) = 2.d0/9.d0*nu-13.d0/45.d0
        ke(16,3) = 1.d0/18.d0*nu+1.d0/18.d0
    
        ke(1,4) = -1.d0/8.d0*nu+1.d0/24.d0
        ke(2,4) = -7.d0/45.d0*nu+31.d0/90.d0
        ke(3,4) = -17.d0/72.d0*nu-17.d0/72.d0
        ke(4,4) = -13.d0/45.d0*nu+13.d0/15.d0
        ke(5,4) = 1.d0/8.d0*nu-1.d0/24.d0
        ke(6,4) = -17.d0/180.d0*nu+73.d0/180.d0
        ke(7,4) = -7.d0/72.d0*nu-7.d0/72.d0
        ke(8,4) = -23.d0/180.d0*nu+23.d0/60.d0
        ke(9,4) = 11.d0/18.d0*nu-1.d0/18.d0
        ke(10,4) = 4.d0/9.d0*nu-17.d0/45.d0
        ke(11,4) = -7.d0/18.d0*nu+5.d0/18.d0
        ke(12,4) = -1.d0/30.d0*nu-77.d0/90.d0
        ke(13,4) = 1.d0/18.d0*nu+1.d0/18.d0
        ke(14,4) = 2.d0/9.d0*nu-13.d0/45.d0
        ke(15,4) = 1.d0/18.d0*nu+1.d0/18.d0
        ke(16,4) = 1.d0/30.d0*nu-43.d0/90.d0
    
        ke(1,5) = -23.d0/180.d0*nu+23.d0/60.d0
        ke(2,5) = 7.d0/72.d0*nu+7.d0/72.d0
        ke(3,5) = -7.d0/45.d0*nu+31.d0/90.d0
        ke(4,5) = 1.d0/8.d0*nu-1.d0/24.d0
        ke(5,5) = -13.d0/45.d0*nu+13.d0/15.d0
        ke(6,5) = 17.d0/72.d0*nu+17.d0/72.d0
        ke(7,5) = -17.d0/180.d0*nu+73.d0/180.d0
        ke(8,5) = -1.d0/8.d0*nu+1.d0/24.d0
        ke(9,5) = 1.d0/30.d0*nu-43.d0/90.d0
        ke(10,5) = -1.d0/18.d0*nu -1.d0/18.d0
        ke(11,5) = 4.d0/9.d0*nu-17.d0/45.d0
        ke(12,5) = -11.d0/18.d0*nu+1.d0/18.d0
        ke(13,5) = -1.d0/30.d0*nu-77.d0/90.d0
        ke(14,5) = 7.d0/18.d0*nu -5.d0/18.d0
        ke(15,5) = 2.d0/9.d0*nu - 13.d0/45.d0
        ke(16,5) = -1.d0/18.d0*nu-1.d0/18.d0
    
        ke(1,6) = 7.d0/72.d0*nu+7.d0/72.d0 
        ke(2,6) = -23.d0/180.d0*nu+23.d0/60.d0
        ke(3,6) =  -1.d0/8.d0*nu +1.d0/24.d0
        ke(4,6) = -17.d0/180.d0*nu +73.d0/180.d0
        ke(5,6) = 17.d0/72.d0*nu+17.d0/72.d0
        ke(6,6) = -13.d0/45.d0*nu+13.d0/15.d0
        ke(7,6) = 1.d0/8.d0*nu-1.d0/24.d0
        ke(8,6) = -7.d0/45.d0*nu+31.d0/90.d0
        ke(9,6) = -1.d0/18.d0*nu-1.d0/18.d0
        ke(10,6) = 2.d0/9.d0*nu - 13.d0/45.d0
        ke(11,6) = 7.d0/18.d0*nu - 5.d0/18.d0
        ke(12,6) = -1.d0/30.d0*nu - 77.d0/90.d0
        ke(13,6) = -11.d0/18.d0*nu+1.d0/18.d0
        ke(14,6) = 4.d0/9.d0*nu -17.d0/45.d0
        ke(15,6) = -1.d0/18.d0*nu - 1.d0/18.d0
        ke(16,6) =  1.d0/30.d0*nu-43.d0/90.d0
    
        ke(1,7) =  -7.d0/45.d0*nu + 31.d0/90.d0
        ke(2,7) =  -1.d0/8.d0*nu+1.d0/24.d0
        ke(3,7) = -23.d0/180.d0*nu+23.d0/60.d0
        ke(4,7) = -7.d0/72.d0*nu-7.d0/72.d0
        ke(5,7) = -17.d0/180.d0*nu+73.d0/180.d0
        ke(6,7) = 1.d0/8.d0*nu-1.d0/24.d0
        ke(7,7) = -13.d0/45.d0*nu+13.d0/15.d0
        ke(8,7) = -17.d0/72.d0*nu-17.d0/72.d0
        ke(9,7) = 1.d0/30.d0*nu-43.d0/90.d0
        ke(10,7) =1.d0/18.d0*nu+1.d0/18.d0
        ke(11,7) =2.d0/9.d0*nu-13.d0/45.d0
        ke(12,7) =1.d0/18.d0*nu+1.d0/18.d0
        ke(13,7) =-1.d0/30.d0*nu-77.d0/90.d0
        ke(14,7) =-7.d0/18.d0*nu+5.d0/18.d0
        ke(15,7) =4.d0/9.d0*nu - 17.d0/45.d0
        ke(16,7) = 11.d0/18.d0*nu-1.d0/18.d0
    
        ke(1,8) = 1.d0/8.d0*nu-1.d0/24.d0
        ke(2,8) = -17.d0/180.d0*nu+73.d0/180.d0
        ke(3,8) = -7.d0/72.d0*nu-7.d0/72.d0
        ke(4,8) = -23.d0/180.d0*nu+23.d0/60.d0
        ke(5,8) = -1.d0/8.d0*nu+1.d0/24.d0
        ke(6,8) = -7.d0/45.d0*nu+31.d0/90.d0
        ke(7,8) = -17.d0/72.d0*nu-17.d0/72.d0
        ke(8,8) = -13.d0/45.d0*nu+13.d0/15.d0
        ke(9,8) = 1.d0/18.d0*nu+1.d0/18.d0
        ke(10,8) = 2.d0/9.d0*nu-13.d0/45.d0
        ke(11,8) = 1.d0/18.d0*nu+1.d0/18.d0
        ke(12,8) = 1.d0/30.d0*nu-43.d0/90.d0
        ke(13,8) = 11.d0/18.d0*nu-1.d0/18.d0
        ke(14,8) = 4.d0/9.d0*nu-17.d0/45.d0
        ke(15,8) = -7.d0/18.d0*nu+5.d0/18.d0
        ke(16,8) = -1.d0/30.d0*nu-77.d0/90.d0
    
        ke(1,9) = -1.d0/30.d0*nu-77.d0/90.d0 
        ke(2,9) = -11.d0/18.d0*nu+1.d0/18.d0
        ke(3,9) = -1.d0/30.d0*nu-77.d0/90.d0
        ke(4,9) = 11.d0/18.d0*nu-1.d0/18.d0
        ke(5,9) = 1.d0/30.d0*nu-43.d0/90.d0
        ke(6,9) = -1.d0/18.d0*nu-1.d0/18.d0
        ke(7,9) = 1.d0/30.d0*nu-43.d0/90.d0
        ke(8,9) = 1.d0/18.d0*nu+1.d0/18.d0
        ke(9,9) = -4.d0/15.d0*nu+92.d0/45.d0
        ke(10,9) = 0.d0
        ke(11,9) = 0.d0
        ke(12,9) = -2.d0/9.d0*nu - 2.d0/9.d0
        ke(13,9) = 4.d0/15.d0*nu + 28.d0/45.d0
        ke(14,9) = 0.d0
        ke(15,9) = 0.d0
        ke(16,9) = 2.d0/9.d0*nu + 2.d0/9.d0
    
        ke(1,10) = 7.d0/18.d0*nu - 5.d0/18.d0
        ke(2,10) = 4.d0/9.d0*nu -17.d0/45.d0
        ke(3,10) = -7.d0/18.d0*nu+5.d0/18.d0
        ke(4,10) = 4.d0/9.d0*nu-17.d0/45.d0
        ke(5,10) = -1.d0/18.d0*nu-1.d0/18.d0
        ke(6,10) = 2.d0/9.d0*nu-13.d0/45.d0
        ke(7,10) = 1.d0/18.d0*nu+1.d0/18.d0
        ke(8,10) = 2.d0/9.d0*nu-13.d0/45.d0
        ke(9,10) = 0.d0
        ke(10,10) =  -8.d0/9.d0*nu+64.d0/45.d0
        ke(11,10) = -2.d0/9.d0*nu - 2.d0/9.d0
        ke(12,10) = 0.d0
        ke(13,10) = 0.d0
        ke(14,10) = -4.d0/9.d0*nu - 4.d0/45.d0
        ke(15,10) = 2.d0/9.d0*nu + 2.d0/9.d0
        ke(16,10) = 0.d0
    
        ke(1,11) = 2.d0/9.d0*nu-13.d0/45.d0
        ke(2,11) = -1.d0/18.d0*nu-1.d0/18.d0
        ke(3,11) = 4.d0/9.d0*nu-17.d0/45.d0
        ke(4,11) = -7.d0/18.d0*nu+5.d0/18.d0
        ke(5,11) = 4.d0/9.d0*nu - 17.d0/45.d0
        ke(6,11) = 7.d0/18.d0*nu-5.d0/18.d0
        ke(7,11) = 2.d0/9.d0*nu - 13.d0/45.d0
        ke(8,11) = 1.d0/18.d0*nu+1.d0/18.d0
        ke(9,11) = 0.d0
        ke(10,11) =  -2.d0/9.d0*nu - 2.d0/9.d0
        ke(11,11) = -8.d0/9.d0*nu + 64.d0/45.d0
        ke(12,11) =  0.d0
        ke(13,11) =  0.d0
        ke(14,11) =  2.d0/9.d0*nu + 2.d0/9.d0
        ke(15,11) = -4.d0/9.d0*nu-4.d0/45.d0
        ke(16,11) = 0.d0
    
        ke(1,12) = -1.d0/18.d0*nu-1.d0/18.d0
        ke(2,12) = 1.d0/30.d0*nu-43.d0/90.d0
        ke(3,12) = 11.d0/18.d0*nu-1.d0/18.d0
        ke(4,12) = -1.d0/30.d0*nu-77.d0/90.d0
        ke(5,12) = -11.d0/18.d0*nu+1.d0/18.d0
        ke(6,12) = -1.d0/30.d0*nu-77.d0/90.d0
        ke(7,12) = 1.d0/18.d0*nu+1.d0/18.d0
        ke(8,12) = 1.d0/30.d0*nu-43.d0/90.d0
        ke(9,12) = -2.d0/9.d0*nu-2.d0/9.d0
        ke(10,12) = 0.d0
        ke(11,12) = 0.d0
        ke(12,12) = -4.d0/15.d0*nu+92.d0/45.d0
        ke(13,12) = 2.d0/9.d0*nu+2.d0/9.d0
        ke(14,12) = 0.d0
        ke(15,12) = 0.d0
        ke(16,12) = 4.d0/15.d0*nu+28.d0/45.d0
    
        ke(1,13) = 1.d0/30.d0*nu-43.d0/90.d0
        ke(2,13) = -1.d0/18.d0*nu-1.d0/18.d0
        ke(3,13) = 1.d0/30.d0*nu-43.d0/90.d0
        ke(4,13) = 1.d0/18.d0*nu+1.d0/18.d0
        ke(5,13) = -1.d0/30.d0*nu-77.d0/90.d0
        ke(6,13) = -11.d0/18.d0*nu+1.d0/18.d0
        ke(7,13) = -1.d0/30.d0*nu-77.d0/90.d0
        ke(8,13) = 11.d0/18.d0*nu-1.d0/18.d0
        ke(9,13) = 4.d0/15.d0*nu+28.d0/45.d0
        ke(10,13) = 0.d0
        ke(11,13) = 0.d0
        ke(12,13) = 2.d0/9.d0*nu+2.d0/9.d0
        ke(13,13) = -4.d0/15.d0*nu+92.d0/45.d0
        ke(14,13) = 0.d0
        ke(15,13) = 0.d0
        ke(16,13) = -2.d0/9.d0*nu-2.d0/9.d0
    
        ke(1,14) = -1.d0/18.d0*nu-1.d0/18.d0
        ke(2,14) = 2.d0/9.d0*nu-13.d0/45.d0
        ke(3,14) = 1.d0/18.d0*nu+1.d0/18.d0
        ke(4,14) = 2.d0/9.d0*nu-13.d0/45.d0
        ke(5,14) = 7.d0/18.d0*nu-5.d0/18.d0
        ke(6,14) = 4.d0/9.d0*nu-17.d0/45.d0
        ke(7,14) = -7.d0/18.d0*nu+5.d0/18.d0
        ke(8,14) = 4.d0/9.d0*nu-17.d0/45.d0
        ke(9,14) = 0.d0
        ke(10,14) = -4.d0/9.d0*nu-4.d0/45.d0
        ke(11,14) = 2.d0/9.d0*nu+2.d0/9.d0
        ke(12,14) = 0.d0
        ke(13,14) = 0.d0
        ke(14,14) = -8.d0/9.d0*nu+64.d0/45.d0
        ke(15,14) = -2.d0/9.d0*nu-2.d0/9.d0
        ke(16,14) = 0.d0
    
        ke(1,15) = 4.d0/9.d0*nu-17.d0/45.d0
        ke(2,15) = 7.d0/18.d0*nu-5.d0/18.d0
        ke(3,15) = 2.d0/9.d0*nu-13.d0/45.d0
        ke(4,15) = 1.d0/18.d0*nu+1.d0/18.d0
        ke(5,15) = 2.d0/9.d0*nu-13.d0/45.d0
        ke(6,15) = -1.d0/18.d0*nu-1.d0/18.d0
        ke(7,15) = 4.d0/9.d0*nu-17.d0/45.d0
        ke(8,15) = -7.d0/18.d0*nu+5.d0/18.d0
        ke(9,15) = 0.d0
        ke(10,15) = 2.d0/9.d0*nu+2.d0/9.d0
        ke(11,15) = -4.d0/9.d0*nu-4.d0/45.d0
        ke(12,15) = 0.d0
        ke(13,15) = 0.d0
        ke(14,15) = -2.d0/9.d0*nu-2.d0/9.d0
        ke(15,15) = -8.d0/9.d0*nu+64.d0/45.d0
        ke(16,15) = 0.d0
    
        ke(1,16) = -11.d0/18.d0*nu+1.d0/18.d0
        ke(2,16) = -1.d0/30.d0*nu-77.d0/90.d0
        ke(3,16) = 1.d0/18.d0*nu+1.d0/18.d0
        ke(4,16) = 1.d0/30.d0*nu-43.d0/90.d0
        ke(5,16) = -1.d0/18.d0*nu-1.d0/18.d0
        ke(6,16) = 1.d0/30.d0*nu-43.d0/90.d0
        ke(7,16) = 11.d0/18.d0*nu-1.d0/18.d0
        ke(8,16) = -1.d0/30.d0*nu-77.d0/90.d0
        ke(9,16) = 2.d0/9.d0*nu+2.d0/9.d0
        ke(10,16) = 0.d0
        ke(11,16) = 0.d0
        ke(12,16) = 4.d0/15.d0*nu+28.d0/45.d0
        ke(13,16) = -2.d0/9.d0*nu-2.d0/9.d0
        ke(14,16) = 0.d0
        ke(15,16) = 0.d0
        ke(16,16) = -4.d0/15.d0*nu+92.d0/45.d0
        
        k = k*ke(i,j)
       
        return
        end subroutine
        subroutine stresprint_k(wm, x, n)
        implicit         none
        integer n, i, j, te, wm
        integer nex, ney, prob
        double precision x(n), l, p
        character (len=100) :: temp


        nex = 15*wm
        ney = 5*wm

        j = 1
        open(unit=2, file="stres.dat", action="write",status="replace")
        do i = 1,ney
        write(2, '(1600E20.7)') (SNGL(x((i-1)*nex + j)),j =1,nex)
        enddo
        close(2)

        write(temp, 999) nex, ney

        CALL SYSTEM(temp)
        j = 1
        open(unit=2, file="stres.dat", action="write",status="replace")
        do i = 1,ney*nex
        write(2, '(1600E20.7)') x(i)
        enddo
        close(2)

        return;

!999     format ('gnuplot -e "xi =', I0,';yi =',I0, 
!     &  '" topoprint.p 2>&1')
999     format ('gnuplot -e "xi =', I0,';yi =',I0, 
     &  '" stresprint.p 2>&1 | grep -v font' )

        end subroutine
        subroutine ltresprint_k(wm, x, n)
        implicit         none
        integer n, i, j, te, wm
        integer nex, ney, prob
        double precision x(n), l, p
        character (len=100) :: temp


        nex = 25*wm
        ney = 25*wm

        j = 1
        open(unit=2, file="stres.dat", action="write",status="replace")
        do i = 1,ney
        write(2, '(1600E20.7)') (SNGL(x((i-1)*nex + j)),j =1,nex)
        enddo
        close(2)

        write(temp, 999) nex, ney

        CALL SYSTEM(temp)
        j = 1
        open(unit=2, file="stres.dat", action="write",status="replace")
        do i = 1,ney*nex
        write(2, '(1600E20.7)') x(i)
        enddo
        close(2)

        return;

!999     format ('gnuplot -e "xi =', I0,';yi =',I0, 
!     &  '" topoprint.p 2>&1')
999     format ('gnuplot -e "xi =', I0,';yi =',I0, 
     &  '" stresprint.p 2>&1 | grep -v font' )

        end subroutine
        subroutine stressroutine_1(xprimal,xdual,
     &  map,nu,sy,l,nsta,nmat,nnz,ep,
     &  con,jac,hes,acol,stres)
        implicit none
        integer nsta, ep, acol(*), map(16), nmat, nnz
        double precision xprimal(*), xdual(*)
        double precision u5,v5,u6,v6,u7,v7,u8,v8,nu,E,sy,l
        double precision con(*), jac(*), hes(*)
        double precision sig(3), stres(*)

        E = 1d2
        u5=xprimal(map(9))
        v5=xprimal(map(10))
        u6=xprimal(map(11))
        v6=xprimal(map(12))
        u7=xprimal(map(13))
        v7=xprimal(map(14))
        u8=xprimal(map(15))
        v8=xprimal(map(16))

        sig(1) = nu*(v5-v7)-u6+u8
        sig(2) = nu*(u6-u8)-v5+v7
        sig(3) = u5-u7-v6+v8
        sig(1) = E/(nu**2 - 1.d0)*sig(1)/l
        sig(2) = E/(-1.0d0*(nu**2 - 1.d0))*sig(2)/l
        sig(3) = E/(-2.d0*(nu + 1.d0))*sig(3)/l

        stres(ep)= 
     &  sqrt(sig(1)**2.d0-sig(1)*sig(2)+sig(2)**2.d0+
     &  3.d0*sig(3)**2.d0)

        con(nsta+ep) =xprimal(nsta+ep)*
     &  (sqrt(sig(1)**2.d0-sig(1)*sig(2)+sig(2)**2.d0+
     &  3.d0*sig(3)**2.d0)/(sy)
     &  -1.0d0)
      
        !to u5
        jac(nnz-9*nmat+(ep-1)*9+1) = 
     &  3.0d0*(u5-u7-v6+v8)*E**2.0d0
     &  /(2.0d0*(nu+1.0d0)**2.0d0*l**2.0d0)/sy
     &  *xprimal(nsta+ep)
     &  /stres(ep)/2.0d0

        !to v5
        jac(nnz-9*nmat+(ep-1)*9+2)=(
     &  2.0d0*(nu*v5-nu*v7-u6+u8)*E**2.0d0*nu/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0)
     &  +(u6*nu-u8*nu-v5+v7)*E**2.0d0*nu/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0)
     &  -(nu*v5-nu*v7-u6+u8)*E**2.0d0/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0)
     &  -2.0d0*(nu*u6-nu*u8-v5+v7)*E**2.0d0/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0))/sy
     &  *xprimal(nsta+ep)
     &  /stres(ep)/2.0d0
   
        !to u8
        jac(nnz-9*nmat+(ep-1)*9+3) = (
     &  -1.0d0*(nu*v5-nu*v7-u6+u8)*E**2.0d0*nu/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0)
     &  -2.0d0*(u6*nu-u8*nu-v5+v7)*E**2.0d0*nu/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0)
     &  +2.0d0*(nu*v5-nu*v7-u6+u8)*E**2.0d0/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0)
     &  +(nu*u6-nu*u8-v5+v7)*E**2.0d0/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0))/sy
     &  *xprimal(nsta+ep)
     &  /stres(ep)/2.0d0

        !to v8
        jac(nnz-9*nmat+(ep-1)*9+4) = 
     &  3.0d0*(u5-u7-v6+v8)*E**2.0d0
     &  /(2.0d0*(nu+1.0d0)**2.0d0*l**2.0d0)/sy
     &  *xprimal(nsta+ep)
     &  /stres(ep)/2.0d0
     
        !to u6
        jac(nnz-9*nmat+(ep-1)*9+5) = (
     &  (nu*v5-nu*v7-u6+u8)*E**2.0d0*nu/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0)
     &  +2.0d0*(u6*nu-u8*nu-v5+v7)*E**2.0d0*nu/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0)
     &  -2.0d0*(nu*v5-nu*v7-u6+u8)*E**2.0d0/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0)
     &  -(nu*u6-nu*u8-v5+v7)*E**2.0d0/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0))/sy
     &  *xprimal(nsta+ep)
     &  /stres(ep)/2.0d0
     
        !to v6
        jac(nnz-9*nmat+(ep-1)*9+6) = 
     &  -3.0d0*(u5-u7-v6+v8)*E**2.0d0
     &  /(2.0d0*(nu+1.0d0)**2.0d0*l**2.0d0)/sy
     &  *xprimal(nsta+ep)
     &  /stres(ep)/2.0d0

        !to u7
        jac(nnz-9*nmat+(ep-1)*9+7) = 
     &  -3.0d0*(u5-u7-v6+v8)*E**2.0d0
     &  /(2.0d0*(nu+1.0d0)**2.0d0*l**2.0d0)/sy
     &  *xprimal(nsta+ep)
     &  /stres(ep)/2.0d0

        !to v7
        jac(nnz-9*nmat+(ep-1)*9+8) = (
     &  -2.0d0*(nu*v5-nu*v7-u6+u8)*E**2.0d0*nu/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0)
     &  -(u6*nu-u8*nu-v5+v7)*E**2.0d0*nu/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0)
     &  +(nu*v5-nu*v7-u6+u8)*E**2.0d0/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0)
     &  +2.0d0*(nu*u6-nu*u8-v5+v7)*E**2.0d0/
     &  ((nu**2.0d0-1.0d0)**2.0d0*l**2.0d0))/sy
     &  *xprimal(nsta+ep)
     &  /stres(ep)/2.0d0
        
        !to x
        jac(nnz-9*nmat+(ep-1)*9+9) =
!    &  ((sig(1)**2.d0-sig(1)*sig(2)+sig(2)**2.d0+
!    &  3.d0*sig(3)**2.d0)/sy**2.0d0
!    &  -1.0d0)
     &  stres(ep)/sy-1.0d0
!       else
!       jac(nnz-9*nmat+(ep-1)*9+9) = 0.0d0
!       endif

!       if(xprimal(nsta+ep).gt.0.0d0) then
!    &  (sig(1)**2.d0-sig(1)*sig(2)+sig(2)**2.d0+
!    &  3.d0*sig(3)**2.d0)/l**2.0d0)
!       else
!       stres(ep) = 0.0d0
!       endif

        acol(nnz-9*nmat+(ep-1)*9+1) = map(9)  !u5
        acol(nnz-9*nmat+(ep-1)*9+2) = map(10) !v5
        acol(nnz-9*nmat+(ep-1)*9+3) = map(15) !u8
        acol(nnz-9*nmat+(ep-1)*9+4) = map(16) !v8
        acol(nnz-9*nmat+(ep-1)*9+5) = map(11) !u6
        acol(nnz-9*nmat+(ep-1)*9+6) = map(12) !v6
        acol(nnz-9*nmat+(ep-1)*9+7) = map(13) !u7
        acol(nnz-9*nmat+(ep-1)*9+8) = map(14) !v7
        acol(nnz-9*nmat+(ep-1)*9+9) = nsta+ep !x
    
!       goto 1234
        !-u5
        hes(map(9)) = hes(map(9))+3.0d0*E**2.0d0
     &  /(2.0d0*(nu+1.0d0)**2.0d0*l**2.0d0)
     &  /sy**2.0d0
     &  *xprimal(nsta+ep)
     &  *xdual(nsta+ep)
        !-v5
        hes(map(10)) = hes(map(10))+2.0d0*
     &  (nu**2.0d0-nu+1.0d0)*E**2.0d0
     &  /(nu**2.0d0-1.0d0)**2.0d0/l**2.0d0
     &  /sy**2.0d0
     &  *xprimal(nsta+ep)
     &  *xdual(nsta+ep)
        !-u8
        hes(map(15)) = hes(map(15))+2.0d0*
     &  (nu**2.0d0-nu+1.0d0)*E**2.0d0
     &  /(nu**2.0d0-1.0d0)**2.0d0/l**2.0d0
     &  /sy**2.0d0
     &  *xprimal(nsta+ep)
     &  *xdual(nsta+ep)
        !-v8
        hes(map(16)) = hes(map(16))+3.0d0*
     &  E**2.0d0
     &  /(2.0d0*(nu+1.0d0)**2.0d0*l**2.0d0)
     &  /sy**2.0d0
     &  *xprimal(nsta+ep)
     &  *xdual(nsta+ep)
        !-u6
        hes(map(11)) = hes(map(11))+2.0d0*
     &  (nu**2.0d0-nu+1.0d0)*E**2.0d0
     &  /(nu**2.0d0-1.0d0)**2.0d0/l**2.0d0
     &  /sy**2.0d0
     &  *xprimal(nsta+ep)
     &  *xdual(nsta+ep)
        !-v6
        hes(map(12)) = hes(map(12)) + 3.0d0*
     &  E**2.0d0
     &  /(2.0d0*(nu+1.0d0)**2.0d0*l**2.0d0)
     &  /sy**2.0d0
     &  *xprimal(nsta+ep)
     &  *xdual(nsta+ep)
        !-u7
        hes(map(13)) = hes(map(13)) + 3.0d0*
     &  E**2.0d0
     &  /(2.0d0*(nu+1.0d0)**2.0d0*l**2.0d0)
     &  /sy**2.0d0
     &  *xprimal(nsta+ep)
     &  *xdual(nsta+ep)
        !-v7
        hes(map(14))= hes(map(14)) + 2.0d0*
     &  (nu**2.0d0-nu+1.0d0)*E**2.0d0
     &  /(nu**2.0d0-1.0d0)**2.0d0/l**2.0d0
     &  /sy**2.0d0
     &  *xprimal(nsta+ep)
     &  *xdual(nsta+ep)

1234    continue
        return
        end subroutine
        subroutine stressroutine_11(xprimal,xdual,
     &  map,nu,sy,l,nsta,nmat,nnz,ep,
     &  con,jac,hes,acol,stres,bs)
        implicit none
        integer nsta, ep, acol(*), map(16), nmat, nnz, bs
        integer i, scal
        double precision xprimal(*), xdual(*)
        double precision u5,v5,u6,v6,u7,v7,u8,v8,nu,E,sy,l
        double precision con(*), jac(*), hes(*)
        double precision sig(3), stres(*), Eb,sigma

        E = 1d2

        u5=xprimal(map(9))
        v5=xprimal(map(10))
        u6=xprimal(map(11))
        v6=xprimal(map(12))
        u7=xprimal(map(13))
        v7=xprimal(map(14))
        u8=xprimal(map(15))
        v8=xprimal(map(16))

        scal=1
        sig(1) = nu*(v5-v7)-u6+u8
        sig(2) = nu*(u8-u6)+v5-v7
        sig(3) = u7-u5+v6-v8
        if(scal.eq.0) then
        Eb=E/(1.0d0 - nu**2.0d0)/l
        endif
        if(scal.eq.1) then
        Eb=E/(1.0d0 - nu**2.0d0)
        endif

        sigma=sig(1)**2.0d0-sig(1)*sig(2)+sig(2)**2.0d0
     &  +3.0d0/4.0d0*(nu-1.0d0)**2.0d0*sig(3)**2.0d0

        if(scal.eq.0) then
        stres(ep)= Eb*sqrt(sigma)!/l
        endif
        if(scal.eq.1) then
        stres(ep)= Eb*sqrt(sigma)/l
        endif
      
        if(scal.eq.0) then
        con(nsta+ep) =xprimal(nsta+ep)*(Eb*sqrt(sigma)/sy-1.0d0)!-1.0d0*l)
        endif
        if(scal.eq.1) then
        con(nsta+ep) =xprimal(nsta+ep)*(Eb*sqrt(sigma)/sy-1.0d0*l)
        endif

        !to u5
        jac(nnz-4*bs-9*nmat+(ep-1)*9+1) = 
     &  -3.0d0/4.0d0*Eb*(nu-1.0d0)**2.0d0*sig(3)/sqrt(sigma)
     &  /sy
     &  *xprimal(nsta+ep)

        !to v5
        jac(nnz-4*bs-9*nmat+(ep-1)*9+2)=
     &  Eb/2.0d0/sqrt(sigma)*((2.0d0*sig(1)-sig(2))*nu
     &  +2.0d0*sig(2)-sig(1))
     &  /sy
     &  *xprimal(nsta+ep)
   
        !to u8
        jac(nnz-4*bs-9*nmat+(ep-1)*9+3) = 
     &  Eb/2.0d0/sqrt(sigma)*((2*sig(2)-sig(1))*nu
     &  +2.0d0*sig(1)-sig(2))
     &  /sy
     &  *xprimal(nsta+ep)

        !to v8
        jac(nnz-4*bs-9*nmat+(ep-1)*9+4) = 
     &  -3.0d0/4.0d0*Eb*(nu-1.0d0)**2.0d0*sig(3)/sqrt(sigma)
     &  /sy
     &  *xprimal(nsta+ep)
     
        !to u6
        jac(nnz-4*bs-9*nmat+(ep-1)*9+5) = 
     &  -Eb/2.0d0/sqrt(sigma)*((2*sig(2)-sig(1))*nu
     &  +2.0d0*sig(1)-sig(2))
     &  /sy
     &  *xprimal(nsta+ep)
     
        !to v6
        jac(nnz-4*bs-9*nmat+(ep-1)*9+6) = 
     &  3.0d0/4.0d0*Eb*(nu-1.0d0)**2.0d0*sig(3)/sqrt(sigma)
     &  /sy
     &  *xprimal(nsta+ep)

        !to u7
        jac(nnz-4*bs-9*nmat+(ep-1)*9+7) = 
     &  3.0d0/4.0d0*Eb*(nu-1.0d0)**2.0d0*sig(3)/sqrt(sigma)
     &  /sy
     &  *xprimal(nsta+ep)

        !to v7
        jac(nnz-4*bs-9*nmat+(ep-1)*9+8) = 
     &  -Eb/2.0d0/sqrt(sigma)*((2.0d0*sig(1)-sig(2))*nu
     &  +2.0d0*sig(2)-sig(1))
     &  /sy
     &  *xprimal(nsta+ep)
        
        !to x
        if(scal.eq.0) then
        jac(nnz-4*bs-9*nmat+(ep-1)*9+9) =
     &  Eb*sqrt(sigma)/sy -1.0d0!*l
        endif
        if(scal.eq.1) then
        jac(nnz-4*bs-9*nmat+(ep-1)*9+9) =
     &  Eb*sqrt(sigma)/sy -1.0d0*l
        endif

        acol(nnz-4*bs-9*nmat+(ep-1)*9+1) = map(9)  !u5
        acol(nnz-4*bs-9*nmat+(ep-1)*9+2) = map(10) !v5
        acol(nnz-4*bs-9*nmat+(ep-1)*9+3) = map(15) !u8
        acol(nnz-4*bs-9*nmat+(ep-1)*9+4) = map(16) !v8
        acol(nnz-4*bs-9*nmat+(ep-1)*9+5) = map(11) !u6
        acol(nnz-4*bs-9*nmat+(ep-1)*9+6) = map(12) !v6
        acol(nnz-4*bs-9*nmat+(ep-1)*9+7) = map(13) !u7
        acol(nnz-4*bs-9*nmat+(ep-1)*9+8) = map(14) !v7
        acol(nnz-4*bs-9*nmat+(ep-1)*9+9) = nsta+ep !x
    
!       do this:
!       goto 1234
        !-u5
        hes(map(9)) = hes(map(9))
     &  -3.0d0/4.0d0*Eb*(nu-1.0d0)**2.0d0*
     &  (3.0d0/4.0d0*sig(3)**2.0d0/sigma**(3.0d0/2.0d0)
     &  *(nu-1.0d0)**2.0d0-1.0d0/sqrt(sigma))
     &  *xprimal(nsta+ep)    
     &  *xdual(nsta+ep)

        !-v5
        hes(map(10)) = hes(map(10))
     &  +Eb/2.0d0*(2.0d0*(nu**2.0d0-nu+1.0d0)/sqrt(sigma)
     &  -((2.0d0*sig(1)-sig(2))*nu
     &  +(2.0d0*sig(2)-sig(1)))**2.0d0
     &  /sigma**(3.0d0/2.0d0)/2.0d0)
     &  *xprimal(nsta+ep)    
     &  *xdual(nsta+ep)

        !-u8
        hes(map(15)) = hes(map(15))
     &  +Eb/2.0d0*(2.0d0*(nu**2.0d0-nu+1.0d0)/sqrt(sigma)
     &  -((2.0d0*sig(2)-sig(1))*nu
     &  +(2.0d0*sig(1)-sig(2)))**2.0d0
     &  /sigma**(3.0d0/2.0d0)/2.0d0)
     &  *xprimal(nsta+ep)    
     &  *xdual(nsta+ep)

        !-v8
        hes(map(16)) = hes(map(16))
     &  -3.0d0/4.0d0*Eb*(nu-1.0d0)**2.0d0*
     &  (3.0d0/4.0d0*sig(3)**2.0d0/sigma**(3.0d0/2.0d0)
     &  *(nu-1.0d0)**2.0d0-1.0d0/sqrt(sigma))
     &  *xprimal(nsta+ep)    
     &  *xdual(nsta+ep)

        !-u6
        hes(map(11)) = hes(map(11))
     &  -Eb/2.0d0*(2.0d0*(-nu**2.0d0+nu-1.0d0)/sqrt(sigma)
     &  +((2.0d0*sig(2)-sig(1))*nu
     &  +(2.0d0*sig(1)-sig(2)))**2.0d0
     &  /sigma**(3.0d0/2.0d0)/2.0d0)
     &  *xprimal(nsta+ep)    
     &  *xdual(nsta+ep)

        !-v6
        hes(map(12)) = hes(map(12))
     &  +3.0d0/4.0d0*Eb*(nu-1.0d0)**2.0d0*
     &  (-3.0d0/4.0d0*sig(3)**2.0d0/sigma**(3.0d0/2.0d0)
     &  *(nu-1.0d0)**2.0d0+1.0d0/sqrt(sigma))
     &  *xprimal(nsta+ep)    
     &  *xdual(nsta+ep)
        
        !-u7
        hes(map(13)) = hes(map(13))
     &  +3.0d0/4.0d0*Eb*(nu-1.0d0)**2.0d0*
     &  (-3.0d0/4.0d0*sig(3)**2.0d0/sigma**(3.0d0/2.0d0)
     &  *(nu-1.0d0)**2.0d0+1.0d0/sqrt(sigma))
     &  *xprimal(nsta+ep)    
     &  *xdual(nsta+ep)

        !-v7
        hes(map(14))= hes(map(14)) 
     &  +Eb/2.0d0*(2.0d0*(nu**2.0d0-nu+1.0d0)/sqrt(sigma)
     &  -((2.0d0*sig(1)-sig(2))*nu
     &  +(2.0d0*sig(2)-sig(1)))**2.0d0
     &  /sigma**(3.0d0/2.0d0)/2.0d0)
     &  *xprimal(nsta+ep)    
     &  *xdual(nsta+ep)

!       if(xprimal(nsta+ep).le.1d-3) then
!       stres(ep) = 0.0d0
!       endif

1234    continue
        return
        end subroutine
! SAOi: 
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Timer 
! SAOi: 


      real function seconds()
      implicit         none
      integer          icount, icount_rate, icount_max 
      real             time
      real             dummy(2), etime       ! single precision, i.e. real*4 !

      call system_clock(icount, icount_rate, icount_max)
      time = real(icount)/real(icount_rate)
      seconds = time
!
      return
      end function seconds
!-----------------------------------------------------------------------
        subroutine topoprint_k(wm, x, n, nex, ney)
        implicit         none
        integer n, i, j, te, wm
        integer nex, ney, tn, prob
        double precision x(n), l, p
        character (len=100) :: temp

        tn = (2*nex+1)*(ney+1)+ney*(nex+1)
        j = 1
        open(unit=2, file="topo.dat", action="write",status="replace")
        do i = 1,ney
        write(2, '(1600E20.7)') (SNGL(x(2*tn + (i-1)*nex + j)),j =1,nex)
        enddo
        close(2)

        write(temp, 999) nex, ney

        CALL SYSTEM(temp)

        return;

!999     format ('gnuplot -e "xi =', I0,';yi =',I0, 
!     &  '" topoprint.p 2>&1')
999     format ('gnuplot -e "xi =', I0,';yi =',I0, 
     &  '" topoprint_k.p 2>&1 | grep -v Warning: | grep -v font')

        end subroutine
        subroutine topoprint_star(wm, x, n, nex, ney)
        implicit         none
        integer n, i, j, te, wm
        integer nex, ney, tn, prob
        double precision x(n), l, p
        character (len=100) :: temp

        tn = (2*nex+1)*(ney+1)+ney*(nex+1)

        j = 1
        open(unit=2, file="topo.dat", action="write",status="replace")
        do i = 1,ney
        write(2, '(1600E20.7)') (SNGL(x(2*tn + (i-1)*nex + j)),j =1,nex)
        enddo
        close(2)

        write(temp, 999) nex, ney

        CALL SYSTEM(temp)

        return;

!999     format ('gnuplot -e "xi =', I0,';yi =',I0, 
!     &  '" topoprint.p 2>&1')
999     format ('gnuplot -e "xi =', I0,';yi =',I0, 
     &  '" topoprint_star.p 2>&1 | grep -v Warning:')

        end subroutine
        subroutine lopoprint_k(wm, x, n, nex, ney)
        implicit         none
        integer n, i, j, te, wm
        integer nex, ney, tn, prob
        double precision x(n), l, p
        character (len=100) :: temp

        tn = (2*nex+1)*(ney+1)+ney*(nex+1)

        j = 1
        open(unit=2, file="topo.dat", action="write",status="replace")
        do i = 1,ney
        write(2, '(1600E20.7)') (SNGL(x(2*tn + (i-1)*nex + j)),j =1,nex)
        enddo
        close(2)

        write(temp, 999) nex, ney

        CALL SYSTEM(temp)

        return;

!999     format ('gnuplot -e "xi =', I0,';yi =',I0, 
!     &  '" topoprint.p 2>&1')
999     format ('gnuplot -e "xi =', I0,';yi =',I0, 
     &  '" lopoprint_k.p 2>&1 | grep -v Warning:')

        end subroutine
        subroutine lopoprint_star(wm, x, n, nex, ney)
        implicit         none
        integer n, i, j, te, wm
        integer nex, ney, tn, prob
        double precision x(n), l, p
        character (len=100) :: temp

        tn = (2*nex+1)*(ney+1)+ney*(nex+1)

        j = 1
        open(unit=2, file="topo.dat", action="write",status="replace")
        do i = 1,ney
        write(2, '(1600E20.7)') (SNGL(x(2*tn + (i-1)*nex + j)),j =1,nex)
        enddo
        close(2)

        write(temp, 999) nex, ney

        CALL SYSTEM(temp)

        return;

!999     format ('gnuplot -e "xi =', I0,';yi =',I0, 
!     &  '" topoprint.p 2>&1')
999     format ('gnuplot -e "xi =', I0,';yi =',I0, 
     &  '" lopoprint_star.p 2>&1 | grep -v Warning:')

        end subroutine

! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: the sparse driver for the galahad lsqp QP solver; incomplete...
! SAOi:


       subroutine drive_galQPs (n,x,m,ne,f_a,c_a,
     &                          iact,nact,xlam,x_h,
     &                          acurv,bcurv,x_l,x_u,
     &                          f,c,gf,gc,ksubiter,
     &                          xkkt,flam,cstage0,ostring,
     &                          ictrl,lctrl,rctrl,cctrl,outerloop,
     &                          nnz,nnzh,Acol,Aptr,yhi,message,
     &                          eqn,lin,z,Hne,Hval,Hrow,Hcol,
     &                          iuser,luser,cuser,ruser)

! This driver calls the Galahad QP solvers (double)
      implicit          none
      include           'ctrl.h'
      character*60      task,csave
      character*90      ostring
      logical           kkt_local,cstage0
      logical           eqn(*), lin(*)
      integer           message,ifree
      integer           ne,n,ksubiter,m,loop,outerloop,i1,ij
      integer           i,j,nact,iact(nact),nnz
      integer           Arow(nnz),Acol(nnz),Aptr(m+1)
      double precision  flambda,factr,pgtol,yhi
      double precision  time1,time2
      double precision  x(n),flam
      double precision  x_h(n),x_l(n),x_u(n),xi(n),x_li(n),x_ui(n)
      double precision  xlam(m)
      double precision  acurv,bcurv(m) 
      double precision  xlamsml,xlambig
      double precision  f,c(m),gf(n),gc(nnz)
      double precision  f_a,c_a(m) 
      double precision  xkkt,gf_a(n),gc_a(nnz),gcji ! gc = Aval
!
      character*16      Astring
      integer           ierr
      double precision  w(n),xlamj,bji,htol
      double precision  y(m),z(n),xlam_h(m)
      double precision  c_l(m),c_u(m),g(n)
!
      character*16      Hstring
      integer           Hrow(nnzh),Hcol(nnzh),Hptr(n+1)
      integer           Hne,nnzh
      double precision  Hval(nnzh)
!
      data              kkt_local /.false./
!
      include           'ctrl_get.inc'
!
      if (n.gt.nmax) stop ' n.gt.nmax in drive_galQPs'
!      
      ierr = 0 
!
      QPform_A = 2 ! 0 = dense                ! not supported; not sensible in sparse
                   ! 1 = coordinate           ! not yet supported
                   ! 2 = sparse_by_rows     
!
      QPform_H = 3 ! 0 = dense                                   ! for debugging only
                    ! 1 = coordinate, diagonal, i.e. pseudo sparse  
                    ! 2 = sparse_by_rows, diagonal, i.e. pseudo sparse 
                    ! 3 = truly diagonal 
                    ! 4 = coordinate, not diagonal
                    ! 5 = sparse_by_rows, not diagonal
                    ! 6 = sparse_by_rows,User supplied
!                   
      htol = 1.d-6
!
      Arow=0
!      
      Hptr=0
!
! pre-process
      call QP_pres (n,m,c,acurv,
     &              bcurv,xlam,gc,x,xi,x_l,x_u,x_li,x_ui,
     &              z,c_l,c_u,w,Arow,Acol,Aptr,
     &              Astring,ictrl,
     &              lctrl,rctrl,cctrl,gf,x_h,nnz,nnzh,eqn,lin,
     &              Hstring,Hrow,Hcol,Hptr,
     &              Hval,Hne,QPform_A,QPform_H,htol,
     &              iuser,luser,cuser,ruser)

! call the Galahad QP solvers (double)
      if     (subsolver.eq.25) then !linear or separable convex quadratic program
        call galahad_lsqp (n,m,nnz,f_a,c_l,c_u,x_li,
     &                     x_ui,Arow,Acol,Aptr,xi,gf,
     &                     w,xlam,z,gc,Astring,ierr,
     &                     ksubiter,f)
!
      elseif (subsolver.eq.26) then !use working/active set method
        call galahad_qpa  (n,m,nnz,f_a,c_l,c_u,x_li,
     &                     x_ui,Arow,Acol,Aptr,xi,gf,
     &                     w,xlam,z,gc,Astring,ierr,
     &                     ksubiter,Hstring,Hrow,Hcol,
     &                     Hptr,Hval,Hne,f)
!
      elseif (subsolver.eq.27) then !primal-dual interior-point trust-region
        call galahad_qpb  (n,m,nnz,f_a,c_l,c_u,x_li,
     &                     x_ui,Arow,Acol,Aptr,xi,gf,
     &                     w,xlam,z,gc,Astring,ierr,
     &                     ksubiter,Hstring,Hrow,Hcol,
     &                     Hptr,Hval,Hne,f)
!
      elseif (subsolver.eq.28) then !uses a crossover method
        call galahad_qpc  (n,m,nnz,f_a,c_l,c_u,x_li,
     &                     x_ui,Arow,Acol,Aptr,xi,gf,
     &                     w,xlam,z,gc,Astring,ierr,
     &                     ksubiter,Hstring,Hrow,Hcol,
     &                     Hptr,Hval,Hne,f)
!
      endif

! update multipliers
      do j=1,m
        xlam(j)= max(-biglam,min(-xlam(j),biglam))
      enddo
      do i=1,n
        z(i)= max(-biglam,min(z(i),biglam))
      enddo

! update x in current point
      do i=1,n
        x(i)=min(x_u(i),max(x_l(i),x(i)+xi(i)))
      enddo
!
      message = ierr
!
      return
      end
!----------------------------------------------------------------------!
      subroutine QP_pres (n,m,c,acurv,
     &                   bcurv,xlam,gc,x,xi,x_l,x_u,x_li,x_ui,
     &                   z,c_l,c_u,w,Arow,Acol,Aptr,
     &                   Astring,
     &                   ictrl,lctrl,rctrl,cctrl,gf,x_h,nnz,nnzh,
     &                   eqn,lin,Hstring,Hrow,Hcol,Hptr,Hval,
     &                   Hne,QPform_A,QPform_H,htol,
     &                   iuser,luser,cuser,ruser)
      implicit none
      include          'ctrl.h'
      logical          eqn(*), lin(*)
      integer          n,m,ierr,i,j,k
      integer          ij,i1,nnz
      integer          Arow(*),Acol(*),Aptr(*)
      double precision w(*),xlamj,bji,x(*),xi(*),x_li(*),x_ui(*)
      double precision x_l(*),x_u(*),xlam(*),gc(*),z(*),c(*)
      double precision c_l(*),c_u(*),acurv,bcurv(*),bcurvn(m,n)
      double precision gcji,gf(*),x_h(*),htol
      character*16     Astring
!
      character*16     Hstring
      integer          Hrow(nnzh),Hcol(nnzh),Hptr(n+1)
      integer          nnzh,Hne
      double precision Hval(nnzh) 
!
      include          'ctrl_get.inc'
      
      QPform_H = 3
! create the Jacobian related entries
      if (QPform_A.ne.2) stop 'QPform_A.ne.2) in QP_pres'
      Astring = 'SPARSE_BY_ROWS' ! the only currently supported option !

! initialize n side constraints
!       do i=1,n
!         z(i) = 0.d0
!       enddo

! m constraint lower bounds
      do j=1,m
        if (eqn(j)) then
          c_l(j) = -c(j)
        else
          c_l(j) = -1.d20
        endif
      enddo

! m constraint upper bounds
      do j=1,m
        c_u(j) = -c(j)
      enddo

! Construct the diagonal Hessian for SAOi
      do i=1,n
        if (approx_f.eq.1) then
          !w(i) = max(atol1,acurv) ! SQ - formed in diaHess_sq1s
        elseif (approx_f.eq.4) then
          w(i) = max(atol1,2.d0/x_h(i)*dabs(gf(i))) ! T2:R
        elseif (approx_f.eq.14) then
          w(i) = max(atol1,-2.d0/x_h(i)*(gf(i))) ! T2:R
        endif
      enddo
      ij=0
      bji=0.d0 !used for "Conditional jump or move depends on uninitialised value(s)"
      xlamj=0.d0 !used for "Conditional jump or move depends on uninitialised value(s)"
      do j=1,m
        xlamj = xlam(j)
        do i1=1,Aptr(j+1)-Aptr(j)
          ij=ij+1
          i=Acol(ij)
          gcji=gc(ij)
          if (approx_c.eq.1) then
            !bji=max(btol1,bcurv(j)) ! SQ - formed in diaHess_sq1s
          elseif (approx_c.eq.4) then
            bji=max(btol1,2.d0/x_h(i)*dabs(gcji)) ! T2:R
          elseif (approx_c.eq.14) then
            bji=max(btol1,-2.d0/x_h(i)*(gcji)) ! T2:R
          endif
          if (bji.ne.0.d0.and.xlamj.ne.0.d0)then
            w(i) =  w(i) + bji*xlamj
          endif
        enddo
      enddo
      
!  Construct the diagonal Hessian for SAOi - we DO NOT work with acurvn, bcurvn
      if (approx_f.eq.1.and.approx_c.eq.1) then
        do i=1,n
          w(i) = acurv
          do j=1,m
            bji   = bcurv(j)
            xlamj = xlam(j)
            if (sph_dia_hess.eq.0) then     ! all approx incl. SQ form 0 - every term positive
              w(i) = max(htol,acurv)
              if (bji.ge.0.d0.and.xlamj.gt.0.d0) then
                w(i) =  w(i) + bji*xlamj
              endif  
            elseif (sph_dia_hess.eq.1.or.sph_dia_hess.eq.2) then
              if (xlamj.ne.0.d0) then
                w(i) =  w(i) + bji*xlamj
              endif
            else
              stop ' sph_dia_hess not defined in dgalQPs.f'
            endif
          enddo
          if (sph_dia_hess.eq.1) w(i) = max(htol,w(i))        ! all approx incl. SQ form 1 - final sum positive
          if (sph_dia_hess.eq.2) w(i) = max(htol,dabs(w(i)))  ! all approx incl. SQ form 2 - with abs operator 
        enddo
      elseif (approx_f.eq.100.and.approx_c.eq.100) then
        call Sep_gg(x,n,m,xlam,w)! Put in Thu 16 Nov 2017 12:23:45 SAST
      elseif (approx_f.eq.101.and.approx_c.eq.101) then
        QPform_H = 6
      endif  
      
! Transform the Hessian to one of four QP forms.
! 1, 2 are not really sensible for diagonal; only given below for illustration
!      
      if (QPform_H.eq.0) then
        if (nnzh.ne.n*(n+1)/2) then
          write(*,*) cname1(1:6),QPform_H, ' nnzh.ne.n*(n+1)/2',
     &                   ' in QP_pres',nnzh,n*(n+1)/2
          stop
        endif
      elseif (QPform_H.ge.1.and.QPform_H.le.3) then
        if (nnzh.ne.n) then
          write(*,*) cname1(1:6),QPform_H, ' nnzh.ne.n in QP_pres',
     &                   nnzh,n*(n+1)/2
          stop
        endif
      else 
        if (nnzh.lt.0) then
          write(*,*) cname1(1:6),QPform_H, ' nnzh.lt.0 in QP_pres',
     &                   nnzh 
          stop
        endif  
        if (nnzh.gt.n*(n+1)/2) then
          write(*,*) cname1(1:6),QPform_H, ' nnzh.gt.n*(n+1)/2',
     &                   ' in QP_pres',nnzh,n*(n+1)/2
          stop 
        endif   
      endif
!
!
      if (QPform_H.eq.0) then ! pseudo diagonal
        Hstring = 'DENSE'
        Hne = nnzh
!
      elseif (QPform_H.eq.1) then ! pseudo diagonal
        Hstring = 'COORDINATE'
        Hne = n
        do i=1,n
          Hval(i) = w(i)
          Hrow(i) = i
          Hcol(i) = i
        enddo
!
      else if (QPform_H.eq.2) then ! pseudo diagonal
        Hstring = 'SPARSE_BY_ROWS'
        Hne = n
        do i=1,n
          Hval(i) = w(i)
          Hcol(i) = i
          Hptr(i) = i
        enddo
        Hptr(n+1) = n+1
!
      else if (QPform_H.eq.3) then ! truly diagonal
        Hstring = 'DIAGONAL'
        Hne = n
        do i=1,n
          Hval(i) = w(i)
          Hrow(i) = i
          Hcol(i) = i
        enddo
!
      else if (QPform_H.eq.4) then ! coordinate - form in SparseHessian.f
        Hstring = 'COORDINATE'
!
      else if (QPform_H.eq.5) then ! sparse_by_rows - form in SparseHessian.f
        Hstring = 'SPARSE_BY_ROWS'

      else if (QPform_H.eq.6) then ! Supplied by user
        Hstring = 'COORDINATE'
        call Fem_gg(x,n,m,xlam,nnzh,Hval,Hptr,Hcol,Hne)
        k=1
        do i=1,n
          do j=1,Hptr(i+1)-Hptr(i)
          Hrow(k) = i
          k=k+1
          enddo
        enddo
      else
       stop 'QPform_H.ne.0,1,2,3,4,5,6) in QP_pres'
!
      endif      

! Initialise move limits
      do i=1,n
        x_li(i)=x_l(i)-x(i)
        x_ui(i)=x_u(i)-x(i)
      enddo

! create copy of x
      call dcopy(n,x,1,xi,1)
!
      return
      end subroutine QP_pres
!----------------------------------------------------------------------!
      subroutine FormHessSparseD (n,m,x,xlam,Hne,Hval,Hrow,Hcol,
     &                            ictrl,lctrl,rctrl,cctrl,iuser,luser,
     &                            cuser,ruser,nnzh,eqn,lin,acurv)
      implicit none
      include           'ctrl.h'
      logical           eqn(*),lin(*)
      integer           n,m,i,j,nnzh 
      integer           Hrow(nnzh),Hcol(nnzh),Hptr(n+1),Hne
      double precision  x(n),xlam(m),Hval(nnzh),acurv
!
      integer           allocatestatus, l, k, p, d, k1, a, i1
      integer           iscale1, iscale2
      double precision, dimension(:),    allocatable :: xs
      double precision, dimension(:, :), allocatable :: toeplitz
      include           'ctrl_get.inc'
!
      k       = iuser(1)  !  problem specific - number of requested eigen values
      p       = iuser(2)  !  problem specific - length of each eigen vector column  
      iscale1 = iuser(3)  !  problem specific 
      iscale2 = iuser(4)  !  problem specific 
!
      allocate ( xs(p), stat = allocatestatus)
      if (allocatestatus /= 0) 
     &          STOP " not enough memory to allocate xs"
!
      allocate ( toeplitz(p, p), stat = allocatestatus)
      if (allocatestatus /= 0) 
     &          STOP " not enough memory to allocate toeplitz"
!
      Hne = nnzh
!
      xs    = 0.d0
      xs(1) = 4.d0
      xs(2) = 1.d0
!
      do i = 1,p
        toeplitz(i,i:p) = xs(1:(p+1-i))
        toeplitz(i:p,i) = xs(1:(p+1-i))
      end do  
!      
!       do i=1,p 
!         write(*,*) i,(toeplitz(i,j),j=1,p)
!       enddo
!       stop
!
!
      Hval = 0.d0
!
!  Construct objective Hessian
!
      i = 0
      do l = 1,k
        do j = 1,p
          i = i + 1
          Hrow(i) = i
          Hcol(i) = i
!          
          Hval(i) = -2.d0*toeplitz(j,j)*(k+iscale1*l-iscale2)       ! -2 since we maximize toeplitz...
          !Hval(i) =  1.d0*max(atol1,acurv) ! *(k+iscale1*l-iscale2) 
          !Hval(i) =  1.d0*max(0.d0,acurv)  ! *(k+iscale1*l-iscale2) 
          !Hval(i) =  1.d0*max(-1.d3,acurv) ! *(k+iscale1*l-iscale2) 
!            
!
!           gf(i) = -2.d0*temp(j)*dble(k+iscale1*l-iscale2)
!           
!           write(*,*) l,i,i,
!      &               nint(-2.d0*toeplitz(j,j)*(k+iscale1*l-iscale2))
!          
        enddo  
      enddo 
!
!  Construct constraint gradients
!
!      gc = 0.d0
!
!  Enforce off-diagonal terms = 0
!       
      i1 = 0

      i=0 
      do l=1,k-1
        do d = l+1,k
          i = i+1
          do j = 1,p
          
          i1 = i1 + 1
          
          Hrow(i1) = (d-1)*p+j  
          Hcol(i1) = (l-1)*p+j 
          Hval(i1) = 1.d0*xlam(i)
!          
!             gc(i,(l-1)*p+j) = temp(j)
!             gc(i,(d-1)*p+j) = mode(j)
!             
              !write(*,1000) i1,i,(l-1)*p+j,(d-1)*p+j,1
!             
          enddo
        enddo  
      enddo  
!
!  Enforce diagonal terms - 1 = 0
!
      k1=k*(k-1)/2
      do i=1,k
        do j=1,p
          
          i1 = i1 + 1
          
          Hrow(i1) = (i-1)*p+j 
          Hcol(i1) = (i-1)*p+j 
          Hval(i1) = 2.d0*xlam(i+k1)
!
!           gc(i+k1,(i-1)*p+j) = 2.d0*mode(j)
!           
            !write(*,1000)i1,i+k1,(i-1)*p+j,(i-1)*p+j,2
!           
        enddo
      enddo      
!      
 1000 format(1000i4)
!            

       do i = 1,nnzh
         !write(*,*) i,Hrow(i),Hcol(i),Hval(i)
       enddo  

       !stop

!            
!  Deallocate some stuff here ...
!
      end subroutine FormHessSparseD
!----------------------------------------------------------------------!
      subroutine FormHessSparseN (n,m,x,xlam,Hne,Hval,Hrow,Hcol,
     &                            ictrl,lctrl,rctrl,cctrl,iuser,luser,
     &                            cuser,ruser,nnzh,eqn,lin,acurv)
      implicit none
      include           'ctrl.h'
      logical           eqn(*),lin(*)
      integer           n,m,i,j,nnzh 
      integer           Hrow(nnzh),Hcol(nnzh),Hptr(n+1),Hne
      double precision  x(n),xlam(m),Hval(nnzh),acurv,temp
!
      integer           allocatestatus, l, k, p, d, k1, a, i1
      integer           iscale1, iscale2, nnzh_off
      double precision, dimension(:),    allocatable :: xs
      double precision, dimension(:, :), allocatable :: toeplitz
      include           'ctrl_get.inc'
!
      k       = iuser(1)  !  problem specific - number of requested eigen values
      p       = iuser(2)  !  problem specific - length of each eigen vector column  
      iscale1 = iuser(3)  !  problem specific 
      iscale2 = iuser(4)  !  problem specific 
!
      allocate ( xs(p), stat = allocatestatus)
      if (allocatestatus /= 0) 
     &          STOP " not enough memory to allocate xs"
!
      allocate ( toeplitz(p, p), stat = allocatestatus)
      if (allocatestatus /= 0) 
     &          STOP " not enough memory to allocate toeplitz"
!
      xs    = 0.d0
      xs(1) = 4.d0
      xs(2) = 1.d0
!
      do i = 1,p
        toeplitz(i,i:p) = xs(1:(p+1-i))
        toeplitz(i:p,i) = xs(1:(p+1-i))
      end do  
!      
!       do i=1,p 
!         write(*,*) i,(toeplitz(i,j),j=1,p)
!       enddo
!       stop
!
!
      Hval = 0.d0
!
!  Construct objective Hessian
!

!
!  Construct constraint gradients
!
!      gc = 0.d0
!
!  Enforce off-diagonal terms = 0
!       

      !do i = 1,m
      !  xlam(i) = dble(i) ! 1.d0
      !enddo  

      i1 = 0

      i=0 
      do l=1,k-1
        do d = l+1,k
          i = i+1
          do j = 1,p
          
          i1 = i1 + 1
          
          Hrow(i1) = (d-1)*p+j  
          Hcol(i1) = (l-1)*p+j 
          Hval(i1) = 1.d0*xlam(i)
!          
!             gc(i,(l-1)*p+j) = temp(j)
!             gc(i,(d-1)*p+j) = mode(j)
!             
              !write(*,1000) i1,i,(l-1)*p+j,(d-1)*p+j,1
!             
          enddo
        enddo  
      enddo  
!
!  Enforce diagonal terms - 1 = 0
!
      k1=k*(k-1)/2
      do i=1,k
        do j=1,p
          
          i1 = i1 + 1
          
          Hrow(i1) = (i-1)*p+j 
          Hcol(i1) = (i-1)*p+j 
          Hval(i1) = 2.d0*xlam(i+k1) 
!
!           gc(i+k1,(i-1)*p+j) = 2.d0*mode(j)
!           
            !write(*,1000)i1,i+k1,(i-1)*p+j,(i-1)*p+j,2
!           
        enddo
      enddo      
!      
!  1000 format(1000i4)
!            


!
      nnzh_off = 0
!
      !i1 = 0
      do l=1,k
        do a=1,p 
          do j=1,a
            temp = toeplitz(a,j)
            if (a+(l-1)*p.eq.j+(l-1)*p) then
            
              i1 = i1 + 1
              Hrow(i1) = a+(l-1)*p
              Hcol(i1) = j+(l-1)*p 
              Hval(i1) = - 2.d0*temp*(k+iscale1*l-iscale2)
            
            else
            
              if (temp.ne.0.d0) then
                nnzh_off = nnzh_off + 1
              
                i1 = i1 + 1
                Hrow(i1) = a+(l-1)*p
                Hcol(i1) = j+(l-1)*p 
                Hval(i1) = - 2.d0*temp*(k+iscale1*l-iscale2)
              
              endif  
              
!               write(*,*) 
!      &        i1,nnzh_off,a+(l-1)*p,j+(l-1)*p,l,
!      &        -2.d0*temp*(k+iscale1*l-iscale2)

            endif
          enddo
        enddo  
      enddo      
!      
      Hne = nnzh
!

      !nnzh = nnzh + nnzh_off
      
      !do i = 1,nnzh
      !  write(*,*) i,Hrow(i),Hcol(i),Hval(i)
      !enddo  

      !stop
      
!            
!  Deallocate some stuff here ...
!
      end subroutine FormHessSparseN
!----------------------------------------------------------------------!

! SAOi: 
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Finite-differences, random numbers, kkt conditions, and many
! SAOi: other utilities. Includes a few badly commented, badly
! SAOi: constructed BLAS hacks, not to be used lightly
! SAOi:


      double precision function conv(it,ir)
      implicit         none
      integer          it,ir,i
      double precision c1,c2,one,tmp
      data             one /1.d0/
!
      c1=dble(it-ir)
      c2=dble(it+1)
      tmp=one
      do i=1,it
        c1=c1+one
        c2=c2+one
        tmp=tmp*c1/c2
      enddo
!
      conv=tmp
!
      return
      end function conv
!----------------------------------------------------------------------!
      subroutine open_warn_file (warn)
!     include          'ctrl.h' - do not use here
      implicit         none
      logical          warn
!
      if (warn) return
!
      warn = .true.
      open (14,file='Warnings.out')
      write(14,1000)
!
      return
 1000 format (/,26x,'SAOi algorithm: warning and errors',//)
      end subroutine open_warn_file
!----------------------------------------------------------------------!
      subroutine SAOi_open (ictrl,lctrl,rctrl,cctrl,
     &                      iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i
!
      open (8 ,file='Variables.out')
      open (9 ,file='History.out')
      open (10,file='Tolerance-X.out')
      open (11,file='Tolerance-KKT.out')
      open (12,file='Subproblems.out')
      open (13,file='Constraints.out')
!
!     unit 14 is reserved for possible warnings
!
      open (15,file='Specialized.out') ! reserved for the test bed
      open (16,file='test.out')
      close (16)
!
      ictrl = 0
      lctrl = .false.
      rctrl = 0.d0
      cctrl = '------------------------'
!
      iuser = 0
      luser = .false.
      ruser = 0.d0
      cuser = '////////////////////////'
!
      return
      end subroutine SAOi_open
!----------------------------------------------------------------------!
      subroutine SAOi_close ()
!     include          'ctrl.h' - do not use here
      implicit         none
      integer          i
!
      do i=8,17    !     no need to close 16
        close (i)
      enddo
!
      return
      end subroutine SAOi_close
!----------------------------------------------------------------------c
      subroutine smma (n,ni,ne,x,x_h,x_h2,x_u,x_l,s,ictrl,lctrl,rctrl,
     &                 cctrl,outerloop,iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,outerloop
      double precision x(n),x_h(n),x_h2(n),delx(n),x_u(n),x_l(n),s(n)
      double precision diff,osc,gammai
      double precision zero,one,two
      data             zero /0.d0/
      data             one  /1.d0/
      data             two  /2.d0/
      include          'ctrl_get.inc'
!
      do i = 1,n
        if (outerloop.lt.2) then
          s(i) = (x_u(i)-x_l(i))/two
        else
          osc = (x(i)-x_h(i))*(x_h(i)-x_h2(i))
          if (osc.lt.zero) then
            gammai = mma_lo
          elseif (osc.gt.zero) then
            gammai = mma_hi
          else
            gammai = one
          endif
          s(i) = s(i)*gammai
        endif
        s(i) = max(s(i),(x_u(i)-x_l(i))/1.d2)
        s(i) = min(s(i),(x_u(i)-x_l(i))*1.d1)
      enddo
!
      return
      end subroutine smma
!----------------------------------------------------------------------!
      subroutine SAOi_time (time0,timef,timeg,times,
     &                      ictrl,lctrl,rctrl,cctrl,
     &                      iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,ipr
      double precision time0, timef, timeg, times, timet
      double precision t1, t2, t3, t4, t5, SAOi_seconds
      include          'ctrl_get.inc'
!
      timet = SAOi_seconds()
      t1=0.d0
      t2=0.d0
      t3=0.d0
      t4=0.d0
      t5=0.d0
      if (timet-time0.gt.0.d0) then
        t1=timef/(timet-time0)*100.d0
        t2=timeg/(timet-time0)*100.d0
        t3=(timef+timeg)/(timet-time0)*100.d0
        t4=times/(timet-time0)*100.d0
        t5=100.d0-t4-t3
      endif
      do ipr=6,9,3
        if ((ipr.eq.6.and.iprint.ge.2).or.ipr.eq.9) then
          write(ipr,1000) timet-time0,
     &                    timef,t1,timeg,t2,timef+timeg,t3,times,t4,
     &                    (timet-time0)-(timef+timeg)-times,t5
        endif
      enddo
!
!     open  (16,file='test.out',status='old',position='append')
!     write (16,*) timet-time0
!     close (16)
!
      return
!
 1000 format(/,'  Elapsed time (wall clock seconds) = ',1f10.2,' ',/,
     &         30x,1f18.2,' (',1f6.2,'%) in evaluating f_j ',/,
     &         30x,1f18.2,' (',1f6.2,'%) in evaluating df_j ',/,/,
     &         30x,1f18.2,' (',1f6.2,'%) in evaluating f_j and df_j ',/,
     &         30x,1f18.2,' (',1f6.2,'%) in solving the subproblems',/,
     &         30x,1f18.2,' (',1f6.2,'%) in overheads',/)
!
      end subroutine SAOi_time
!----------------------------------------------------------------------!
      subroutine reset_conv (ictrl,lctrl,rctrl,cctrl,
     &                       iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          ifc
      equivalence      (ifc,force_converge)
      include          'ctrl_get.inc'
!
      conservative    = .false.
      trust_region    = .false.
      classic_conserv = .true.
      classic_trust   = .true.
!
      if (ifc.lt.0.or.ifc.gt.3) ifc = 0
!
      if (ifc.eq.1) conservative = .true.
!
      if (ifc.eq.2.or.ifc.eq.3) trust_region = .true.
!
      if (ifc.eq.3) classic_trust = .false.
!
      include          'ctrl_set.inc'
      return
      end subroutine reset_conv
!----------------------------------------------------------------------!
!                                                                      !
      subroutine formgrad (n,x,ni,ne,f,c,gf,gc,ictrl,lctrl,rctrl,
     &                     cctrl,nfe,nge,do_diff,iuser,luser,cuser,
     &                     ruser,nnz,Acol,Aptr,eqn,lin)
!                                                                      !
!     compute the gradient vectors                                     !
!                                                                      !
!----------------------------------------------------------------------!
      implicit          none
      include           'ctrl.h'
      logical           eqn(*), lin(*)
      logical           do_diff
      integer           i,j,k,l,icnt,m,n,ni,ne,nnz,nfe,nge
      integer           Acol(*),Aptr(*)
      integer           allocatestatus,ierr,go
      integer           ir(nnz),jc(nnz)        ! for coo
      double precision  gco(nnz)               ! for coo
      double precision  x(n),tmp,tmp2,glarge
      double precision  dx(n),c_x(ni),c_d(ni),f_x,f_d
      double precision  f,c(ni),gf(n),gc(nnz)
      double precision, dimension(:), allocatable :: gcd
!
!
!
! ni here is still the total no of constraints...       
! testing sparse gradients is only sensible for small instances of the presumably scalable test problems...  
!
! changing the logic here to only check specified non-zero entries is non-sensical, since entries may be missed...
! catch 22 - contact  -  albertg@sun.ac.za  -  if you really need to do the above, but I cannot promise anything...
!
!
!      
      include          'ctrl_get.inc'                  
!
      do_diff=.false.
!
      if (.not.finite_diff.or.check_grad) then
        call SAOi_grads (n,ni,ni-ne,ne,x,gf,gc,nnz,Acol,Aptr,iuser,
     &                   luser,cuser,ruser,eqn,lin,ictrl,lctrl,
     &                   rctrl,cctrl)
      endif
!
      if (finite_diff.or.check_grad) then
!
        do_diff=.true.
!
        if (check_grad) then
          glarge=0.d0
          write(17,1000)
        endif
!
        if (structure.lt.3) then
!
          call SAOi_funcs(n,ni,ni-ne,ne,x,f_x,c_x,iuser,luser,
     &                  cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
          nfe=nfe+1          
          dx=x
!
          do i=1,n
            dx(i)=x(i)+deltx
            call SAOi_funcs(n,ni,ni-ne,ne,dx,f_d,c_d,iuser,luser,
     &                    cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
            dx(i)=x(i)
            nfe=nfe+1
            tmp=gf(i)
            gf(i)=(f_d-f_x)/deltx
            if (check_grad) then
              tmp2=dabs((tmp-gf(i))/(dabs(tmp)+1.d0))
              glarge=max(glarge,tmp2)
              write(17,1001)i,tmp,gf(i),tmp2,glarge
              gf(i)=tmp
            endif
            do j=1,ni
              k=(i-1)*ni+j
              tmp=gc(k)
              gc(k)=(c_d(j)-c_x(j))/deltx
              if (check_grad) then
                tmp2=dabs((tmp-gc(k))/(dabs(tmp)+1.d0))
                glarge=max(glarge,tmp2)
                write(17,1002)i,j,k,tmp,gc(k),tmp2,glarge
                !gc(k)=tmp
              endif
            end do
          end do
!
        elseif (structure.eq.3.and.strict_struct) then
!
          if (finite_diff) then
            write(*,1003) 
            stop
          endif
!
          if (Aptr(ni+1).ne.(nnz+1)) then
            write(*,*) ' Terminal: Aptr(ni+1).ne.(nnz+1) in CSR',
     &                 ' in formgrad',Aptr(ni+1),nnz+1
            stop
          endif
          if (Aptr(1).ne.(1)) then
            write(*,*) ' Terminal: Aptr(1).ne.(1) in CSR in formgrad'
            stop
          endif
!
          allocate ( gcd(ni*n), stat = allocatestatus)
          if (allocatestatus /= 0)
     &           STOP " not enough memory to allocate gcd in formgrad"
!
          call SAOi_funcs(n,ni,ni-ne,ne,x,f_x,c_x,iuser,luser,
     &                  cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
          nfe=nfe+1          
          dx=x
!
          gcd=0.d0
          icnt=0
!
          do i=1,n
            dx(i)=x(i)+deltx
            call SAOi_funcs(n,ni,ni-ne,ne,dx,f_d,c_d,iuser,luser,
     &                    cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
            dx(i)=x(i)
            nfe=nfe+1
            tmp=gf(i)
            gf(i)=(f_d-f_x)/deltx
            if (check_grad) then
              tmp2=dabs((tmp-gf(i))/(dabs(tmp)+1.d0))
              glarge=max(glarge,tmp2)
              write(17,1001)i,tmp,gf(i),tmp2,glarge
              gf(i)=tmp
            endif
!            
            do j=1,ni
              k=(i-1)*ni+j
              !tmp=gc(k)
              gcd(k)=(c_d(j)-c_x(j))/deltx
              if(gcd(k).ne.0.d0) icnt=icnt+1
              !if (check_grad) then
                !tmp2=dabs((tmp-gc(k))/(dabs(tmp)+1.d0))
                !glarge=max(glarge,tmp2)
                !write(17,1002)i,j,k,tmp,gc(k),tmp2,glarge
                !gc(k)=tmp
              !endif
            end do
          end do
!
          call csrcoo (ni,3,nnz,gc,Acol,Aptr,nnz,gco,ir,jc,ierr) ! can be done partly in-place, but helps little, 
!                                                                !    since gcd is so expensive...
!
          if (ierr.ne.0) 
     &        STOP "csrcoo format conversion failed in formgrad"
!           
          do i = 1,nnz
            tmp=gcd((jc(i)-1)*ni+ir(i))
            tmp2=dabs((tmp-gco(i))/(dabs(tmp)+1.d0))
            glarge=max(glarge,tmp2)
            write(17,1002) jc(i),ir(i),i,gco(i),tmp,tmp2,glarge
          enddo
          write(17,1004)nnz,icnt,nnz-icnt
          if (nnz-icnt.ne.0) then
            write(17,1005)glarge
            write(*,1005)glarge
            go=0
            !write(*,1006)
            !read(*,*) go
            if (go.ne.0) stop
          endif
!
        elseif (.not.strict_struct) then
!
          stop '.not.strict_struct or structure not defined' ! maybe add this; easy (pointers simply change each iteration) !
!
        else
!
          stop 'strict_struct or structure not defined' ! keep; this is non-sensical ! 
!
        endif

      endif
!
      if (check_grad) then
        write (*,*)  ' The largest relative error with delta x = ',
     &               deltx,' is ',glarge
        if (glarge.gt.5.d-4) then
          stop
          ! read (*,*) ! maybe do something here ... !
        endif
      endif
!
      do_diff=.false.
!
      return
!
 1000 format(/,'       n       m       k    ',
     &      '    df & dc coded       df & dc fin-diff   ',
     &      '      rel. error        largest rel. error',/)
 1001 format(1i8,16x,4e22.10)
 1002 format(3i8,4e22.10)
 1003 format(/,' Requesting finite differences in sparse algebra',
     &         ' is not sensical.',/,' Checking sparse gradients', 
     &         ' however may be, if only for small instances of',
     &         ' a presumably scalable test problem.',/)
 1004 format(/,' nnz coded      = ',i9,
     &       /,' nnz calculated = ',i9,
     &       /,' difference     = ',i9,/,/)
 1005 format(/,' WARNING: nnz coded and nnz calculated differ.',
     &         ' This is not necessarily terminal.',
     &         ' The largest relative error is ',e14.4,/)
 1006 format(/,' Press 0 to continue ',/)
!   
      end subroutine formgrad
!----------------------------------------------------------------------c
      subroutine SAOistatus(n,ni,iactblo,iactbhi,iactc,cplus_norm,xdual,
     &                      violation,c,x,x_l0,x_u0,ic,ig,hsum,
     &                      ictrl,lctrl,rctrl,cctrl,feasible,
     &                      iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      logical          feasible
      integer          i,j,n,ni,iactblo,iactbhi,iactc,ic,ig
      double precision x(n),x_l0(n),x_u0(n),c(ni),violation
      double precision cplus_norm,act_lim_prnt,xdual(ni),hsum
      include          'ctrl_get.inc'
!
      feasible=.true.
      iactblo=0
      iactbhi=0
      iactc=0
      violation=0.d0
      act_lim_prnt=feaslim ! / 1.d2
      cplus_norm=0.d0
!
      if (ic.lt.5) then
!
        violation=-big
        hsum=0.d0
        do i=1,ni
          hsum=hsum+max(0.d0,c(i))
          cplus_norm=cplus_norm+(max(0.d0,c(i)))**2
          if (c(i).gt.violation) violation=c(i)
          if (c(i).gt.-act_lim_prnt) iactc=iactc+1
        end do
        cplus_norm=dsqrt(cplus_norm)
        if (violation.gt.feaslim) feasible=.false.
      endif
!
      if (ig.lt.0) return
!
      do i=1,n
        if (dabs(x(i)-x_l0(i)).lt.1.d-8) iactblo=iactblo+1
        if (dabs(x_u0(i)-x(i)).lt.1.d-8) iactbhi=iactbhi+1
      end do
!
      return
      end subroutine SAOistatus
!----------------------------------------------------------------------c
      subroutine SAOistatusE(n,ni,iactblo,iactbhi,iactc,cplus_norm, 
     &                       xdual,violation,c,x,x_l0,x_u0,ic,ig,hsum,
     &                       ictrl,lctrl,rctrl,cctrl,feasible,eqn,lin,
     &                       iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      logical          feasible
      logical          eqn(*),lin(*)
      integer          i,j,n,ni,iactblo,iactbhi,iactc,ic,ig
      double precision x(n),x_l0(n),x_u0(n),c(ni),violation
      double precision cplus_norm,act_lim_prnt,xdual(ni),hsum
      include          'ctrl_get.inc'
!
      violation=0.d0
      do j=1,ni
        if (eqn(j)) then
          violation=max(violation,dabs(c(j)))
        else
          violation=max(violation,c(j))
        endif
      enddo
      if (violation.gt.feaslim) feasible=.false.
!
      return
      end subroutine SAOistatusE
!----------------------------------------------------------------------c
      subroutine store_iter (n,ni,ne,f_h,f,v_h,viol,x_h,x_h2,x,c_h,c,
     &                       cplus_norm,cplus_norm_h,gf_h,gf,gc_h,gc,
     &                       h_h,h,gh_h,gh,xlam,xlam_h,ic,ig,z,z_h,
     &                       iuser,luser,cuser,ruser,ictrl,lctrl,
     &                       rctrl,cctrl)
      implicit         none
      include          'ctrl.h'
      integer          i,j,n,ni,ne,ic,ig
      double precision f_h,f,v_h,viol,x_h(n),x_h2(n),x(n),c_h(ni),c(ni)
      double precision gf_h(n),gf(n),gc_h(ni,n),gc(ni,n),h_h(ne),h(ne)
      double precision cplus_norm,cplus_norm_h,gh_h(ne,n),gh(ne,n)
      double precision xlam(ni),xlam_h(ni),z(n),z_h(n)
      include          'ctrl_get.inc'
!
      if (ic.eq.1) then
        cplus_norm_h=cplus_norm
        f_h=f
        v_h=viol
        call dcopy (n,x_h,1,x_h2,1)
        call dcopy (n,x,1,x_h,1)
        call dcopy (ni,c,1,c_h,1)
        call dcopy (ni,xlam,1,xlam_h,1)
        call dcopy (n,z,1,z_h,1)
!
        if (ig.eq.0) return
!
        call dcopy (n,gf,1,gf_h,1)
        call dcopy (n*ni,gc,1,gc_h,1)
!
      elseif (ic.eq.-1) then
        cplus_norm=cplus_norm_h
        f=f_h
        viol=v_h
        call dcopy (n,x_h,1,x,1)
        call dcopy (ni,c_h,1,c,1)
        call dcopy (ni,xlam_h,1,xlam,1)
        call dcopy (n,z_h,1,z,1)
!
        if (ig.eq.0) return
!
        call dcopy (n,gf_h,1,gf,1)
        call dcopy (n*ni,gc_h,1,gc,1)
!
      else
        stop ' undefined in store_it'
      endif
!
      return
      end subroutine store_iter
!----------------------------------------------------------------------c
      subroutine sparser (n,ni,gf,gc,ns,ictrl,lctrl,rctrl,cctrl,warn,
     &                    iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      logical          warn
      integer          i,j,n,ni,ns
      double precision gf(n),gc(ni,n),zero
      data             zero /0.d0/
      include          'ctrl_get.inc'
!
      if (structure.eq.2) then         ! pseudo
        ns=0
        do i=1,n
          do j=1,ni
            if (gc(j,i).ne.zero) ns=ns+1
          end do
        end do
      elseif (structure.eq.1) then     ! dense
        ns=n*ni
      elseif (structure.eq.3) then     ! sparse
        stop ' you should not be able to get here :-) '
      else
        stop ' illegal structure was defined '
      endif
!
!       if (1.d0-dble(ns)/dble(n*ni).lt.sfracmin.and.
!      &                                    structure.eq.2) then
!         !call open_warn_file (warn)
!         !write (14,1000) 1.d0-dble(ns)/dble(n*ni)
!         structure=1
!         ns=n*ni
!       endif
!
      return
 1000 format (' Pseudo-sparsity was killed, sparsity is only ',1f16.6)
      end subroutine sparser
!----------------------------------------------------------------------c
      subroutine sparset (n,ni,gf,gc,ns,ncol,nrow,gcvec,bvec,bcurvn,
     &                    ictrl,lctrl,rctrl,cctrl,
     &                    iuser,luser,cuser,ruser) !  still uses the COO storage scheme
      implicit         none
      include          'ctrl.h'
      integer          i,j,n,ni
      integer          ns,ncol(ns),nrow(ns),icnt
      double precision gcvec(ns),bvec(ns),zero
      double precision gf(n),gc(ni,n),bcurvn(ni,n)
      data             zero /0.d0/
      include          'ctrl_get.inc'
!
      if (structure.eq.2) then
        icnt=0
        do i=1,n
          do j=1,ni
            if (gc(j,i).ne.zero) then
              icnt=icnt+1
              ncol(icnt)=j
              nrow(icnt)=i
              gcvec(icnt)=gc(j,i)
              bvec(icnt)=bcurvn(j,i)
            endif
          enddo
        enddo
        if (icnt.ne.ns) then
          write(*,*) icnt,ns
          stop ' error 1 in sparset '
        endif
      endif
!
      return
      end subroutine sparset
!----------------------------------------------------------------------c
      subroutine form_act(ni,c,iact,nact,ictrl,lctrl,rctrl,cctrl,
     &                    iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          j,ni,nact,iact(ni)
      double precision c(ni)
      include          'ctrl_get.inc'
!
      nact=0
      if (use_active_set) then
        stop ' Active sets are not supported in the current version '
        do j=1,ni
          if (c(j).gt.-actlim) then
            nact=nact+1
            iact(nact)=j
          end if
        end do
      else
        nact=ni
        do j=1,ni
          iact(j)=j
        end do
      endif
!
      return
      end subroutine form_act
!----------------------------------------------------------------------c
      double precision function norm2b (n,a,b)
!
!  calculate the Euclidian or 2-norm between a[n] and b[n]
!
      implicit         none
      !include          'ctrl.h'
      integer          i, n
      double precision norm2, a(n), b(n), zero
      data             zero /0.d0/
!
      norm2 = zero
      do i = 1,n
        norm2 = norm2 + (a(i) - b(i))**2
      end do
      norm2b=dsqrt(norm2)
!
      return
      end function norm2b
!----------------------------------------------------------------------c
      double precision function norm2f (n,a,b)
!
!  calculate the infinity norm between a[n] and b[n]
!
      implicit         none
      !include          'ctrl.h'
      integer          i, n
      double precision a(n), b(n), zero
      data             zero /0.d0/
!
      norm2f = zero
      do i = 1,n
        norm2f = max(norm2f,dabs(a(i) - b(i)))
      end do
!
      return
      end function norm2f
!----------------------------------------------------------------------!
      subroutine form_kkt (xkkt,n,ni,ne,x,iact,nact,x_l0,x_u0,
     &                     gf,gc,gh,xlam,xlamsml,xlambig,ifree,
     &                     ictrl,lctrl,rctrl,cctrl,z,
     &                     iuser,luser,cuser,ruser)
!
! kkt conditions including the bound constraints via z               !
!
      implicit         none
      include          'ctrl.h'
      integer          i,j,ij,j1,n,ni,ne,nact,iact(nact),ifree
      double precision x(n),x_l0(n),x_u0(n),gf(n),gc(ni,n),gh(ne,n),gcji
      double precision xkkt,kkt(n),xlam(*),xlamsml,xlambig,xlamj,z(n)
!
      include 'ctrl_get.inc'
!
      ifree   = 0
      xkkt    = 0.d0
      xlamsml = 1.d20
      xlambig =-1.d20
!
      do j=1,ni
        xlamj=xlam(j)
        if (xlamj.gt.xlambig) xlambig=xlamj
        if (xlamj.lt.xlamsml) xlamsml=xlamj
      enddo
!
      call dcopy (n,gf,1,kkt,1)
      if (ni.gt.0) then
        call dgemv ('t',ni,n,1.0d0,gc,ni,xlam,1,1.d0,kkt,1)
      endif
!
      if (subsolver.ge.25.and.subsolver.le.28) then
        do i=1,n
          kkt(i)=kkt(i)-z(i)
          xkkt=xkkt+kkt(i)**2
          ifree=ifree+1
        enddo
      else  
        do i=1,n
          if     (x(i).le.x_l0(i)+0.d-6) then
!
          elseif (x(i).ge.x_u0(i)-0.d-6) then
!
          else
            xkkt=xkkt+kkt(i)**2
            ifree=ifree+1
          endif
        enddo
      endif
!
      xkkt=dsqrt(xkkt)
!
      return
      end subroutine form_kkt
!----------------------------------------------------------------------!
      double precision function drandu1(iseed)
      implicit none
      integer iseed
c     **********
c
c     function rand
c
c     Rand is the portable random number generator of L. Schrage.
c
c     The generator is full cycle, that is, every integer from
c     1 to 2**31 - 2 is generated exactly once in the cycle.
c     It is completely described in TOMS 5(1979),132-138.
c
c     The function statement is
c
c       real function rand(iseed)
c
c     where
c
c       iseed is a positive integer variable which specifies
c         the seed to the random number generator. Given the
c         input seed, rand returns a random number in the
c         open interval (0,1). On output the seed is updated.
c
c     Argonne National Laboratory. MINPACK Project. March 1981.
c
!
!     temporarily hacked to produce a double ...
!
c     **********
      integer a,b15,b16,fhi,k,leftlo,p,xhi,xalo
      real c,rand
c
c     Set a = 7**5, b15 = 2**15, b16 = 2**16, p = 2**31-1, c = 1/p.
c
      data a/16807/, b15/32768/, b16/65536/, p/2147483647/
      data c/4.656612875e-10/
c
c     There are 8 steps in rand.
c
c     1. Get 15 hi order bits of iseed.
c     2. Get 16 lo bits of iseed and form lo product.
c     3. Get 15 hi order bits of lo product.
c     4. Form the 31 highest bits of full product.
c     5. Get overflo past 31st bit of full product.
c     6. Assemble all the parts and pre-substract p.
c        The parentheses are essential.
c     7. Add p back if necessary.
c     8. Multiply by 1/(2**31-1).
c
      xhi = iseed/b16
      xalo = (iseed - xhi*b16)*a
      leftlo = xalo/b16
      fhi = xhi*a + leftlo
      k = fhi/b15
      iseed = (((xalo-leftlo*b16)-p) + (fhi-k*b15)*b16) + k
      if (iseed .lt. 0) iseed = iseed + p
      rand = c*float(iseed)
!
      drandu1=dble(rand)
!
      return
c
c     Last card of function rand. ! drandu1 !
c
      end
!----------------------------------------------------------------------!
      subroutine sort (a,k,n)
      implicit         none
      integer          nl1,n,i,ip1,j,k
      double precision a(k,n)
!
! presumably from numerical methods in fortran
!
c
c --- sort real array a from large to small --------------------------
c
      nl1=n-1
      do 20 i=1,nl1
        ip1=i+1
        do 10 j=ip1,n
          if (a(1,i).lt.a(1,j)) then
            call swap(a(1,i),a(1,j))
            call swap(a(2,i),a(2,j))
          endif
 10     continue
 20   continue
      return
      end
!-----------------------------------------------------------------------
      subroutine bubsort3a (a,k,n)
      implicit         none
      integer          i,j,k,n
      double precision a(k,n)
      logical          sort
!
! presumably from numerical methods in fortran
!
c
c --- bubble sort real matrix a from large to small w.r.t. a(1,*)
c
      sort=.true.
 10   if (sort) then
        sort=.false.
        do 20 i=1,n-1
           if (a(1,i).lt.a(1,i+1)) then
              call swap(a(1,i),a(1,i+1))
c ........... swop columns 2 through m .................................
              do 15 j=2,2
                 call swap(a(j,i),a(j,i+1))
 15           continue
c ........... end swop matrix a ........................................
              sort=.true.
           endif
 20     continue
        goto 10
      endif
!
      return
      end
!-----------------------------------------------------------------------
      subroutine swap (arg1,arg2)
      implicit         none
      double precision arg1,arg2,temp
!
      temp=arg1
      arg1=arg2
      arg2=temp
!
      return
      end
!-----------------------------------------------------------------------
      subroutine iswap (iarg1,iarg2)
      implicit         none
      integer          iarg1,iarg2,itemp
!
      itemp=iarg1
      iarg1=iarg2
      iarg2=itemp
!
      return
      end
!-----------------------------------------------------------------------
      SUBROUTINE DGERX1(M,N,ALPHA,X,INCX,Y,INCY,A,LDA,BTOL1)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BTOL1
      INTEGER INCX,INCY,LDA,M,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(LDA,*)
*     ..
*
*  Purpose
*  =======
*
*  DGERX1   performs the operation
*
*     A(J,I) := max(BTOL1,ALPHA/X(I)*|Y(J,I)|)
*
*  where alpha is a scalar, x is an n element vector, y is an m by n
*  matrix and A is an m by n matrix.
*
*  Based on the level 2 Blas routine DGER.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JY,KX
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGERX1  ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (INCX.NE.1) STOP ' INCX.NE.1 IN DGERX1'
      IF (INCY.NE.1) STOP ' INCY.NE.1 IN DGERX1'
*
      DO 20 J = 1,N
          DO 10 I = 1,M
              !IF (Y(I,J).NE.ZERO) THEN
                  TEMP = ALPHA*DABS(Y(I,J))
                  A(I,J) = max(BTOL1,TEMP/X(J))
              !END IF
   10     CONTINUE
   20 CONTINUE
*
      RETURN
*
*     End of DGERX1  .
*
      END
!--------------------------------------------------------------------------
       SUBROUTINE DGERX2(M,N,ALPHA,X,INCX,Y,INCY,A,LDA,BP,BTOL1)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BTOL1
      INTEGER INCX,INCY,LDA,M,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(LDA,*),BP(LDA,*)
*     ..
*
*  Purpose
*  =======
*
*  DGERX2   performs the operation
*
*     A(J,I) := max(btol1,-(BP(J,I)-1.D0)/X(I)*|Y(J,I)|)
*
*  where alpha (unused) is a scalar, x is an n element vector, y is an m by n
*  matrix, A is an m by n matrix and BP is an m by n matrix.
*
*
*  Based on the level 2 Blas routine DGER.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JY,KX
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGERX2  ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (INCX.NE.1) STOP ' INCX.NE.1 IN DGERX2'
      IF (INCY.NE.1) STOP ' INCY.NE.1 IN DGERX2'
*
      DO 20 J = 1,N
          DO 10 I = 1,M
              !IF (Y(I,J).NE.ZERO) THEN
                  TEMP = -(bp(i,j)-1.d0)*DABS(Y(I,J))
                  A(I,J) = max(BTOL1,TEMP/X(J))
              !END IF
   10     CONTINUE
   20 CONTINUE
*
      RETURN
*
*     End of DGERX2  .
*
      END
!--------------------------------------------------------------------------
       SUBROUTINE DGERX3(M,N,ALPHA,X,INCX,Y,INCY,A,LDA,BP)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BTOL1
      INTEGER INCX,INCY,LDA,M,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(LDA,*),BP(LDA,*)
*     ..
*
*  Purpose
*  =======
*
*  DGERX3   performs the operation
*
*     A(J,I) := (BP(J,I)-1.D0)/X(I)*Y(J,I)
*
*  where alpha (unused) is a scalar, x is an n element vector, y is an m by n
*  matrix, A is an m by n matrix and BP is an m by n matrix.
*
*
*  Based on the level 2 Blas routine DGER.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JY,KX
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGERX2  ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (INCX.NE.1) STOP ' INCX.NE.1 IN DGERX3'
      IF (INCY.NE.1) STOP ' INCY.NE.1 IN DGERX3'
*
      DO 20 J = 1,N
          DO 10 I = 1,M
              !IF (Y(I,J).NE.ZERO) THEN
                  TEMP = (bp(i,j)-1.d0)*Y(I,J)
                  A(I,J) = TEMP/X(J)
              !END IF
   10     CONTINUE
   20 CONTINUE
*
      RETURN
*
*     End of DGERX3  .
*
      END 
!--------------------------------------------------------------------------
      SUBROUTINE DGERX1_NC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA,BTOL1)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BTOL1
      INTEGER INCX,INCY,LDA,M,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(LDA,*)
*     ..
*
*  Purpose
*  =======
*
*  DGERX1_NC   performs the operation
*
*     A(J,I) := max(BTOL1,-ALPHA/X(I)*Y(J,I))
*
*  where alpha is a scalar, x is an n element vector, y is an m by n
*  matrix and A is an m by n matrix.
*
*  Based on the level 2 Blas routine DGER.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JY,KX
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGERX1  ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (INCX.NE.1) STOP ' INCX.NE.1 IN DGERX1_NC'
      IF (INCY.NE.1) STOP ' INCY.NE.1 IN DGERX1_NC'
*
      DO 20 J = 1,N
          DO 10 I = 1,M
              !IF (Y(I,J).NE.ZERO) THEN
                  TEMP = -ALPHA*Y(I,J)
                  A(I,J) = max(BTOL1,TEMP/X(J))
              !END IF
   10     CONTINUE
   20 CONTINUE
*
      RETURN
*
*     End of DGERX1_NC  .
*
      END
!----------------------------------------------------------------------!
      subroutine FinDiffGradHessDense (n,x,m,ni,ne,f,c,ictrl,lctrl,
     &           rctrl,cctrl,do_diff,iuser,luser,cuser,ruser,eqn,lin,
     &           sloop)      
!----------------------------------------------------------------------!
!                                                                      !
!  Compute the Jacobian and Hessian using differences. Coded for       !
!  clarity, not efficiency.                                            !  
!  Please see the users manual for type declarations and comments      !
!                                                                      !
!----------------------------------------------------------------------!
      implicit          none
      include           'ctrl.h'
      logical           eqn(*), lin(*)
      logical           do_diff
      logical           debug_local,form_grad_here,form_hess_here
      integer           i,j,k,l,m,n,ni,ne,nfe,sloop
      double precision  x(n),dx1(n),dx2(n),temp
      double precision  dx12(n),dx21(n)
      double precision  c_x(m),f_x
      double precision  c_d1(m),f_d1
      double precision  c_d12(m),f_d12,c_d21(m),f_d21
      double precision  c_d2(m),f_d2
      double precision  f,c(m),gf(n),gc(m,n)
      double precision  gf2(n,n),gc2(m,n,n) 
!
      data              debug_local    /.false./
      data              form_grad_here /.false./
      data              form_hess_here /.true./
!      
      include          'ctrl_get.inc'                  
!
      if (sloop.ne.1) return ! only applicable for mode extraction problems
!
      do_diff=.true.
!
      !if (structure.ne.1) stop ' structure.ne.1 in finDiffGradHessDense'
!
      call SAOi_funcs(n,m,ni,ne,x,f_x,c_x,iuser,luser,
     &                cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
!  Construct the Jacobian using forward differences   
!
      if (form_grad_here) then
!
        nfe = 0          
        dx1 = x
!
        do i=1,n
          dx1(i)=x(i)+deltx
          call SAOi_funcs(n,m,ni,ne,dx1,f_d1,c_d1,iuser,luser,
     &                  cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
          nfe=nfe+1
          gf(i)=(f_d1-f_x)/deltx
!        
          if (debug_local) write(35,1001)i,gf(i) ! same format as in CheckSAOi_grad.out
!                        
          do j=1,m
            k=(i-1)*ni+j                     ! for debugging only
            gc(j,i)=(c_d1(j)-c_x(j))/deltx
!
            if (debug_local) write(35,1002)i,j,k,gc(j,i) ! same format as in CheckSAOi_grad.out 
!
          end do
          dx1(i)=x(i)
        end do  
!
      endif
!
!  Construct the diagonal (lower triangle) Hessian terms using forward differences   
!
      if (form_hess_here) then
!
        !deltxH = 1.d-4  
!
        gf2 = 0.d0
        gc2 = 0.d0
!
        dx1 = x
        dx2 = x
        do i=1,n
!
          dx1(i)=x(i)+deltxH
          dx2(i)=x(i)+2.d0*deltxH        
!        
          call SAOi_funcs(n,m,ni,ne,dx1,f_d1,c_d1,iuser,luser,
     &                  cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
          call SAOi_funcs(n,m,ni,ne,dx2,f_d2,c_d2,iuser,luser,
     &                  cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
          nfe=nfe+2
          gf2(i,i)=(f_d2-2.d0*f_d1+f_x)/deltxH**2
!        
          do j=1,m
            gc2(j,i,i)=(c_d2(j)-2.d0*c_d1(j)+c_x(j))/deltxH**2
          end do
!
          dx1(i)=x(i)
          dx2(i)=x(i)
!
        end do
!
!  Construct the off-diagonal lower triangle Hessian terms using forward differences   
!
        k = 0
        dx12 = x
        dx21 = x
        dx2 = x
        do i=1,n
          do l=i+1,n
!     
!  f{(x(i) + h)(x(l) + h)}
!
            dx2(i)=x(i)+deltxH
            dx2(l)=x(l)+deltxH
!
            call SAOi_funcs(n,m,ni,ne,dx2,f_d2,c_d2,iuser,luser,
     &                      cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!     
            dx2(i)=x(i)
            dx2(l)=x(l)
!
!  f{(x(i) + h)(x(l)}
!     
            dx12(i)=x(i)+deltxH  
! 
            call SAOi_funcs(n,m,ni,ne,dx12,f_d12,c_d12,iuser,luser,
     &                      cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
            dx12(i)=x(i)
!
!  f{(x(i))(x(l) + h}
!
            dx21(l)=x(l)+deltxH  

            call SAOi_funcs(n,m,ni,ne,dx21,f_d21,c_d21,iuser,luser,
     &                    cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
            dx21(l)=x(l)   
!
            nfe=nfe+3
!
!  Now construct the lower triangle
!          
            gf2(l,i)=(f_d2-f_d12-f_d21+f_x)/deltxH**2
!        
            do j=1,m
              gc2(j,l,i)=(c_d2(j)-c_d12(j)-c_d21(j)+c_x(j))/deltxH**2
            end do
!          
          enddo
        end do
!
        if (debug_local) then 
          write (36,*) ' '
          write (36,*) ' Lower Triangle of the Objective'
          write (36,*) ' '
          do i=1,n
            write(36,1000) (nint(gf2(i,l)),l=1,i)
          enddo
!
          do j=1,m 
            write (36,*) ' '
            write (36,*) ' Lower Triangle of Constraint no ',j
            write (36,*) ' '
            do i=1,n
              write(36,1000) (nint(gc2(j,i,l)),l=1,i)
            enddo
          enddo  
          write (36,*) ' '
        endif
!
        do i=1,n
          do l=1,i
            k=k+1
            ruser(k) = gf2(i,l)
          enddo
        enddo
!
        do j=1,m 
          do i=1,n
            do l=1,i
              k=k+1
              ruser(k) = gc2(j,i,l)
            enddo
          enddo
        enddo
!
        if (k.gt.rmax) then
          write(*,*) 'k.gt.rmax in FinDiffGradHessDense;',
     &               ' rmax is defined in size.h',k,rmax
          stop
        endif  
!
      endif
!      
      do_diff=.false.
!
      return
      !stop
!
 1000 format(1000i4)
 1001 format(1i8,16x,1e22.10)
 1002 format(3i8,1e22.10)
!   
      end subroutine FinDiffGradHessDense
!--------------------------------------------------------------------------
      subroutine FinDiffGradHessDense1 (n,x,m,ni,ne,f,c,ictrl,lctrl,
     &           rctrl,cctrl,do_diff,iuser,luser,cuser,ruser,eqn,lin,
     &           sloop,acurv,nnzh,xlam,H_ne,H_val,H_row,H_col,H_ptr)      
!----------------------------------------------------------------------!
!                                                                      !
!  Compute the Jacobian and Hessian using differences. Coded for       !
!  clarity, not efficiency.                                            !  
!  Please see the users manual for type declarations and comments      !
!                                                                      !
!----------------------------------------------------------------------!
      implicit          none
      include           'ctrl.h'
      logical           eqn(*), lin(*)
      logical           do_diff
      logical           debug_local
      integer           i,j,k,l,m,n,ni,ne,nfe,sloop
      double precision  x(n),dx1(n),dx2(n),temp
      double precision  dx12(n),dx21(n)
      double precision  c_x(m),f_x
      double precision  c_d1(m),f_d1
      double precision  c_d12(m),f_d12,c_d21(m),f_d21
      double precision  c_d2(m),f_d2
      double precision  f,c(m),gf(n),gc(m,n)
      double precision  gf2(n,n),gc2(m,n,n),acurv
!
      
      integer           H_row(nnzh),H_col(nnzh),H_ptr(n+1)
      integer           H_ne,nnzh,k1,l1
      double precision  H_val(nnzh),xlam(m),xlamj
!      
      data              debug_local /.true./
!      
      include          'ctrl_get.inc'                  
!
      if (sloop.eq.1) then
!
        do_diff=.true.
!
        !if (structure.ne.1) stop ' structure.ne.1 in finDiffGradHessDense'
!
        call SAOi_funcs(n,m,ni,ne,x,f_x,c_x,iuser,luser,
     &                cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
!  Construct the Jacobian using forward differences   
!
        nfe = 0          
        dx1 = x
!
        do i=1,n
          dx1(i)=x(i)+deltx
          call SAOi_funcs(n,m,ni,ne,dx1,f_d1,c_d1,iuser,luser,
     &                    cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
          nfe=nfe+1
          gf(i)=(f_d1-f_x)/deltx
!          
          if (debug_local) write(35,1001)i,gf(i) ! same format as in CheckSAOi_grad.out
!                          
          do j=1,m
            k=(i-1)*ni+j                     ! for debugging only
            gc(j,i)=(c_d1(j)-c_x(j))/deltx
!
            if (debug_local) write(35,1002)i,j,k,gc(j,i) ! same format as in CheckSAOi_grad.out 
!
          end do
          dx1(i)=x(i)
        end do        
!
!  Construct the diagonal (lower triangle) Hessian terms using forward differences   
!
        !deltxH = 1.d-4  
!
        gf2 = 0.d0
        gc2 = 0.d0
!
        dx1 = x
        dx2 = x
        do i=1,n
!
          dx1(i)=x(i)+deltxH
          dx2(i)=x(i)+2.d0*deltxH        
!          
          call SAOi_funcs(n,m,ni,ne,dx1,f_d1,c_d1,iuser,luser,
     &                    cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
          call SAOi_funcs(n,m,ni,ne,dx2,f_d2,c_d2,iuser,luser,
     &                    cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
          nfe=nfe+2
          gf2(i,i)=(f_d2-2.d0*f_d1+f_x)/deltxH**2
!          
          do j=1,m
            gc2(j,i,i)=(c_d2(j)-2.d0*c_d1(j)+c_x(j))/deltxH**2
          end do
!
          dx1(i)=x(i)
          dx2(i)=x(i)
!
        end do
!
!    Construct the off-diagonal lower triangle Hessian terms using forward differences   
!
        k = 0
        dx12 = x
        dx21 = x
        dx2 = x
        do i=1,n
          do l=i+1,n
!       
!    f{(x(i) + h)(x(l) + h)}
!
            dx2(i)=x(i)+deltxH
            dx2(l)=x(l)+deltxH
!
            call SAOi_funcs(n,m,ni,ne,dx2,f_d2,c_d2,iuser,luser,
     &                      cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!       
            dx2(i)=x(i)
            dx2(l)=x(l)
!
!    f{(x(i) + h)(x(l)}
!       
            dx12(i)=x(i)+deltxH  
! 
            call SAOi_funcs(n,m,ni,ne,dx12,f_d12,c_d12,iuser,luser,
     &                      cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
            dx12(i)=x(i)
!
!    f{(x(i))(x(l) + h}
!
            dx21(l)=x(l)+deltxH  

            call SAOi_funcs(n,m,ni,ne,dx21,f_d21,c_d21,iuser,luser,
     &                      cuser,ruser,eqn,lin,ictrl,lctrl,rctrl,cctrl)
!
            dx21(l)=x(l)   
!
            nfe=nfe+3
!
!  Now construct the lower triangle
!            
            gf2(l,i)=(f_d2-f_d12-f_d21+f_x)/deltxH**2
!          
            do j=1,m
              gc2(j,l,i)=(c_d2(j)-c_d12(j)-c_d21(j)+c_x(j))/deltxH**2
            end do
!            
          enddo
        end do
!
        if (debug_local) then 
          write (36,*) ' '
          write (36,*) ' Lower Triangle of the Objective'
          write (36,*) ' '
          do i=1,n
            write(36,1000) (nint(gf2(i,l)),l=1,i)
          enddo
!
          do j=1,m 
            write (36,*) ' '
            write (36,*) ' Lower Triangle of Constraint no ',j
            write (36,*) ' '
            do i=1,n
              write(36,1000) (nint(gc2(j,i,l)),l=1,i)
            enddo
          enddo  
          write (36,*) ' '
        endif
!
        do i=1,n
          do l=1,i
            k=k+1
            ruser(k) = gf2(i,l)
          enddo
        enddo
!
        do j=1,m 
          do i=1,n
            do l=1,i
              k=k+1
              ruser(k) = gc2(j,i,l)
            enddo
          enddo
        enddo
!
        if (k.gt.rmax) then
          write(*,*) 'k.gt.rmax in FinDiffGradHessDense;',
     &               ' rmax is defined in size.h',k,rmax
          stop
        endif  
!
        do_diff=.false.
!        
      endif
!
      H_val = 0.d0

      H_ne = nnzh
      k = 0
      k1 = 0
      do i=1,n
        do l1=1,i
          k=k+1
          k1=k1+1
          H_val(k1) = ruser(k)
          !if (i.eq.l1) H_val(k1) = acurv
          H_row(k1) = i
          H_col(k1) = l1
        enddo
      enddo
!
      do j=1,m 
        k1 = 0
        xlamj = xlam(j)      
        do i=1,n
          do l1=1,i
            k=k+1
            k1=k1+1
            H_val(k1) = H_val(k1) + xlamj*ruser(k)
          enddo
        enddo
      enddo
!
      !write(*,*)xlam,H_row,H_col,H_val
!
      return
!
 1000 format(1000i4)
 1001 format(1i8,16x,1e22.10)
 1002 format(3i8,1e22.10)
!   
      end subroutine FinDiffGradHessDense1
!--------------------------------------------------------------------------
      subroutine FinDiffGradHessDense0 (n,x,m,ni,ne,f,c,ictrl,lctrl,
     &           rctrl,cctrl,do_diff,iuser,luser,cuser,ruser,eqn,lin,
     &           sloop,acurv,nnzh,xlam,H_ne,H_val,H_row,H_col,H_ptr)      
!----------------------------------------------------------------------!
!                                                                      !
!  Compute the Jacobian and Hessian using differences. Coded for       !
!  clarity, not efficiency.                                            !  
!  Please see the users manual for type declarations and comments      !
!                                                                      !
!----------------------------------------------------------------------!
      implicit          none
      include           'ctrl.h'
      logical           eqn(*), lin(*)
      logical           do_diff
      logical           debug_local
      integer           i,j,k,l,m,n,ni,ne,nfe,sloop
      double precision  x(n),temp
      double precision  f,c(m),gf(n),gc(m,n)
      double precision  acurv
!
      integer           H_row(nnzh),H_col(nnzh),H_ptr(n+1)
      integer           H_ne,nnzh,k1,l1
      double precision  H_val(nnzh),xlam(m),xlamj
!      
      data              debug_local /.true./
!      
      include          'ctrl_get.inc'                  
!
!
      H_val = 0.d0

      H_ne = nnzh
      k = 0
      k1 = 0
      do i=1,n
        do l1=1,i
          k=k+1
          k1=k1+1
          H_val(k1) = ruser(k)
          !if (i.eq.l1) H_val(k1) = acurv
          H_row(k1) = i
          H_col(k1) = l1
        enddo
      enddo
!
      do j=1,m 
        k1 = 0
        xlamj = xlam(j)      
        do i=1,n
          do l1=1,i
            k=k+1
            k1=k1+1
            H_val(k1) = H_val(k1) + xlamj*ruser(k)
          enddo
        enddo
      enddo
!
      !write(*,*)xlam,H_row,H_col,H_val
!
      return
!
 1000 format(1000i4)
 1001 format(1i8,16x,1e22.10)
 1002 format(3i8,1e22.10)
!   
      end subroutine FinDiffGradHessDense0
!--------------------------------------------------------------------------

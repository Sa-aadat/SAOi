! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: objective and constraint functions gradients
! SAOi:


      subroutine SAOi_grads (n, m, ni, ne, x, gf, gc, nnz, Acol, Aptr,
     &                       iuser, luser, cuser, ruser, eqn, lin,
     &                       ictrl, lctrl, rctrl, cctrl)
!----------------------------------------------------------------------!
!                                                                      !
!  Compute the gradients gf(i) of the objective function f, and the    !
!  derivatives gc(k) of the inequality constraint functions c          !
!  w.r.t. the variables x(i)                                           !
!                                                                      !
!  Please see the users manual for type declarations and comments      !
!                                                                      !
!                                                                      !
!  The default storage scheme is the dense storage scheme. In          !
!      this scheme, it is only required to specify the vector          !
!      gf and the matrix gc. The dimensions of gf and gc are           !
!      gf(n) and gc(ni,n) respectively. For details, please            !
!      see the users manual                                            !
!                                                                      !
!  The default storage scheme if the sparse implementation is          !
!      selected is the compressed sparse row (CSR) storage scheme.     !
!      The dimensions of Acol and Aptr in the CSR storage sheme        !
!      are Acol(nnz) and Aptr(ni+1) respectively. For details,         !
!      please see the users manual                                     !
!                                                                      !
!  nnz is the number of non-zero entries in gc. If the standard        !
!      sparse storage structure is used, nnz should be declared        !
!      in Initialize.f, and its value should not be changed            !
!                                                                      !
!  Acol and Aptr are to be declared here if the sparse storage         !
!      structure is used. For details, please see the users manual     !
!                                                                      !
!  iuser, luser, cuser and ruser are user arrays, which may be used    !
!      at will to pass arbitrary data around between the user          !
!      routines                                                        !
!                                                                      !
!----------------------------------------------------------------------!
      implicit         none
      include          'ctrl.h'
      logical          eqn(*), lin(*)
      integer          i, j, n, ni, ne, m, nnz, Acol(nnz), Aptr(ni+1)
      double precision x(n), gf(n), gc(nnz), gfold(n)
      double precision tmp0, tmp1, tmp2, tmp3, tmp4
!
      logical          loaded_knee, dooo
      integer          nelx, nely, istrgrad, iadjoint, mcompliance
      integer          ind, icount, loop, k, l, ii, jj, jni
      integer          isdense, jacobfilt, jblockwrite
      integer          lentot, actemp(n), aptrstr(n+1)
      double precision penal, volfrac, rmin, p1, p2, penalq, ds
      double precision rho_min, sall, conlim, gcvol(n), theta
      double precision eps, eta, penal_bruns, penal2, glo, ghi
      double precision gradlim, divval, gctemp1, gctemp(n), fscale
!
      include          'ctrl_get.inc'
!
      write(*,*) ' Entering G '
!
      rho_min        = ruser(1)
      penal          = ruser(2)
      volfrac        = ruser(3)
      rmin           = ruser(4)
      p1             = ruser(5)
      p2             = ruser(6)
      penalq         = ruser(7)
      penal_bruns    = ruser(8)
      penal2         = ruser(9)
      sall           = ruser(10)
      conlim         = ruser(11)
      gradlim        = ruser(12)
      divval         = ruser(13)
      theta          = ruser(14)
      eps            = ruser(15)  
      fscale         = ruser(16)
      nelx           = iuser(1)
      nely           = iuser(2)
      istrgrad       = iuser(3)
      iadjoint       = iuser(4)
      isdense        = iuser(5)
      jacobfilt      = iuser(6)
      jblockwrite    = iuser(7)
      mcompliance    = iuser(8) 
      loaded_knee    = luser(1)
!
! ----------------------------------------------------------------------
! --------------------- Get sensitivities of objective function --------
! ----------------------------------------------------------------------    
!
      lentot = 0
      if (structure.ge.3) then 
         Aptr(1) = 1
         k = 2
      endif
!
      if (istrgrad.ne.0) then 
!     
         OPEN (41,FILE='Sgrads.str',STATUS='old',
     &        FORM='UNFORMATTED',ERR=1009)
         rewind(41)
         OPEN (43,FILE='Jacobindex.ind',STATUS='old',
     &        FORM='UNFORMATTED',ERR=1009)
         rewind(43)
         OPEN (44,FILE='Jacobpointer.ind',STATUS='old',
     &        FORM='UNFORMATTED',ERR=1009)
         rewind(44)
         read(44) aptrstr
         close(44)
!     
         lentot = aptrstr(n+1) - 1
         gctemp1 = 0.0d0
         if (structure.ge.3) then
            do ii=2,n+1
               Aptr(ii) = aptrstr(ii)
            enddo
            k = n + 2
         endif
!     
         if (jblockwrite.ne.1) then
! Read/write sequentially
            do i=1,lentot
               if (structure.ge.3) read(43) acol(i)
               read(41) gc(i)
               if (dabs(gc(i)).gt.gctemp1) gctemp1 = dabs(gc(i))
            enddo
         else 
! Block read/write
            dooo = .true.
            icount = 0
            jj = 0
            do while (dooo)
               read(43) actemp
               read(41) gctemp
               do ii = 1,n
                  jj = jj + 1
                  if (structure.ge.3) acol(jj) = actemp(ii)
                  gc(jj) = gctemp(ii)
                  if (dabs(gc(jj)).gt.gctemp1) gctemp1 = dabs(gc(jj))
               enddo
               icount = icount + n
               if ((lentot - icount).lt.n) then 
                  dooo = .false.
               endif
            enddo
            if ((lentot - icount).gt.0) then
               read(43) actemp  
               read(41) gctemp
               jj = 0
               do ii = (icount+1),lentot
                  jj = jj + 1
                  if (actemp(jj).eq.(-9999)) 
     &                 stop 'Indexing error 1 in Gradients.f'
                  if (structure.ge.3) acol(ii) = actemp(jj)
                  gc(ii) = gctemp(jj)
                  if (dabs(gc(ii)).gt.gctemp1) gctemp1 = dabs(gc(ii))
               enddo
               if (actemp(jj+1).ne.(-9999))
     &              stop 'Indexing error 2 in Gradients.f'
            endif
         endif
!
         if (jacobfilt.eq.1) then
            gradlim = gctemp1/divval
            ruser(12) = gradlim
         endif
!               
         close(43)
         close(41)
!          
         if (.not.strict_struct) then
            jni = 0
            do i=1,n
               if (aptrstr(i+1)-aptrstr(i).eq.0) then
                  jni=jni+1
               endif
            enddo
         endif
!
      endif
!
      if (mcompliance.eq.1) then
         OPEN (24,FILE='TowerDervs',STATUS='OLD',
     &        FORM='UNFORMATTED',ERR=1009)
         rewind(24)
         read(24) gf
         read(24) gcvol
         close(24)
         do i=1,n
            gc(i + lentot) = gcvol(i)
            if (structure.ge.3) acol(i + lentot) = i
         enddo
         if (structure.ge.3) Aptr(k) = Aptr(k-1) + n
      else
         do i=1,n
            gfold(i) = fscale*penal2*x(i)**(penal2-1.d0) ! For weight min problem
         enddo
         if (rmin.gt.1.d0) then
           write(*,*) rmin 
           call filter1(nelx,nely,rmin,x,gf,gfold,n)
         else
           do i=1,n
             gf(i)=gfold(i)
           enddo
        endif
      endif
      call system('rm -f TowerDervs')
!
      if (structure.le.2) call convert_gcstruc(n, ni, nnz, gc)
!
! ---------------- Debug stuff -----------------------------------------  
!       call system('rm -f Dgrads.chk')
!       open(file='Dgrads.chk',unit=49,status='unknown')
!       icount = 0
!       do i=1,n
!          loop = aptr(i+1) - aptr(i)
!          do l = 1,loop
!             icount = icount + 1
!             ind = acol(icount)
!             ds = gc(icount)
!             if (iadjoint.eq.1) then
!                ii = i
!                jj = ind
!             else
!                ii = ind
!                jj = i
!             endif
!             k=(jj-1)*n + ii
!             write(49,1245) k,ds
!          enddo
!       enddo
!       close(49)
!       stop 'Gradients.f'
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! ------- No checks on gradients yet (like in finite difference) ------
! ---------------------------------------------------------------------
!
      write(*,*) penal,eps! ,Ruser(13),Ruser(14)
!
      glo= 1.d16
      ghi=-1.d16
      do j=1,ni*n
        glo=min(glo,gc(j))
        ghi=max(ghi,gc(j))
      enddo
!
      glo= 1.d16
      ghi=-1.d16
      do j=1,ni*n
        if (gc(j).ne.0.d0) glo=min(glo,dabs(gc(j)))
                           ghi=max(ghi,dabs(gc(j)))
      enddo
!
      write(23,*) ' Leaving G - |lo| and |hi| gradients are',glo,ghi
!      
      return
!      
 1009 stop ' fatal i/o error in Gradients.f. period. '
 1200 format (/,' p = ',f5.3,'; lo = ',f7.5,'; sparsity = ',f5.3,
     &        '; selected = ',f5.3,/)
 1245 format (i7,d24.15)
!
      end subroutine SAOi_grads
!----------------------------------------------------------------------!
      subroutine convert_gcstruc(n, ni, nnz, gc)
      integer          i, j, jj, k, n, ni
      double precision gc(nnz), gctemp(nnz)
!            
      jj = 0
      do i=1,n
         do j=1,ni
            jj = jj + 1
            k=(j-1)*ni + i
            gctemp(k) = 0.d0
            gctemp(k) = gc(jj)
         enddo
      enddo
!
      do i=1,nnz
         gc(i) = 0.d0
         gc(i) = gctemp(i)
      enddo
!      
      return
!
      end subroutine convert_gcstruc
!----------------------------------------------------------------------!

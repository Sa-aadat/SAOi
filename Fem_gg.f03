subroutine Fem_gg(x,n,ni,xlam,nnzh,Hval,Hptr,Hcol,Hne)
!-------------------------------------------------------------------------
! Wed 11 Oct 2017 07:51:45 SAST
! Modified to be sparse
! Calculating the exact curvature information for SAND constraints
! Called from cplexqp.f
! Adapted from Program 5.1 in Programming the Finite Element Method-Smith 
!-------------------------------------------------------------------------
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15),rik=selected_int_kind(12)
 INTEGER::i,iel,k,ndim=2,ndof,j,nlen,nodof=2,knels,layerv
 INTEGER::nels,nn,nod,neq,np_types,hels,hcount,Hrow(nnzh)
 INTEGER,INTENT(IN)::ni,n,nnzh
 REAL(iwp),INTENT(IN)::x(n),xlam(ni)
 REAL(iwp),INTENT(OUT)::Hval(nnzh)
 INTEGER,INTENT(OUT)::Hcol(nnzh),Hptr(n+1),Hne
 REAL(iwp)::zero=0.0_iwp,xlamj(0:ni)
 INTEGER(rik)::ij
!-----------------------dynamic arrays------------------------------------
 INTEGER,ALLOCATABLE::g(:),g_g(:,:),cntrnels(:),cntrneq(:),cntrneq2(:)
 INTEGER(rik),ALLOCATABLE::kdiag(:)
 REAL(iwp),ALLOCATABLE::km(:,:),disp(:),loads(:),kv(:),g_c(:)

!-----------------------input and initialisation--------------------------
 rewind(102)
 read(102) nels,nod,neq,layerv,nn,knels
 Hne=nels+knels
 Hval=zero
 xlamj(0)=zero
 xlamj(1:ni)=xlam
 ndof=nod*nodof
 ALLOCATE(g(ndof),g_c(0:neq),km(ndof,ndof),g_g(ndof,nels),disp(0:neq))
 ALLOCATE(kdiag(neq))
 disp(0) = zero
 do i = 1,neq
  disp(i)=x(nels+i)
 enddo
 rewind(101)
 read(101) g_g,kdiag
!  print*,g_g(6,6)
!  print*,g_g(7,6)
 ALLOCATE(kv(kdiag(neq)),cntrnels(nels),cntrneq(knels),cntrneq2(knels))
 rewind(104)
 read(104) cntrnels,cntrneq,cntrneq2,Hcol,Hptr
!-----------------------constraint gradients----------------
 kv=zero
 ij=kdiag(neq)
 do i = 1,neq
  disp(i)=x(nels+i)
 enddo
 k=1
 elements_1: DO iel=1,nels
   do i = 1,ndof
    g(i)=g_g(i,iel)
   enddo
   km=zero
   SELECT CASE(nod)
    CASE(4)
        CALL Q4_ortho_dd(km,x(iel))
!     CASE(8)
!         CALL Q4_ortho_dd(km,x(i))
    CASE DEFAULT
        WRITE(*,'(A)') "# Element type not recognised"
   END SELECT
   g_c(g) = MATMUL(km,disp(g))!works because of first few elements
   Hval(cntrnels(iel))=dot_product(g_c(g),xlamj(g))!for cplex optimiser*2.d0
   km=zero
    SELECT CASE(nod)
    CASE(4)
        CALL Q4_ortho_d(km,x(iel))
!     CASE(8)
!         CALL Q4_ortho_dd(km,x(i))
    CASE DEFAULT
        WRITE(*,'(A)') "# Element type not recognised"
   END SELECT
     g_c(g) = MATMUL(km,xlamj(g))
   do i = 1,ndof
    if ((g(i) /= 0 ).and. (abs(g_c(g(i)))>1.d-9))then!
      Hval(cntrneq(k))=g_c(g(i))
      Hval(cntrneq2(k))=g_c(g(i))
      k=k+1
    endif
   enddo
 ENDDO elements_1
DEALLOCATE(g,g_c,km,g_g,disp,kdiag,kv,cntrnels,cntrneq,cntrneq2)

END subroutine Fem_gg

subroutine Sep_gg(x,n,ni,xlam,w)
!-------------------------------------------------------------------------
! Mon 09 Jan 2017 07:07:11 SAST Sa-aadat Parker
! Tue 07 Feb 2017 09:20:28 SAST Sa-aadat Parker
! Modified to be sparse
! Calculating the separable curvature information for SAND constraints
! Called from cplexqp.f
! Adapted from Program 5.1 in Programming the Finite Element Method-Smith
!-------------------------------------------------------------------------
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15),rik=selected_int_kind(12)
 INTEGER::i,iel,k,ndim=2,ndof,nr,j,nlen,nodof=2
 INTEGER::nels,nn,nod,neq,np_types,hels,hcount
 INTEGER,INTENT(IN)::ni,n
 REAL(iwp),INTENT(IN)::x(n),xlam(ni)
 REAL(iwp),INTENT(OUT)::w(n)
 REAL(iwp)::zero=0.0_iwp,xlamj(0:ni)
!-----------------------dynamic arrays------------------------------------
 INTEGER,ALLOCATABLE::g(:),g_g(:,:)
 INTEGER(rik),ALLOCATABLE::kdiag(:)
 INTEGER(rik)::ij
 REAL(iwp),ALLOCATABLE::km(:,:),disp(:),loads(:),kv(:),g_c(:)

!-----------------------input and initialisation--------------------------
 rewind(102)
 read(102) nels,nod,neq,nn
 xlamj(0)=zero
 xlamj(1:ni)=xlam
 ndof=nod*nodof
 ALLOCATE(g(ndof),g_c(0:neq),km(ndof,ndof),g_g(ndof,nels),disp(0:neq))
 ALLOCATE(kdiag(neq))
 disp(0) = zero
 do i = 1,neq
  disp(i)=x(nels+i)
 enddo
 rewind(101)
 read(101) g_g,kdiag
 ALLOCATE(kv(kdiag(neq)))
!-----------------------constraint gradients----------------
 kv=zero
 ij=kdiag(neq)
 do i = 1,neq
  disp(i)=x(nels+i)
 enddo
 elements_1: DO iel=1,nels
   do i = 1,ndof
    g(i)=g_g(i,iel)
   enddo
   km=zero
    SELECT CASE(nod)
    CASE(4)
        CALL Q4_ortho_dd(km,x(iel))
!     CASE(8)
!         CALL Q4_ortho_dd(km,x(i))
    CASE DEFAULT
        WRITE(*,'(A)') "# Element type not recognised"
   END SELECT
   g_c(g) = MATMUL(km,disp(g))
   w(iel) = dot_product(g_c(g),xlamj(g))
 END DO elements_1

 do i = 1,n
   w(i) = max(5.d-6,(w(i)))!dabs
 enddo
END subroutine Sep_gg

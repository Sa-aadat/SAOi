subroutine FemInitJ(n,x,nnz)
!-------------------------------------------------------------------------
! Wed 01 Feb 2017 07:27:43 SAST Sa-aadat Parker
! This file calculates some invariant gradients in sparse form
! Plane or axisymmetric strain analysis of an elastic solid
! using 3-, 6-, 10- or 15-node right-angled triangles or
! 4-, 8- or 9-node rectangular quadrilaterals. Mesh numbered
! in x(r)- or y(z)- direction.  
! Adapted from Program 5.1 in Programming the Finite Element Method-Smith 
!-------------------------------------------------------------------------
 
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15),rik=selected_int_kind(12)
 INTEGER::iel,ndim=2,ndof,kneq,knels,nodof=2
 INTEGER::nn,nels,nod,neq,layerv
 REAL(iwp)::x(n)
 INTEGER,INTENT(IN)::n,nnz
 REAL(iwp)::zero=0.0_iwp

!-----------------------dynamic arrays------------------------------------
 INTEGER,ALLOCATABLE::etype(:),g(:),g_g(:,:),cntrneq(:),Acolnels(:),      &
  Arowneq(:),Acolsl(:),Arowsl(:),Arow(:),cntrnels(:),Aptr(:),Arowcneq(:)  &
  ,Arownels(:),Acolneq(:),Arowcnels(:),ncount(:),Acol(:),Acolcoo(:)
 INTEGER(rik)::i,nnzf,j,ij,k
 INTEGER(rik),ALLOCATABLE::kdiag(:)
 REAL(iwp),ALLOCATABLE::km(:,:),kv(:),acoo(:),anels(:),aneq(:),xl(:)
 
!-----------------------input and initialisation--------------------------
 rewind(102)
 read(102) nels,nod,neq,layerv,nn,knels,kneq

 ndof=nod*nodof
 ALLOCATE(g(ndof),km(ndof,ndof),g_g(ndof,nels))
 ALLOCATE(cntrnels(knels),Acol(nnz),Aptr(nels+neq+1),kdiag(neq))
 ALLOCATE(cntrneq(kneq),xl(layerv))
 kdiag=0 
!-----------------------loop the elements to find global arrays sizes-----
 ij=ndof*nels
 rewind(101)
 read(101) g_g,kdiag
 ALLOCATE(kv(kdiag(neq)),Acolcoo(ij),Arowcnels(ij),acoo(ij))
!-----------------------element stiffness integration and assembly--------
 kv=zero 
 k=1
 nnzf = kdiag(neq)

 elements_1: DO iel=1,nels
  do i = 1,layerv
    xl(i) = x(iel+nels*(i-1))
  enddo
  g=g_g(:,iel)
  do i = 1,ndof !looking at the column vector related to each angle variable
    if (g(i) /= 0)then!  .and. g_c(g(i))/= 0
!--------equality constraint gradients wrt angle x----
      acoo(k) = k
      Acolcoo(k) = iel
      Arowcnels(k) = g(i)
      k=k+1
    endif
  enddo
  km=zero
  SELECT CASE(nod)
    CASE(4)
        CALL Q4_ortho_n(km,xl,layerv)
    CASE DEFAULT
        WRITE(*,'(A)') "# Element type not recognised"
  END SELECT
  CALL fsparv(kv,km,g,kdiag,ndof,neq,nnzf)
 END DO elements_1
 ALLOCATE(Acolnels(knels),Arownels(neq+1),anels(knels))!
 call coocsr(neq,knels,acoo,Arowcnels,Acolcoo,anels,Acolnels,Arownels)
!--------equality constraint gradients wrt displacement coo format-
!--------this is equal to the global stiffness matrix--------------
 deallocate(Acolcoo,acoo)
 ij = kdiag(neq)*2-neq
 ALLOCATE(Acolcoo(ij),Arowcneq(ij))
 ALLOCATE(acoo(ij))
 acoo = 0
 acoo(1) = 1
 Acolcoo(1) = 1
 Arowcneq(1) = 1
 k=2! start from second row
 DO i = 2,neq
    DO j = kdiag(i-1)+1,kdiag(i)
        if (kv(j) /= 0.d0) then
          acoo(k) = k
          Acolcoo(k) = i
          Arowcneq(k) = i-(kdiag(i)-j)
          k=k+1
          if (i /= i-(kdiag(i)-j)) then
            acoo(k) = k
            Acolcoo(k) = i-(kdiag(i)-j)
            Arowcneq(k) = i
            k=k+1
          endif
        endif
    ENDDO
 ENDDO
 kneq= k-1
 ALLOCATE(Acolneq(kneq),Arowneq(neq+1),aneq(kneq))
 call coocsr(neq,kneq,acoo,Arowcneq,Acolcoo,aneq,Acolneq,Arowneq)
! !place equality constraint gradients wrt thickness x in csr 

 ij = 0

 Aptr = 0
 Aptr(1) = 1

 Acol = 0

 do i=1,neq !travelling from row to row
 !adding row pointers for equality constraints
  j = Arownels(i+1)-Arownels(i)
  Acol(ij+1:ij+j) = Acolnels(Arownels(i):Arownels(i+1)-1)
  ij = ij+j
  k = Arowneq(i+1)-Arowneq(i)
  Acol(ij+1:ij+k) = Acolneq(Arowneq(i):Arowneq(i+1)-1)+nels
  ij= ij+k
  Aptr(i+1) = Aptr(i)+j+k
 enddo
 deallocate(Acolcoo,acoo)
 allocate(Arow(knels))
 
 do i= 1,knels
 Arow(int(anels(i))) = i
 enddo
 
 k=1
 j=0
 do i = 1,neq
  do j = 1,knels
    if (Arowcnels(j) == i) then 
      cntrnels(j) = int(Arow(j)+arowneq(i)-arowneq(1))
    endif
  enddo
 enddo 
 
 deallocate(Arow)
 allocate(Arow(kneq))
 
 do i= 1,kneq
 Arow(int(aneq(i))) = i
 enddo
 
 k=1
 j=0
 do i = 1,neq
  do j = 1,kneq
    if (Arowcneq(j) == i) then 
      cntrneq(j) = int(Arow(j)+arownels(i+1)-arownels(1))
    endif
  enddo
 enddo 

 open(unit=103,FORM='unformatted',file='control.dat')
 write(103) cntrnels,cntrneq,Acol,Aptr
return
END subroutine FemInitJ

subroutine FemInitGG(nnzh)
!-------------------------------------------------------------------------
! Tue 10 Oct 2017 07:31:54 SAST Sa-aadat Parker
! This file calculates some hessian structure in sparse form
! Plane or axisymmetric strain analysis of an elastic solid
! using 3-, 6-, 10- or 15-node right-angled triangles or
! 4-, 8- or 9-node rectangular quadrilaterals. Mesh numbered
! in x(r)- or y(z)- direction.
! Adapted from Program 5.1 in Programming the Finite Element Method-Smith
!-------------------------------------------------------------------------

 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15),rik=selected_int_kind(12)
 INTEGER::i,iel,loaded_nodes,ndim=2,ndof,nstep,npri,kneq,ksl,knels,      &
   nlen,nodof=2,nprops=2,np_types,j,nnzf,ij,k
 INTEGER::n,ni,ne,nn,nels,nod,neq,nnzh
 REAL(iwp)::one=1.0_iwp,penalty=1.0e20_iwp,zero=0.0_iwp,dtim
 CHARACTER(LEN=15)::argv,element,dir,type_2d
 LOGICAL::solid=.TRUE.
!-----------------------dynamic arrays------------------------------------
 INTEGER,ALLOCATABLE::etype(:),g(:),g_g(:,:),g_num(:,:),nf(:,:),Hptr(:),  &
   no(:),node(:),num(:),Arowneq(:),Acolsl(:),Arowsl(:),Acolcoo(:),Arow(:) &
   ,Acolnels(:),Arownels(:),Acolneq(:),Arowcnels(:),Arowcneq(:),ncount(:) &
   ,cntrnels(:),cntrneq(:),cntrneq2(:),Hcol(:)
 INTEGER(rik),ALLOCATABLE::kdiag(:)
 REAL(iwp),ALLOCATABLE::coord(:,:),g_coord(:,:),acoo(:),prop(:,:)         &
   ,anels(:),aneq(:)

!-----------------------input and initialisation--------------------------
 rewind(102)
 read(102) nels,nod,neq,nn,knels!,knels,kneq,cntrnels(knels),cntrneq(kneq)

 ndof=nod*nodof
 ALLOCATE(kdiag(neq),g_g(ndof,nels),g(ndof),Arowcnels(nels),Acolcoo(nels) &
  ,cntrnels(nels),cntrneq(knels),cntrneq2(knels),acoo(nels))

 kdiag=0
!-----------------------loop the elements to find global arrays sizes-----
 rewind(101)
 read(101) g_g,kdiag

!-----------------------element stiffness integration and assembly--------
 k=1
 nnzf = kdiag(neq)

elements_1: DO iel=1,nels
 !------equality constraint hessian wrt element angle x----
  acoo(iel) = iel
  Acolcoo(iel) = iel
  Arowcnels(iel) = iel
 END DO elements_1
 ALLOCATE(Acolnels(nels),Arownels(nels+1),anels(nels))!
 call coocsr(nels,nels,acoo,Arowcnels,Acolcoo,anels,Acolnels,Arownels)
 deallocate(Acolcoo,acoo)

 k=1
 ij=ndof*nels
 ALLOCATE(Acolcoo(ij),Arowcneq(ij),acoo(ij))
 elements_2: DO iel=1,nels
  g=g_g(:,iel)
  do i = 1,ndof !looking at the row vector related to each angle variable
    if (g(i) /= 0)then!  .and. g_c(g(i))/= 0
!--------equality constraint gradients wrt angle x----
      acoo(k) = k
      Acolcoo(k) = g(i)
      Arowcneq(k) = iel
      k=k+1
    endif
  enddo
 END DO elements_2
 ALLOCATE(Acolneq(knels),Arowneq(nels+1),aneq(knels))
 call coocsr(nels,knels,acoo,Arowcneq,Acolcoo,aneq,Acolneq,Arowneq)
!--------equality constraint hessian wrt displacement coo format-

 ij = 0

 nnzh = knels*2 + nels

 ALLOCATE(Hcol(nnzh),Hptr(nels+neq+1))
 Hptr(1) = 1
 Hcol = 0

 do i=1,nels !travelling from row to row
 !adding row pointers for equality constraints
  j = Arownels(i+1)-Arownels(i)
  Hcol(ij+1:ij+j) = Acolnels(Arownels(i):Arownels(i+1)-1)
  ij = ij+j
  k = Arowneq(i+1)-Arowneq(i)
  Hcol(ij+1:ij+k) = Acolneq(Arowneq(i):Arowneq(i+1)-1)+nels
  ij= ij+k
  Hptr(i+1) = Hptr(i)+j+k
 enddo
 allocate(Arow(nels))

 do i= 1,nels
  Arow(int(anels(i))) = i
 enddo

 k=1
 j=0
 do i = 1,nels
  do j = 1,nels
    if (Arowcnels(j) == i) then
      cntrnels(j) = int(Arow(j)+arowneq(i)-arowneq(1))
    endif
  enddo
 enddo

 deallocate(Arow)

 ! bottom part of the hessian matrix transpose of top right corner
 allocate(Arow(knels))

 do i= 1,knels
 Arow(int(aneq(i))) = i
 enddo

 k=1
 j=0
 do i = 1,nels
  do j = 1,knels
    if (Arowcneq(j) == i) then
      cntrneq(j) = int(Arow(j)+arownels(i+1)-arownels(1))
    endif
  enddo
 enddo

 k=1
 elements_3: DO iel=1,nels
  g=g_g(:,iel)
  do i = 1,ndof !looking at the row vector related to each angle variable
    if (g(i) /= 0)then!  .and. g_c(g(i))/= 0
!--------equality constraint gradients wrt angle x----
      acoo(k) = k
      Acolcoo(k) = iel
      Arowcneq(k) = g(i)
      k=k+1
    endif
  enddo
 END DO elements_3
 deallocate(Arowneq)
 ALLOCATE(Arowneq(neq+1))
 call coocsr(neq,knels,acoo,Arowcneq,Acolcoo,aneq,Acolneq,Arowneq)

 do i=1,neq !travelling from row to row
!---adding row pointers for equality constraints
  k = Arowneq(i+1)-Arowneq(i)
  Hcol(ij+1:ij+k) = Acolneq(Arowneq(i):Arowneq(i+1)-1)
  ij= ij+k
  Hptr(i+nels+1) = Hptr(i+nels)+k
 enddo

 do i= 1,knels
 Arow(int(aneq(i))) = i
 enddo

 k=1
 j=0
 do i = 1,neq
  do j = 1,knels
    if (Arowcneq(j) == i) then
      cntrneq2(j) = int(Arow(j))+nels+knels
    endif
  enddo
 enddo

 open(unit=104,FORM='unformatted',file='hessian.dat')
 write(104) cntrnels,cntrneq,cntrneq2,Hcol,Hptr

 DEALLOCATE(kdiag,g_g,g,Arowcnels,Acolcoo,cntrnels,cntrneq,acoo,cntrneq2  &
  ,Acolnels,Arownels,anels,Arowcneq,Acolneq,Arowneq,aneq,Hcol,Hptr,Arow)

return
END subroutine FemInitGG


subroutine FemInit(n,x)
!-------------------------------------------------------------------------
! Fri 02 Dec 2016 07:54:04 SAST Sa-aadat Parker
! Plane or axisymmetric strain analysis of an elastic solid
! using 3-, 6-, 10- or 15-node right-angled triangles or
! 4-, 8- or 9-node rectangular quadrilaterals. Mesh numbered
! in x(r)- or y(z)- direction.  
! Adapted from Program 5.1 in Programming the Finite Element Method-Smith 
!-------------------------------------------------------------------------
 
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15),rik=selected_int_kind(12)
 INTEGER::i,iel,loaded_nodes,ndim=2,ndof,nstep,npri,nlen,nodof=2        &
   ,k,nnz,hcount,hels,nxe,nye,layerv
 INTEGER,INTENT(OUT)::n
 INTEGER::nels,nn,nod,neq,kneq,knels,ierror
 REAL(iwp),INTENT(OUT)::x(*)
 REAL(iwp)::zero=0.0_iwp,dtim,pi
 CHARACTER(LEN=15)::argv,element,dir,type_2d,loading
 LOGICAL::solid=.TRUE.

!-----------------------dynamic arrays------------------------------------
 INTEGER,ALLOCATABLE::etype(:),g(:),g_g(:,:),g_num(:,:),nf(:,:),          &
   num(:),ncount(:)
 INTEGER(rik),ALLOCATABLE::kdiag(:)
 INTEGER(rik)::ij,nnzf,j
 REAL(iwp),ALLOCATABLE::coord(:,:),g_coord(:,:),km(:,:),kv(:),loads(:),   &
   x_coords(:),y_coords(:),xl(:)
 
!-----------------------input and initialisation--------------------------
      !fem setup element type and number
 argv = 'TopOpt'
 nlen = 6
!  OPEN(31,FILE=argv(1:nlen)//'.res')
 type_2d = 'plane'
 element = 'quadrilateral'
 loading = 'compression'!'tension'!
 nod = 4
 dir = 'y'
 nxe = 128
 nye = 128 !can only be an even number
 hels = 0
 layerv = 1 !number of layer variables in a symmetric layup
 CALL mesh_size_hole(element,nod,nels,hels,nn,nxe,nye,0)
 pi = dacos(0.d0)/2.d0

 allocate(xl(layerv))

 do i = 1,layerv*nels
  x(i) = 1.d-6+pi*(0.99d0*(-1)**(i+int((i-1)/(nye))))!pi*((-1)**floor(real((i-1)/nels)))!
 enddo
 ndof=nod*nodof
 ALLOCATE(nf(nodof,nn),g(ndof),g_coord(ndim,nn),km(ndof,ndof),            &
   coord(nod,ndim),g_num(nod,nels),y_coords(nye+1),etype(nels),num(nod),  &
   g_g(ndof,nels),x_coords(nxe+1),ncount(2+2*nxe))

 !  Setting mesh parameters
 x_coords = (/(1.0*i, i=0, nxe)/)
 y_coords = (/(-1.0*i, i=0, nye)/)  
 nf=1 ! setting the node freedom array to all free
 if (loading == 'tension') then
!  ! boundary conditions for plate in tension
  nf(1,1:nye*nod/4+1) =  0 ! setting symmetry boundary conditions at end
  nf(2,int(nye*nod/8+1)) =  0!+1
 elseif (loading == 'compression')then
 ! boundary conditions for plate in compression
   nf(:,nye+3) = 0!int((nn+1)/2)
 endif
 CALL formnf(nf,nodof,nn) 
 neq = MAXVAL(nf)
 n = layerv*nels!
!  ni = 0!1+2*((nye-1)*nxe+(nxe-1)*nye)
!  ne = 0!neq!
 ALLOCATE(loads(0:neq),kdiag(neq))
 kdiag=0
 !---------------------loop the elements to find global arrays sizes-----
 hcount=0
 ncount=0 !integer string which tracks nodes lost in hole.
 ncount(3)=5
 ncount(4)=nxe+4
 elements_1: DO iel=1,nels
   CALL geom_hole(element,iel,x_coords,nxe,y_coords,nye,coord,num,dir,    &
   nod,ndim,hels,hcount,ncount)
   CALL num_to_g(num,nf,g,ndof,nod,nodof,nn)
   g_num(:,iel)=num
   g_coord(:,num)=TRANSPOSE(coord)
   g_g(:,iel)=g
   CALL fkdiag(kdiag,g,ndof,neq)
 END DO elements_1
!  CALL mesh(g_coord,g_num,argv,6,32,nod,nels,ndim,nn)
 DO i=2,neq 
   kdiag(i)=kdiag(i)+kdiag(i-1) 
 END DO 
 ij=kdiag(neq)
 ALLOCATE(kv(ij),stat=ierror)

!------------------Loading applied in this section---
  loads=zero
 if (loading == 'tension') then
 ! loading for a plate in tension
  loaded_nodes = nye*nod/4+1 ! this will have to come from program
  do i = nn-nye*nod/4+1,nn-1
    loads(nf(:,i))=(/1.d0,0.d0/)
  enddo
  loads(nf(:,nn-nye*nod/4))=(/0.5d0,0.d0/)
  loads(nf(:,nn))=(/0.5d0,0.d0/)
 elseif (loading == 'compression')then
  ! loading for a plate in compression
  loaded_nodes = nn-(nye-1)*(nxe-1) ! this will have to come from program
    !left hand side nodes loaded
  do i = 2,nye
    loads(nf(:,i))=(/1.d0,0.d0/)/(nxe/16)
  enddo
  !right hand side nodes loaded
  do i = nn-nye+1,nn-1
    loads(nf(:,i))=(/-1.d0,0.d0/)/(nxe/16)
  enddo

  do i = 1, nxe-1
  !top nodes loaded
    loads(nf(:,ncount(i+4)))=(/0.d0,-1.d0/)/(nxe/16)
    !bottom nodes loaded
    loads(nf(:,ncount(i+nxe+3)))=(/0.d0,1.d0/)/(nxe/16)
  enddo

  ! bottom left corner
  loads(nf(:,nye+1))=(/0.5d0,0.5d0/)/(nxe/16)
  ! top right corner
  loads(nf(:,nn-nye))=(/-0.5d0,-0.5d0/)/(nxe/16)
  ! bottom right corner
  loads(nf(:,nn))=(/-0.5d0,0.5d0/)/(nxe/16)
  ! top left corner
  loads(nf(:,1))=(/0.5d0,-0.5d0/)/(nxe/16)
 endif
!-------------------element stiffness integration and assembly--------
 kv=zero 
!  nnzf=ij+ne+ni
 k=1

 elements_3: DO iel=1,nels
   do i = 1,layerv
    xl(i) = x(iel+nels*(i-1))
   enddo
   num=g_num(:,iel)
   coord=TRANSPOSE(g_coord(:,num)) 
   g=g_g(:,iel) 
   km=zero
  do i = 1,ndof
    if (g(i) /= 0 )then !.and. g_c(g(i))/= 0
!--------equality constraint gradients wrt thickness x---- 
      k=k+1
    endif
  enddo
  SELECT CASE(nod)
    CASE(4)
      CALL Q4_ortho_n(km,xl,layerv)
    CASE DEFAULT
      WRITE(*,'(A)') "# Element type not recognised"
   END SELECT
  CALL fsparv(kv,km,g,kdiag,ndof,neq,nnzf)
 END DO elements_3
!-------plot mesh with loadings---------------------------
 open(unit=101,FORM='unformatted',file='store.dat')
 write(101) g_g,kdiag,loads,g_num,g_coord
 nstep = 5
 npri = 1
 dtim = 1
 CALL mesh_ensi(argv,nlen,g_coord,g_num,element,x,nf,loads,               &
                nstep,npri,dtim,solid,nod,nels,nodof,nn,ndim,neq)

!---------------aptr, acol and nnz calculation---------------------
 knels = 0 !k-1 !equality constraint nnz
!--------equality constraint gradients wrt displacement coo format-
!--------this is equal to the global stiffness matrix--------------
!  k=2 ! start from second row
!  DO i = 2,neq
!     DO j = kdiag(i-1)+1,kdiag(i)
!       if (kv(j) /= 0.d0) then
!         k=k+1
!         if (i /= i-(kdiag(i)-j)) then
!           k=k+1
!         endif
!       endif  
!     ENDDO
!  ENDDO

 kneq= 0!k-1
 nnz = kneq+knels
! 
! !  normal displacement start points SAND only
!  CALL sparin(kv,kdiag,neq,ij)
!  CALL spabac(kv,loads,kdiag,neq,ij)
! 
!  do i = layerv*nels+1,n
!     x(i) =loads(i-layerv*nels)*(1.d0-1.d-4)!1.d-3!!!!1.d-3*(-1)**(i-1)*modulo(i,5)!
!  enddo

 DEALLOCATE(nf,g,g_coord,km,coord,g_num,y_coords,etype,num,g_g,x_coords,  &
 ncount,loads,kdiag,kv)
 
 open(unit=102,FORM='unformatted',file='integers.dat')
 write(102) nels,nod,neq,layerv,nn
return
END subroutine FemInit

subroutine Fem_g(x,n,gf)
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
 INTEGER::i,iel,k,loaded_nodes,ndim=2,ndofksl,nodof=2        &
   ,ndof,l,nels,nn,nod,neq,layerv
 INTEGER,INTENT(IN)::n
 REAL(iwp),INTENT(IN)::x(n)
 REAL(iwp),INTENT(OUT)::gf(n)
 REAL(iwp)::zero=0.0_iwp,dx
!-----------------------dynamic arrays------------------------------------
 INTEGER,ALLOCATABLE::g(:),g_g(:,:)
 INTEGER(rik),ALLOCATABLE::kdiag(:)
 INTEGER(rik)::ij,j
 REAL(iwp),ALLOCATABLE::km(:,:),kv(:),loads(:),disp(:),g_c(:),xl(:),xc(:)
!-----------------------input and initialisation--------------------------
 rewind(102)
 read(102)nels,nod,neq,layerv,nn

 ndof=nod*nodof
 ALLOCATE(g(ndof),km(ndof,ndof),g_g(ndof,nels))

 ALLOCATE(loads(0:neq),kdiag(neq),disp(0:neq),g_c(0:neq))
!-------------------element stiffness integration and assembly--------
 ij=kdiag(neq)
 disp(0)=zero
 rewind(101)
 read(101) g_g,kdiag,loads
 ALLOCATE(kv(kdiag(neq)),xl(layerv),xc(n))
 kv=zero ! global stiffness matrix
 
!  rewind(103)
!  read(103)xc,disp
!  
!  dx =0.d0
!  do i = 1,10
!   dx = abs(xc(i) - x(i)) + dx
!  enddo
!  
!  if (dx .gt. 1.d-5) then
  !-------------------NAND variable displacements---------------------
  disp=loads
  !-------------------element stiffness integration and assembly------
  kv=zero 
  ij=kdiag(neq) 
  elements_2: DO iel=1,nels
    do i = 1,layerv
      xl(i) = x(iel+nels*(i-1))
    enddo
    g=g_g(:,iel)
    km=zero
      SELECT CASE(nod)
      CASE(4)
          CALL Q4_ortho_n(km,xl,layerv)
      CASE DEFAULT
          WRITE(*,'(A)') "# Element type not recognised"
    END SELECT
    CALL fsparv(kv,km,g,kdiag,ndof,neq,ij)
  END DO elements_2

  ! displacement calculation from u = k^(-1)F
  CALL sparin(kv,kdiag,neq,ij)
  CALL spabac(kv,disp,kdiag,neq,ij)
!  endif
!--------Gradients wrt element angles----
elements_1: DO iel=1,nels
  do i = 1,layerv
    xl(i) = x(iel+nels*(i-1))
  enddo
  g=g_g(:,iel) 
  km=zero
     SELECT CASE(nod)
    CASE(4)
        CALL Q4_ortho_nd(km,x(iel),layerv)
    CASE DEFAULT
        WRITE(*,'(A)') "# Element type not recognised"
   END SELECT
  g_c(g) = MATMUL(km,disp(g))
  gf(iel) = -dot_product(g_c(g),disp(g))
 END DO elements_1

 DEALLOCATE(g,km,g_g,loads,kdiag,disp,g_c,kv)
return
END subroutine Fem_g

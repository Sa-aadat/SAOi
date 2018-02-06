subroutine Fem(x,f,n)

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
 INTEGER::i,iel,ndof,nodof=2,nels,nod,neq,layerv
 INTEGER,INTENT(IN)::n
 REAL(iwp),INTENT(IN)::x(n)
 REAL(iwp),INTENT(OUT)::f
 REAL(iwp)::zero=0.0_iwp

!-----------------------dynamic arrays------------------------------------
 INTEGER,ALLOCATABLE::g(:),g_g(:,:)
 INTEGER(rik)::ij
 INTEGER(rik),ALLOCATABLE::kdiag(:)
 REAL(iwp),ALLOCATABLE::disp(:),kv(:),loads(:),km(:,:),xl(:)

!-----------------------input and initialisation--------------------------
 rewind(102)
 read(102)nels,nod,neq,layerv

 ndof=nod*nodof
 ALLOCATE(g(ndof),km(ndof,ndof),g_g(ndof,nels),loads(0:neq)      &
    ,kdiag(neq),disp(0:neq),xl(layerv))

 rewind(101)
 read(101) g_g,kdiag,loads
 ALLOCATE(kv(kdiag(neq)))
 !-------------------NAND variable displacements---------------------------
 disp=zero!(0)
 disp=loads
 !-----------------------element stiffness integration and assembly--------
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
 
!   open(unit=103,FORM='unformatted',file='u.dat')
!  write(103) x, disp

!----------------------calculate objective--------------------------------
 f = dot_product(loads,disp)

 DEALLOCATE(g,km,g_g,loads,kdiag,disp,kv,xl)
return
END subroutine Fem

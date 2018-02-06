SUBROUTINE formnf(nf,nodof,nn)
!
! This subroutine forms the nf matrix.
!
 IMPLICIT NONE
 INTEGER::i,j,m,nodof, nn
 INTEGER,INTENT(IN OUT)::nf(nodof,nn)
 m=0
 DO j=1,UBOUND(nf,2)!nn
   DO i=1,UBOUND(nf,1)!nodof
     IF(nf(i,j)/=0)THEN
       m=m+1
       nf(i,j)=m
     END IF
   END DO
 END DO
RETURN
END SUBROUTINE formnf

SUBROUTINE geom_rect(element,iel,x_coords,nxe,y_coords,nye,               &
coord,num,dir,nod,ndim)

! This subroutine forms the coordinates and connectivity for a
! rectangular mesh of rt. angled triangular elements (3, 6, 10 or 15-node)
! or quadrilateral elements (4, 8 or 9-node) counting in the
! x- or y-dir. 
!
 IMPLICIT NONE
 INTEGER::ip,iq,jel,fac1,nod,nxe,nye,ndim
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::x_coords(nxe+1),y_coords(nye+1)
 REAL(iwp),INTENT(OUT)::coord(nod,ndim)
 CHARACTER(LEN=15),INTENT(IN)::element
 CHARACTER(LEN=1),INTENT(IN)::dir
 INTEGER,INTENT(IN)::iel
 INTEGER,INTENT(OUT)::num(nod)
 REAL(iwp)::pt5=0.5_iwp,two=2.0_iwp,d3=3.0_iwp 
 nxe=UBOUND(x_coords,1)-1
 nod=UBOUND(num,1)
 IF(element=='triangle')THEN
   nye=(UBOUND(y_coords,1)-1)*2
   IF(dir=='x'.OR.dir=='r')THEN
     jel=2*nxe*((iel-1)/(2*nxe))
     ip=(iel-jel+1)/2
     iq=2*((iel-1)/(2*nxe)+1)-1+((iel/2)*2)/iel
   ELSE  
     jel=(iel-1)/nye
     ip=jel+1
     iq=iel-nye*jel
   END IF
   SELECT CASE(nod)
   CASE(3)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(nxe+1)*(iq-1)/2+ip
         num(2)=num(1)+1              
         num(3)=(nxe+1)*(iq+1)/2+ip
       ELSE
         num(1)=(ip-1)*(nye+2)/2+(iq+1)/2
         num(2)=num(1)+(nye+2)/2
         num(3)=num(1)+1
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(1,2)=y_coords((iq+1)/2)
       coord(2,1)=x_coords(ip+1)   
       coord(2,2)=y_coords((iq+1)/2)
       coord(3,1)=x_coords(ip)   
       coord(3,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(nxe+1)*iq/2+ip+1     
         num(2)=num(1)-1               
         num(3)=(nxe+1)*(iq-2)/2+ip+1
       ELSE
         num(1)=ip*(nye+2)/2+(iq+2)/2
         num(2)=(ip-1)*(nye+2)/2+(iq+1)/2+1
         num(3)=num(1)-1
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(1,2)=y_coords((iq+2)/2)
       coord(2,1)=x_coords(ip)   
       coord(2,2)=y_coords((iq+2)/2)
       coord(3,1)=x_coords(ip+1) 
       coord(3,2)=y_coords(iq/2)
     END IF
   CASE(6)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(iq-1)*(2*nxe+1)+2*ip-1
         num(2)=num(1)+1 
         num(3)=num(1)+2 
         num(4)=(iq-1)*(2*nxe+1)+2*nxe+2*ip+1
         num(5)=(iq+1)*(2*nxe+1)+2*ip-1
         num(6)=num(4)-1 
       ELSE
         num(1)=2*(nye+1)*(ip-1)+iq
         num(2)=2*(nye+1)*(ip-1)+nye+1+iq
         num(3)=2*(nye+1)*ip+iq
         num(4)=num(2)+1
         num(5)=num(1)+2 
         num(6)=num(1)+1
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(1,2)=y_coords((iq+1)/2)
       coord(3,1)=x_coords(ip+1)   
       coord(3,2)=y_coords((iq+1)/2)
       coord(5,1)=x_coords(ip)   
       coord(5,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=iq*(2*nxe+1)+2*ip+1
         num(2)=num(1)-1 
         num(3)=num(1)-2 
         num(4)=(iq-2)*(2*nxe+1)+2*nxe+2*ip+1
         num(5)=(iq-2)*(2*nxe+1)+2*ip+1
         num(6)=num(4)+1 
       ELSE 
         num(1)=2*(nye+1)*ip+iq+1 
         num(2)=2*(nye+1)*(ip-1)+nye+iq+2
         num(3)=2*(nye+1)*(ip-1)+iq+1
         num(4)=num(2)-1 
         num(5)=num(1)-2
         num(6)=num(1)-1
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(1,2)=y_coords((iq+2)/2)
       coord(3,1)=x_coords(ip)   
       coord(3,2)=y_coords((iq+2)/2)
       coord(5,1)=x_coords(ip+1) 
       coord(5,2)=y_coords(iq/2)
     END IF
     coord(2,:)=pt5*(coord(1,:)+coord(3,:))
     coord(4,:)=pt5*(coord(3,:)+coord(5,:))
     coord(6,:)=pt5*(coord(5,:)+coord(1,:))
   CASE(10)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(iq-1)/2*(3*nxe+1)*3+3*ip-2
         num(2)=num(1)+1
         num(3)=num(1)+2
         num(4)=num(1)+3
         num(5)=(iq-1)/2*(3*nxe+1)*3+3*nxe+1+3*ip
         num(6)=(iq-1)/2*(3*nxe+1)*3+6*nxe+2+3*ip-1
         num(7)=(iq-1)/2*(3*nxe+1)*3+9*nxe+3+3*ip-2
         num(8)=num(6)-1
         num(9)=num(5)-2
         num(10)=num(9)+1
       ELSE
         num(1)=(9*(nye-2)/2+12)*(ip-1)+3*(iq-1)/2+1
         num(2)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)/2+4+3*(iq-1)/2+1
         num(3)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)+8+3*(iq-1)/2+1
         num(4)=(9*(nye-2)/2+12)*(ip-1)+9*(nye-2)/2+12+3*(iq-1)/2+1
         num(5)=num(3)+1 
         num(6)=num(2)+2
         num(7)=num(1)+3
         num(8)=num(1)+2
         num(9)=num(1)+1
         num(10)=num(2)+1
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(2,1)=x_coords(ip)+(x_coords(ip+1)-x_coords(ip))/d3
       coord(3,1)=x_coords(ip)+two*(x_coords(ip+1)-x_coords(ip))/d3
       coord(4,1)=x_coords(ip+1)
       coord(4,2)=y_coords((iq+1)/2)
       coord(5,2)=y_coords((iq+1)/2)+                                     &
         (y_coords((iq+3)/2)-y_coords((iq+1)/2))/d3
       coord(6,2)=y_coords((iq+1)/2)+                                     &
         two*(y_coords((iq+3)/2)-y_coords((iq+1)/2))/d3
       coord(7,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(iq-2)/2*(3*nxe+1)*3+9*nxe+3+3*ip+1
         num(2)=num(1)-1
         num(3)=num(1)-2
         num(4)=num(1)-3
         num(5)=(iq-2)/2*(3*nxe+1)*3+6*nxe+2+3*ip-1
         num(6)=(iq-2)/2*(3*nxe+1)*3+3*nxe+1+3*ip
         num(7)=(iq-2)/2*(3*nxe+1)*3+3*ip+1
         num(8)=num(6)+1
         num(9)=num(5)+2
         num(10)=num(9)-1
       ELSE
         num(1)=(9*(nye-2)/2+12)*(ip-1)+9*(nye-2)/2+12+3*iq/2+1
         num(2)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)+8+3*iq/2+1
         num(3)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)/2+4+3*iq/2+1
         num(4)=(9*(nye-2)/2+12)*(ip-1)+3*iq/2+1
         num(5)=num(3)-1
         num(6)=num(2)-2
         num(7)=num(1)-3
         num(8)=num(1)-2
         num(9)=num(1)-1
         num(10)=num(2)-1
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(2,1)=x_coords(ip+1)-(x_coords(ip+1)-x_coords(ip))/d3
       coord(3,1)=x_coords(ip+1)-two*(x_coords(ip+1)-x_coords(ip))/d3
       coord(4,1)=x_coords(ip)
       coord(4,2)=y_coords((iq+2)/2)
       coord(5,2)=y_coords((iq+2)/2)-(y_coords((iq+2)/2)-y_coords(iq/2))/d3
       coord(6,2)=y_coords((iq+2)/2)-                                     &
         two*(y_coords((iq+2)/2)-y_coords(iq/2))/d3
       coord(7,2) =y_coords(iq/2)
     END IF
     coord(5,1)=coord(3,1)
     coord(6,1)=coord(2,1)
     coord(7,1)=coord(1,1)
     coord(8,1)=coord(1,1)
     coord(9,1)=coord(1,1)
     coord(10,1)=coord(2,1)
     coord(1,2)=coord(4,2)
     coord(2,2)=coord(4,2)
     coord(3,2)=coord(4,2)
     coord(8,2)=coord(6,2)
     coord(9,2)=coord(5,2)
     coord(10,2)=coord(5,2)
   CASE(15)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
       fac1=4*(4*nxe+1)*(iq-1)/2
         num(1)=fac1+4*ip-3
         num(2)=num(1)+1
         num(3)=num(1)+2
         num(4)=num(1)+3
         num(5)=num(1)+4
         num(6)=fac1+ 4*nxe+1+4*ip
         num(7)=fac1+ 8*nxe+1+4*ip
         num(8)=fac1+12*nxe+1+4*ip
         num(9)=fac1+16*nxe+1+4*ip
         num(10)=num(8)-1
         num(11)=num(7)-2
         num(12)=num(6)-3
         num(13)=num(12)+1
         num(14)=num(12)+2
         num(15)=num(11)+1
       ELSE
         fac1=4*(2*nye+1)*(ip-1)+2*iq-1 
         num(1)=fac1
         num(2)=fac1+2*nye+1
         num(3)=fac1+4*nye+2 
         num(4)=fac1+6*nye+3 
         num(5)=fac1+8*nye+4
         num(6)=fac1+6*nye+4 
         num(7)=fac1+4*nye+4 
         num(8)=fac1+2*nye+4
         num(9)=fac1+4 
         num(10)=fac1+3 
         num(11)=fac1+2 
         num(12)=fac1+1
         num(13)=fac1+2*nye+2 
         num(14)=fac1+4*nye+3
         num(15)=fac1+2*nye+3  
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(1,2)=y_coords((iq+1)/2)
       coord(5,1)=x_coords(ip+1)   
       coord(5,2)=y_coords((iq+1)/2)
       coord(9,1)=x_coords(ip)   
       coord(9,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         fac1=4*(4*nxe+1)*(iq-2)/2
         num(1)=fac1+16*nxe+5+4*ip
         num(2)=num(1)-1
         num(3)=num(1)-2
         num(4)=num(1)-3
         num(5)=num(1)-4
         num(6)=fac1+12*nxe+1+4*ip
         num(7)=fac1+8*nxe+1+4*ip
         num(8)=fac1+4*nxe+1+4*ip
         num(9)=fac1+4*ip+1
         num(10)=num(8)+1
         num(11)=num(7)+2
         num(12)=num(6)+3
         num(13)=num(12)-1
         num(14)=num(12)-2
         num(15)=num(11)-1
       ELSE
         fac1=4*(2*nye+1)*(ip-1)+2*iq+8*nye+5 
         num(1)=fac1 
         num(2)=fac1-2*nye-1
         num(3)=fac1-4*nye-2 
         num(4)=fac1-6*nye-3 
         num(5)=fac1-8*nye-4
         num(6)=fac1-6*nye-4  
         num(7)=fac1-4*nye-4 
         num(8)=fac1-2*nye-4
         num(9)=fac1-4
         num(10)=fac1-3 
         num(11)=fac1-2 
         num(12)=fac1-1
         num(13)=fac1-2*nye-2  
         num(14)=fac1-4*nye-3
         num(15)=fac1-2*nye-3 
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(1,2)=y_coords((iq+2)/2)
       coord(5,1)=x_coords(ip)   
       coord(5,2)=y_coords((iq+2)/2)
       coord(9,1)=x_coords(ip+1) 
       coord(9,2)=y_coords(iq/2)
     END IF
     coord(3,:)=pt5*(coord(1,:)+coord(5,:))
     coord(7,:)=pt5*(coord(5,:)+coord(9,:))
     coord(11,:)=pt5*(coord(9,:)+coord(1,:))
     coord(2,:)=pt5*(coord(1,:)+coord(3,:))
     coord(4,:)=pt5*(coord(3,:)+coord(5,:))
     coord(6,:)=pt5*(coord(5,:)+coord(7,:))
     coord(8,:)=pt5*(coord(7,:)+coord(9,:))
     coord(10,:)=pt5*(coord(9,:)+coord(11,:))
     coord(12,:)=pt5*(coord(11,:)+coord(1,:))
     coord(15,:)=pt5*(coord(7,:)+coord(11,:))
     coord(14,:)=pt5*(coord(3,:)+coord(7,:))
     coord(13,:)=pt5*(coord(2,:)+coord(15,:))
   CASE DEFAULT
     WRITE(31,'(a)')"Wrong number of nodes for triangular element"
     STOP
   END SELECT
 ELSE
   nye=UBOUND(y_coords,1)-1
   IF(dir=='x'.OR.dir=='r')THEN
     iq=(iel-1)/nxe+1
     ip=iel-(iq-1)*nxe
   ELSE
     ip=(iel-1)/nye+1!iel manipulated
     iq=iel-(ip-1)*nye
   END IF
   SELECT CASE(nod)
   CASE(4)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(nxe+1)+ip
       num(2)=(iq-1)*(nxe+1)+ip
       num(3)=num(2)+1
       num(4)=num(1)+1
     ELSE
       num(1)=(ip-1)*(nye+1)+iq+1
       num(2)=num(1)-1
       num(3)=ip*(nye+1)+iq
       num(4)=num(3)+1
     END IF
!
     coord(1:2,1)=x_coords(ip)
     coord(3:4,1)=x_coords(ip+1)
     coord(1,2)=y_coords(iq+1)
     coord(2:3,2)=y_coords(iq)
     coord(4,2)=coord(1,2)
   CASE(5)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(2*nxe+1)+ip
       num(2)=(iq-1)*(2*nxe+1)+ip
       num(3)=num(2)+1
       num(4)=num(1)+1
       num(5)=iq*(2*nxe+1)+ip-nxe
     ELSE
       num(1)=(ip-1)*(2*nye+1)+iq+1
       num(2)=num(1)-1
       num(3)=ip*(2*nye+1)+iq
       num(4)=num(3)+1
       num(5)=ip*(2*nye+1)+iq-nye
     END IF
!
     coord(1:2,1)=x_coords(ip)
     coord(3:4,1)=x_coords(ip+1)
     coord(1,2)=y_coords(iq+1)
     coord(2:3,2)=y_coords(iq)
     coord(4,2)=coord(1,2)
     coord(5,:)=0.25_iwp*(coord(1,:)+coord(2,:)+coord(3,:)+coord(4,:))
   CASE(8)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(3*nxe+2)+2*ip-1                 
       num(2)=iq*(3*nxe+2)+ip-nxe-1
       num(3)=(iq-1)*(3*nxe+2)+2*ip-1   
       num(4)=num(3)+1
       num(5)=num(4)+1
       num(6)=num(2)+1
       num(7)=num(1)+2
       num(8)=num(1)+1
     ELSE
       num(1)=(ip-1)*(3*nye+2)+2*iq+1
       num(2)=num(1)-1
       num(3)=num(1)-2
       num(4)=(ip-1)*(3*nye+2)+2*nye+iq+1
       num(5)=ip*(3*nye+2)+2*iq-1
       num(6)=num(5)+1
       num(7)=num(5)+2
       num(8)=num(4)+1
     END IF
!
     coord(1:3,1)=x_coords(ip)
     coord(5:7,1)=x_coords(ip+1)
     coord(4,1)=pt5*(coord(3,1)+coord(5,1))
     coord(8,1)=pt5*(coord(7,1)+coord(1,1))
     coord(1,2)=y_coords(iq+1)
     coord(7:8,2)=y_coords(iq+1)
     coord(3:5,2)=y_coords(iq)
     coord(2,2)=pt5*(coord(1,2)+coord(3,2))
     coord(6,2)=pt5*(coord(5,2)+coord(7,2))
   CASE(9)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(4*nxe+2)+2*ip-1
       num(2)=iq*(4*nxe+2)+2*ip-nxe-4
       num(3)= (iq-1)*(4*nxe+2)+2*ip-1
       num(4)=num(3)+1
       num(5)=num(4)+1
       num(6)=num(2)+2
       num(7)=num(1)+2
       num(8)=num(1)+1
       num(9)=num(2)+1
     ELSE
       num(1)=(ip-1)*2*(2*nye+1)+2*iq+1
       num(2)=num(1)-1
       num(3)=num(1)-2
       num(4)=(ip-1)*2*(2*nye+1)+2*nye+2*iq
       num(5)=ip*2*(2*nye+1)+2*iq-1
       num(6)=num(5)+1
       num(7)=num(5)+2
       num(8)=num(4)+2
       num(9)=num(4)+1
     END IF
!
     coord(1:3,1)=x_coords(ip)
     coord(5:7,1)=x_coords(ip+1)
     coord(4,1)=pt5*(coord(3,1)+coord(5,1))
     coord(8,1)=pt5*(coord(7,1)+coord(1,1))
     coord(1,2)=y_coords(iq+1)
     coord(7:8,2)=y_coords(iq+1)
     coord(3:5,2)=y_coords(iq)
     coord(2,2)=pt5*(coord(1,2)+coord(3,2))
     coord(6,2)=pt5*(coord(5,2)+coord(7,2))
     coord(9,:)=pt5*(coord(4,:)+coord(8,:))
   CASE DEFAULT
     WRITE(31,'(a)')"Wrong number of nodes for quadrilateral element"
     STOP
   END SELECT
 END IF
RETURN
END SUBROUTINE geom_rect

SUBROUTINE geom_hole(element,iel,x_coords,nxe,y_coords,nye,               &
coord,num,dir,nod,ndim,hels,hcount,ncount)

! This subroutine forms the coordinates and connectivity for a
! rectangular mesh of rt. angled triangular elements (3, 6, 10 or 15-node)
! or quadrilateral elements (4, 8 or 9-node) counting in the
! x- or y-dir.
!
 IMPLICIT NONE
 INTEGER::ip,iq,jel,fac1,nod,nxe,nye,ndim,i
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::x_coords(nxe+1),y_coords(nye+1)
 REAL(iwp),INTENT(OUT)::coord(nod,ndim)
 CHARACTER(LEN=15),INTENT(IN)::element
 CHARACTER(LEN=1),INTENT(IN)::dir
 INTEGER,INTENT(INOUT)::hcount,ncount(2+2*nxe)
 INTEGER,INTENT(IN)::hels,iel
 INTEGER,INTENT(OUT)::num(nod)
 REAL(iwp)::rectmidpt(1,2),elemidpoint(1,2),dels,dnods
 REAL(iwp)::pt5=0.5_iwp,two=2.0_iwp,d3=3.0_iwp
 nxe=UBOUND(x_coords,1)-1
 nod=UBOUND(num,1)
 IF(element=='triangle')THEN
   nye=(UBOUND(y_coords,1)-1)*2
   IF(dir=='x'.OR.dir=='r')THEN
     jel=2*nxe*((iel-1)/(2*nxe))
     ip=(iel-jel+1)/2
     iq=2*((iel-1)/(2*nxe)+1)-1+((iel/2)*2)/iel
   ELSE
     jel=(iel-1)/nye
     ip=jel+1
     iq=iel-nye*jel
   END IF
   SELECT CASE(nod)
   CASE(3)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(nxe+1)*(iq-1)/2+ip
         num(2)=num(1)+1
         num(3)=(nxe+1)*(iq+1)/2+ip
       ELSE
         num(1)=(ip-1)*(nye+2)/2+(iq+1)/2
         num(2)=num(1)+(nye+2)/2
         num(3)=num(1)+1
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(1,2)=y_coords((iq+1)/2)
       coord(2,1)=x_coords(ip+1)
       coord(2,2)=y_coords((iq+1)/2)
       coord(3,1)=x_coords(ip)
       coord(3,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(nxe+1)*iq/2+ip+1
         num(2)=num(1)-1
         num(3)=(nxe+1)*(iq-2)/2+ip+1
       ELSE
         num(1)=ip*(nye+2)/2+(iq+2)/2
         num(2)=(ip-1)*(nye+2)/2+(iq+1)/2+1
         num(3)=num(1)-1
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(1,2)=y_coords((iq+2)/2)
       coord(2,1)=x_coords(ip)
       coord(2,2)=y_coords((iq+2)/2)
       coord(3,1)=x_coords(ip+1)
       coord(3,2)=y_coords(iq/2)
     END IF
   CASE(6)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(iq-1)*(2*nxe+1)+2*ip-1
         num(2)=num(1)+1
         num(3)=num(1)+2
         num(4)=(iq-1)*(2*nxe+1)+2*nxe+2*ip+1
         num(5)=(iq+1)*(2*nxe+1)+2*ip-1
         num(6)=num(4)-1
       ELSE
         num(1)=2*(nye+1)*(ip-1)+iq
         num(2)=2*(nye+1)*(ip-1)+nye+1+iq
         num(3)=2*(nye+1)*ip+iq
         num(4)=num(2)+1
         num(5)=num(1)+2
         num(6)=num(1)+1
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(1,2)=y_coords((iq+1)/2)
       coord(3,1)=x_coords(ip+1)
       coord(3,2)=y_coords((iq+1)/2)
       coord(5,1)=x_coords(ip)
       coord(5,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=iq*(2*nxe+1)+2*ip+1
         num(2)=num(1)-1
         num(3)=num(1)-2
         num(4)=(iq-2)*(2*nxe+1)+2*nxe+2*ip+1
         num(5)=(iq-2)*(2*nxe+1)+2*ip+1
         num(6)=num(4)+1
       ELSE
         num(1)=2*(nye+1)*ip+iq+1
         num(2)=2*(nye+1)*(ip-1)+nye+iq+2
         num(3)=2*(nye+1)*(ip-1)+iq+1
         num(4)=num(2)-1
         num(5)=num(1)-2
         num(6)=num(1)-1
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(1,2)=y_coords((iq+2)/2)
       coord(3,1)=x_coords(ip)
       coord(3,2)=y_coords((iq+2)/2)
       coord(5,1)=x_coords(ip+1)
       coord(5,2)=y_coords(iq/2)
     END IF
     coord(2,:)=pt5*(coord(1,:)+coord(3,:))
     coord(4,:)=pt5*(coord(3,:)+coord(5,:))
     coord(6,:)=pt5*(coord(5,:)+coord(1,:))
   CASE(10)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(iq-1)/2*(3*nxe+1)*3+3*ip-2
         num(2)=num(1)+1
         num(3)=num(1)+2
         num(4)=num(1)+3
         num(5)=(iq-1)/2*(3*nxe+1)*3+3*nxe+1+3*ip
         num(6)=(iq-1)/2*(3*nxe+1)*3+6*nxe+2+3*ip-1
         num(7)=(iq-1)/2*(3*nxe+1)*3+9*nxe+3+3*ip-2
         num(8)=num(6)-1
         num(9)=num(5)-2
         num(10)=num(9)+1
       ELSE
         num(1)=(9*(nye-2)/2+12)*(ip-1)+3*(iq-1)/2+1
         num(2)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)/2+4+3*(iq-1)/2+1
         num(3)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)+8+3*(iq-1)/2+1
         num(4)=(9*(nye-2)/2+12)*(ip-1)+9*(nye-2)/2+12+3*(iq-1)/2+1
         num(5)=num(3)+1
         num(6)=num(2)+2
         num(7)=num(1)+3
         num(8)=num(1)+2
         num(9)=num(1)+1
         num(10)=num(2)+1
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(2,1)=x_coords(ip)+(x_coords(ip+1)-x_coords(ip))/d3
       coord(3,1)=x_coords(ip)+two*(x_coords(ip+1)-x_coords(ip))/d3
       coord(4,1)=x_coords(ip+1)
       coord(4,2)=y_coords((iq+1)/2)
       coord(5,2)=y_coords((iq+1)/2)+                                     &
         (y_coords((iq+3)/2)-y_coords((iq+1)/2))/d3
       coord(6,2)=y_coords((iq+1)/2)+                                     &
         two*(y_coords((iq+3)/2)-y_coords((iq+1)/2))/d3
       coord(7,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(iq-2)/2*(3*nxe+1)*3+9*nxe+3+3*ip+1
         num(2)=num(1)-1
         num(3)=num(1)-2
         num(4)=num(1)-3
         num(5)=(iq-2)/2*(3*nxe+1)*3+6*nxe+2+3*ip-1
         num(6)=(iq-2)/2*(3*nxe+1)*3+3*nxe+1+3*ip
         num(7)=(iq-2)/2*(3*nxe+1)*3+3*ip+1
         num(8)=num(6)+1
         num(9)=num(5)+2
         num(10)=num(9)-1
       ELSE
         num(1)=(9*(nye-2)/2+12)*(ip-1)+9*(nye-2)/2+12+3*iq/2+1
         num(2)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)+8+3*iq/2+1
         num(3)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)/2+4+3*iq/2+1
         num(4)=(9*(nye-2)/2+12)*(ip-1)+3*iq/2+1
         num(5)=num(3)-1
         num(6)=num(2)-2
         num(7)=num(1)-3
         num(8)=num(1)-2
         num(9)=num(1)-1
         num(10)=num(2)-1
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(2,1)=x_coords(ip+1)-(x_coords(ip+1)-x_coords(ip))/d3
       coord(3,1)=x_coords(ip+1)-two*(x_coords(ip+1)-x_coords(ip))/d3
       coord(4,1)=x_coords(ip)
       coord(4,2)=y_coords((iq+2)/2)
       coord(5,2)=y_coords((iq+2)/2)-(y_coords((iq+2)/2)-y_coords(iq/2))/d3
       coord(6,2)=y_coords((iq+2)/2)-                                     &
         two*(y_coords((iq+2)/2)-y_coords(iq/2))/d3
       coord(7,2) =y_coords(iq/2)
     END IF
     coord(5,1)=coord(3,1)
     coord(6,1)=coord(2,1)
     coord(7,1)=coord(1,1)
     coord(8,1)=coord(1,1)
     coord(9,1)=coord(1,1)
     coord(10,1)=coord(2,1)
     coord(1,2)=coord(4,2)
     coord(2,2)=coord(4,2)
     coord(3,2)=coord(4,2)
     coord(8,2)=coord(6,2)
     coord(9,2)=coord(5,2)
     coord(10,2)=coord(5,2)
   CASE(15)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
       fac1=4*(4*nxe+1)*(iq-1)/2
         num(1)=fac1+4*ip-3
         num(2)=num(1)+1
         num(3)=num(1)+2
         num(4)=num(1)+3
         num(5)=num(1)+4
         num(6)=fac1+ 4*nxe+1+4*ip
         num(7)=fac1+ 8*nxe+1+4*ip
         num(8)=fac1+12*nxe+1+4*ip
         num(9)=fac1+16*nxe+1+4*ip
         num(10)=num(8)-1
         num(11)=num(7)-2
         num(12)=num(6)-3
         num(13)=num(12)+1
         num(14)=num(12)+2
         num(15)=num(11)+1
       ELSE
         fac1=4*(2*nye+1)*(ip-1)+2*iq-1
         num(1)=fac1
         num(2)=fac1+2*nye+1
         num(3)=fac1+4*nye+2
         num(4)=fac1+6*nye+3
         num(5)=fac1+8*nye+4
         num(6)=fac1+6*nye+4
         num(7)=fac1+4*nye+4
         num(8)=fac1+2*nye+4
         num(9)=fac1+4
         num(10)=fac1+3
         num(11)=fac1+2
         num(12)=fac1+1
         num(13)=fac1+2*nye+2
         num(14)=fac1+4*nye+3
         num(15)=fac1+2*nye+3
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(1,2)=y_coords((iq+1)/2)
       coord(5,1)=x_coords(ip+1)
       coord(5,2)=y_coords((iq+1)/2)
       coord(9,1)=x_coords(ip)
       coord(9,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         fac1=4*(4*nxe+1)*(iq-2)/2
         num(1)=fac1+16*nxe+5+4*ip
         num(2)=num(1)-1
         num(3)=num(1)-2
         num(4)=num(1)-3
         num(5)=num(1)-4
         num(6)=fac1+12*nxe+1+4*ip
         num(7)=fac1+8*nxe+1+4*ip
         num(8)=fac1+4*nxe+1+4*ip
         num(9)=fac1+4*ip+1
         num(10)=num(8)+1
         num(11)=num(7)+2
         num(12)=num(6)+3
         num(13)=num(12)-1
         num(14)=num(12)-2
         num(15)=num(11)-1
       ELSE
         fac1=4*(2*nye+1)*(ip-1)+2*iq+8*nye+5
         num(1)=fac1
         num(2)=fac1-2*nye-1
         num(3)=fac1-4*nye-2
         num(4)=fac1-6*nye-3
         num(5)=fac1-8*nye-4
         num(6)=fac1-6*nye-4
         num(7)=fac1-4*nye-4
         num(8)=fac1-2*nye-4
         num(9)=fac1-4
         num(10)=fac1-3
         num(11)=fac1-2
         num(12)=fac1-1
         num(13)=fac1-2*nye-2
         num(14)=fac1-4*nye-3
         num(15)=fac1-2*nye-3
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(1,2)=y_coords((iq+2)/2)
       coord(5,1)=x_coords(ip)
       coord(5,2)=y_coords((iq+2)/2)
       coord(9,1)=x_coords(ip+1)
       coord(9,2)=y_coords(iq/2)
     END IF
     coord(3,:)=pt5*(coord(1,:)+coord(5,:))
     coord(7,:)=pt5*(coord(5,:)+coord(9,:))
     coord(11,:)=pt5*(coord(9,:)+coord(1,:))
     coord(2,:)=pt5*(coord(1,:)+coord(3,:))
     coord(4,:)=pt5*(coord(3,:)+coord(5,:))
     coord(6,:)=pt5*(coord(5,:)+coord(7,:))
     coord(8,:)=pt5*(coord(7,:)+coord(9,:))
     coord(10,:)=pt5*(coord(9,:)+coord(11,:))
     coord(12,:)=pt5*(coord(11,:)+coord(1,:))
     coord(15,:)=pt5*(coord(7,:)+coord(11,:))
     coord(14,:)=pt5*(coord(3,:)+coord(7,:))
     coord(13,:)=pt5*(coord(2,:)+coord(15,:))
   CASE DEFAULT
     WRITE(31,'(a)')"Wrong number of nodes for triangular element"
     STOP
   END SELECT
 ELSE
   nye=UBOUND(y_coords,1)-1 ! element = 'quadrilateral'
   IF(dir=='x'.OR.dir=='r')THEN
     iq=(iel+hcount-1)/nxe+1
     ip=iel+hcount-(iq-1)*nxe
   ELSE !dir = 'y'
     ip=(iel+hcount-1)/nye+1 !dir = 'y' iel manipulated
     iq=(iel+hcount)-(ip-1)*nye
   END IF
   SELECT CASE(nod)
   CASE(4)

   dels = 0.d0
   dnods = 0.d0

 do while (dels .le. hels)
  coord(1:2,1)=x_coords(ip)
  coord(3:4,1)=x_coords(ip+1)
  coord(1,2)=y_coords(iq+1)
  coord(2:3,2)=y_coords(iq)
  coord(4,2)=coord(1,2)

  elemidpoint(1,1)=sum(coord(1:4,1))/4.d0
  elemidpoint(1,2)=sum(coord(1:4,2))/4.d0
  rectmidpt(1,1)=sum(x_coords)/nxe-0.5d0
  rectmidpt(1,2)=sum(y_coords)/nye+0.5d0

  dels = dsqrt((rectmidpt(1,1)-elemidpoint(1,1))**2.d0+    &
  (rectmidpt(1,2)-elemidpoint(1,2))**2.d0)

  if (dels .le. hels) then
    hcount=hcount+1
    ip=(iel+hcount-1)/nye+1
    iq=(iel+hcount)-(ip-1)*nye
  endif

  dnods = dsqrt((abs(rectmidpt(1,1)-coord(3,1))+0.5d0) &
  **2.d0+(abs(rectmidpt(1,2)-coord(3,2))+0.5d0)**2.d0)

  if (dnods .le. hels) then
    ncount(1) = ncount(1)+1
  endif

  dnods = dsqrt((abs(rectmidpt(1,1)-coord(1,1))+0.5)  &
  **2.d0+(abs(rectmidpt(1,2)-coord(1,2))+0.5)**2.d0)

    if (dnods .le. hels) then
      ncount(2) = ncount(2)+1
    endif

 enddo
  coord(1:2,1)=x_coords(ip)
  coord(3:4,1)=x_coords(ip+1)
  coord(1,2)=y_coords(iq+1)
  coord(2:3,2)=y_coords(iq)
  coord(4,2)=coord(1,2)

    IF(dir=='x'.OR.dir=='r')THEN
    num(1)=iq*(nxe+1)+ip
    num(2)=(iq-1)*(nxe+1)+ip
    num(3)=num(2)+1
    num(4)=num(1)+1
  ELSE !dir=='y'
    num(1)=(ip-1)*(nye+1)+iq+1-ncount(2)
    num(2)=num(1)-1
    num(3)=ip*(nye+1)+iq-ncount(1)
    num(4)=num(3)+1

  if ((abs(coord(2,2)) .le. 1.d-5) .and. (num(2).ne.1))  then
    ncount(ncount(3))=num(2)!starts at ncount(3)=5
    ncount(3)=ncount(3)+1
  endif

  if ((abs(coord(1,2)+nye) .le. 1.d-5).and. (num(1).ne.nye+1))  then
    ncount(ncount(4))=num(1)
    ncount(4)=ncount(4)+1
  endif

  END IF
   CASE(5)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(2*nxe+1)+ip
       num(2)=(iq-1)*(2*nxe+1)+ip
       num(3)=num(2)+1
       num(4)=num(1)+1
       num(5)=iq*(2*nxe+1)+ip-nxe
     ELSE
       num(1)=(ip-1)*(2*nye+1)+iq+1
       num(2)=num(1)-1
       num(3)=ip*(2*nye+1)+iq
       num(4)=num(3)+1
       num(5)=ip*(2*nye+1)+iq-nye
     END IF
!
     coord(1:2,1)=x_coords(ip)
     coord(3:4,1)=x_coords(ip+1)
     coord(1,2)=y_coords(iq+1)
     coord(2:3,2)=y_coords(iq)
     coord(4,2)=coord(1,2)
     coord(5,:)=0.25_iwp*(coord(1,:)+coord(2,:)+coord(3,:)+coord(4,:))
   CASE(8)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(3*nxe+2)+2*ip-1
       num(2)=iq*(3*nxe+2)+ip-nxe-1
       num(3)=(iq-1)*(3*nxe+2)+2*ip-1
       num(4)=num(3)+1
       num(5)=num(4)+1
       num(6)=num(2)+1
       num(7)=num(1)+2
       num(8)=num(1)+1
     ELSE
       num(1)=(ip-1)*(3*nye+2)+2*iq+1
       num(2)=num(1)-1
       num(3)=num(1)-2
       num(4)=(ip-1)*(3*nye+2)+2*nye+iq+1
       num(5)=ip*(3*nye+2)+2*iq-1
       num(6)=num(5)+1
       num(7)=num(5)+2
       num(8)=num(4)+1
     END IF
!
     coord(1:3,1)=x_coords(ip)
     coord(5:7,1)=x_coords(ip+1)
     coord(4,1)=pt5*(coord(3,1)+coord(5,1))
     coord(8,1)=pt5*(coord(7,1)+coord(1,1))
     coord(1,2)=y_coords(iq+1)
     coord(7:8,2)=y_coords(iq+1)
     coord(3:5,2)=y_coords(iq)
     coord(2,2)=pt5*(coord(1,2)+coord(3,2))
     coord(6,2)=pt5*(coord(5,2)+coord(7,2))
   CASE(9)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(4*nxe+2)+2*ip-1
       num(2)=iq*(4*nxe+2)+2*ip-nxe-4
       num(3)= (iq-1)*(4*nxe+2)+2*ip-1
       num(4)=num(3)+1
       num(5)=num(4)+1
       num(6)=num(2)+2
       num(7)=num(1)+2
       num(8)=num(1)+1
       num(9)=num(2)+1
     ELSE
       num(1)=(ip-1)*2*(2*nye+1)+2*iq+1
       num(2)=num(1)-1
       num(3)=num(1)-2
       num(4)=(ip-1)*2*(2*nye+1)+2*nye+2*iq
       num(5)=ip*2*(2*nye+1)+2*iq-1
       num(6)=num(5)+1
       num(7)=num(5)+2
       num(8)=num(4)+2
       num(9)=num(4)+1
     END IF
!
     coord(1:3,1)=x_coords(ip)
     coord(5:7,1)=x_coords(ip+1)
     coord(4,1)=pt5*(coord(3,1)+coord(5,1))
     coord(8,1)=pt5*(coord(7,1)+coord(1,1))
     coord(1,2)=y_coords(iq+1)
     coord(7:8,2)=y_coords(iq+1)
     coord(3:5,2)=y_coords(iq)
     coord(2,2)=pt5*(coord(1,2)+coord(3,2))
     coord(6,2)=pt5*(coord(5,2)+coord(7,2))
     coord(9,:)=pt5*(coord(4,:)+coord(8,:))
   CASE DEFAULT
     WRITE(31,'(a)')"Wrong number of nodes for quadrilateral element"
     STOP
   END SELECT
 END IF
RETURN
END SUBROUTINE geom_hole

SUBROUTINE mesh_size(element,nod,nels,nn,nxe,nye,nze)
!
!  This subroutine returns the number of elements (nels) and the number
!  of nodes (nn) in a 2-d geometry-created mesh.
!
 IMPLICIT NONE
 CHARACTER(LEN=15),INTENT(IN)::element
 INTEGER,INTENT(IN)::nod,nxe,nye
 INTEGER,INTENT(IN),OPTIONAL::nze
 INTEGER,INTENT(OUT)::nels,nn
 IF(element=="triangle")THEN
   nels=nxe*nye*2
   IF(nod==3)nn=(nxe+1)*(nye+1)
   IF(nod==6)nn=(2*nxe+1)*(2*nye+1)
   IF(nod==10)nn=(3*nxe+1)*(3*nye+1)
   IF(nod==15)nn=(4*nxe+1)*(4*nye+1)
 ELSE IF(element=="quadrilateral")THEN
   nels=nxe*nye
   IF(nod==4)nn=(nxe+1)*(nye+1)
   IF(nod==5)nn=(nxe+1)*(nye+1)+nxe*nye
   IF(nod==8)nn=(2*nxe+1)*(nye+1)+(nxe+1)*nye
   IF(nod==9)nn=(2*nxe+1)*(2*nye+1)
 ELSE IF(element=="hexahedron")THEN
   nels=nxe*nye*nze
   IF(nod==8)nn=(nxe+1)*(nye+1)*(nze+1)
   IF(nod==14)nn=4*nxe*nye*nze+2*(nxe*nye+nye*nze+nze*nxe)+nxe+nye+nze+1
   IF(nod==20)nn=((2*nxe+1)*(nze+1)+(nxe+1)*nze)*(nye+1)+                 &
     (nxe+1)*(nze+1)*nye
 END IF
RETURN
END SUBROUTINE mesh_size

SUBROUTINE mesh_size_hole(element,nod,nels,hels,nn,nxe,nye,nze)
!
!  This subroutine returns the number of elements (nels) and the number
!  of nodes (nn) in a 2-d geometry-created mesh.
!
 IMPLICIT NONE
 CHARACTER(LEN=15),INTENT(IN)::element
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp)::dels,dnods
 INTEGER::nnout, neout,i,j
 INTEGER,INTENT(IN)::nod,nxe,nye,hels
 INTEGER,INTENT(IN),OPTIONAL::nze
 INTEGER,INTENT(OUT)::nels,nn
 IF(element=="triangle")THEN
   nels=nxe*nye*2
   IF(nod==3)nn=(nxe+1)*(nye+1)
   IF(nod==6)nn=(2*nxe+1)*(2*nye+1)
   IF(nod==10)nn=(3*nxe+1)*(3*nye+1)
   IF(nod==15)nn=(4*nxe+1)*(4*nye+1)
 ELSE IF(element=="quadrilateral")THEN
   nnout=0
   neout=0
   do i = 1,nxe
      do j = 1,nxe
        dels = dsqrt((abs(i)-0.5d0)**2.d0+(abs(j)-0.5d0)**2.d0)
      if (dels .le. hels) then
        neout = neout+4
      endif
    enddo
   enddo
   do i = -nxe,nxe
      do j = -nxe,nxe
        dnods =  dsqrt((abs(i)+0.5)**2.d0+(abs(j)+0.5)**2.d0)
        if (dnods .le. hels) then
          nnout = nnout+1
        endif
      enddo
   enddo
    nels=nxe*nye-neout
   IF(nod==4)nn=(nxe+1)*(nye+1)-nnout
   IF(nod==5)nn=(nxe+1)*(nye+1)+nxe*nye
   IF(nod==8)nn=(2*nxe+1)*(nye+1)+(nxe+1)*nye
   IF(nod==9)nn=(2*nxe+1)*(2*nye+1)
 ELSE IF(element=="hexahedron")THEN
   nels=nxe*nye*nze
   IF(nod==8)nn=(nxe+1)*(nye+1)*(nze+1)
   IF(nod==14)nn=4*nxe*nye*nze+2*(nxe*nye+nye*nze+nze*nxe)+nxe+nye+nze+1
   IF(nod==20)nn=((2*nxe+1)*(nze+1)+(nxe+1)*nze)*(nye+1)+                 &
     (nxe+1)*(nze+1)*nye
 END IF
RETURN
END SUBROUTINE mesh_size_hole

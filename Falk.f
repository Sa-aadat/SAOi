! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Some Falk-like dual statements
! SAOi:


       subroutine falk (nprimal,ni,xprimal,xdual,x_l,x_u,x_h,
     &                 f,c,gf,gc,f_a,c_a,acurv,bcurv,acurvn,bcurvn,g,
     &                 kounts,iact,nact,flambda,ns,ncol,nrow,gcvec,
     &                 bvec,ictrl,lctrl,rctrl,cctrl,yhi,xdual_h,
     &                 iuser,luser,cuser,ruser)
!
      implicit         none
      include          'ctrl.h'
      integer          i,j,nprimal,ni,kounts,nact,iact(nact)
      integer          ns,ncol(ns),nrow(ns),icnt
      double precision gcvec(ns),bvec(ns),xdual(ni),xdual_h(ni)
      double precision x_l(nprimal),x_u(nprimal),x_h(nprimal)
      double precision xprimal(nprimal),gf(nprimal),gc(ni,nprimal)
      double precision acurv,bcurv(ni),acurvn (nprimal),yhi
      double precision bcurvn(ni,nprimal),c(ni),f,f_a,c_a(ni)
      double precision temp1,temp2,beta,flambda,g(ni)
      include          'ctrl_get.inc'

!  calculate the Falk dual - sparse implementation
      if (structure.eq.2) then
        if (ifalk.eq.0) then
          call falk_dq_pseudo (nprimal,ni,xprimal,xdual,x_l,x_u,
     &                         x_h,f,c,gf,f_a,c_a,acurv,bcurv,
     &                         acurvn,g,kounts,iact,nact,flambda,
     &                         ns,ncol,nrow,gcvec,bvec,ictrl,lctrl,
     &                         rctrl,cctrl,yhi,xdual_h,
     &                         iuser,luser,cuser,ruser)
        elseif (ifalk.eq.1) then
          call falk_dq_pseudoq (nprimal,ni,xprimal,xdual,x_l,x_u,
     &                          x_h,f,c,gf,f_a,c_a,acurv,bcurv,
     &                          acurvn,g,kounts,iact,nact,flambda,
     &                          ns,ncol,nrow,gcvec,bvec,ictrl,lctrl,
     &                          rctrl,cctrl,yhi,xdual_h,
     &                          iuser,luser,cuser,ruser)
        else
          stop ' ifalk not defined in falk'
        endif

!  calculate the Falk dual - dense implementation
      elseif (structure.eq.1) then
        if (ifalk.eq.0) then
          call falk_dq (nprimal,ni,xprimal,xdual,x_l,x_u,
     &                  x_h,f,c,gf,gc,f_a,c_a,acurv,bcurv,
     &                  acurvn,bcurvn,g,kounts,iact,nact,
     &                  flambda,ictrl,lctrl,rctrl,cctrl,yhi,xdual_h,
     &                  iuser,luser,cuser,ruser)
        elseif (ifalk.eq.1) then
          call falk_dqq (nprimal,ni,xprimal,xdual,x_l,x_u,
     &                   x_h,f,c,gf,gc,f_a,c_a,acurv,bcurv,
     &                   acurvn,bcurvn,g,kounts,iact,nact,
     &                   flambda,ictrl,lctrl,rctrl,cctrl,yhi,xdual_h,
     &                   iuser,luser,cuser,ruser)
        else
          stop ' ifalk not defined in falk'
        endif
      endif
!
      return
      end
!----------------------------------------------------------------------c
      subroutine falk_dq (nprimal,ni,xprimal,xdual,x_l,x_u,x_h,
     &                    f,c,gf,gc,f_a,c_a,acurv,bcurv,acurvn,
     &                    bcurvn,g,kounts,iact,nact,flambda,
     &                    ictrl,lctrl,rctrl,cctrl,yhi,xdual_h,
     &                    iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,nprimal,ni,kounts,nact,iact(nact),ilo,ihi
      integer          imeth      
      double precision xdual(ni),x_l(nprimal),x_u(nprimal),zero,one,two
      double precision xprimal(nprimal),gf(nprimal),gc(ni,nprimal)
      double precision acurv,bcurv(ni),acurvn (nprimal),gcji,bcji
      double precision bcurvn(ni,nprimal),c(ni),f,f_a,c_a(ni)
      double precision temp1,temp2,beta,flambda,g(ni),xdual_h(ni)
      double precision x_h(nprimal),xdualj,y,sum1,yhi
      double precision temp1n(nprimal),temp2n(nprimal)
      data             zero /0.d0/,one /1.d0/,two /2.d0/
      data             imeth /1/
      include          'ctrl_get.inc'

!  construct the bounded dual
      if (imeth.eq.1) then
        do i=1,nprimal
          temp1 = 0.d0
          temp2 = 0.d0
          do k = 1,nact
            j = iact(k)
            xdualj=xdual(j)
            if (xdualj.gt.0.d0) then
              gcji=gc(j,i)
              if (gcji.ne.0.d0) then
                temp1 = temp1 + xdualj*gcji
              endif
              bcji=bcurvn(j,i)
              if (bcji.gt.0.d0) then
                temp2 = temp2 + xdualj*bcji
              endif
            endif
          enddo
!          
          beta = x_h(i)-(gf(i) + temp1)/max(1.d-6,(acurvn(i)+temp2))
          xprimal(i) = min(x_u(i),max(x_l(i),beta))
        enddo
!      
      elseif (imeth.eq.2) then ! for the GPU; enforce BLAS !
        call dcopy (nprimal,gf,1,temp1n,1)
        call dgemv ('t',ni,nprimal,1.0d0,gc,ni,xdual,1,1.d0,temp1n,1)
        call dcopy (nprimal,acurvn,1,temp2n,1)
        call dgemv ('t',ni,nprimal,1.0d0,bcurvn,ni,xdual,1,1.d0,
     &              temp2n,1)
!
        do i=1,nprimal
          beta = x_h(i)-temp1n(i)/max(1.d-6,temp2n(i))
          xprimal(i) = min(x_u(i),max(x_l(i),beta))
        enddo
      endif  

!  evaluate the single relaxation variable y if needed
      if (relax) then
        sum1 = zero
        do k = 1,nact
          j = iact(k)
          sum1=sum1+xdual(j)
        enddo
        y = (sum1-pen1)/pen2       !   we desire biglam > pen1   !
        y = min(ymax,max(zero,y))
      else
        y = zero
      endif

!  have updated xprimal; now get function value to maximize the dual
      kounts = kounts + 1
      call fun_a (nprimal,xprimal,f_a,x_h,acurv,acurvn,f,gf)
      call conin_a (nprimal,ni,xprimal,c_a,iact,nact,x_h,bcurv,
     &              bcurvn,c,gc)
      flambda = -f_a
      do k = 1,nact
        j = iact(k)
        flambda = flambda - xdual(j)*c_a(j)
      enddo

!  compute gradient g
      do k = 1,nact
        j = iact(k)
        g(j) = -c_a(j)
      enddo

!  do some corrections for relaxation
      if (relax) then
        flambda = flambda - pen1*y - pen2*y**2/two
        do k = 1,nact
          j = iact(k)
          flambda = flambda + xdual(j)*y
          g(j) = g(j) + y
        enddo
        yhi = y
      endif
!
      return
      end
!----------------------------------------------------------------------c
      subroutine falk_dq_pseudo (nprimal,ni,xprimal,xdual,x_l,x_u,
     &                           x_h,f,c,gf,f_a,c_a,acurv,bcurv,
     &                           acurvn,g,kounts,iact,nact,flambda,
     &                           ns,ncol,nrow,gcvec,bvec,ictrl,lctrl,
     &                           rctrl,cctrl,yhi,xdual_h,
     &                           iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,nprimal,ni,kounts,nact,iact(nact)
      integer          ij,i1,ns,ncol(ns),nrow(ns),icnt,ilo,ihi,nga,nga2
      double precision gcvec(ns),bvec(ns)
      double precision xdual(ni),x_l(nprimal),x_u(nprimal),zero,one,two
      double precision xprimal(nprimal),gf(nprimal),xdual_h(ni),glo,ghi
      double precision acurv,bcurv(ni),acurvn (nprimal),gcji,bcji
      double precision c(ni),f,f_a,c_a(ni),xlo,xhi
      double precision temp1(nprimal),temp2(nprimal)
      double precision beta,flambda,g(ni)
      double precision x_h(nprimal),xdualj,y,sum1,yhi
      data             zero /0.d0/,one /1.d0/,two /2.d0/
      include          'ctrl_get.inc'

!  construct the bounded dual,looping over ndual
      do i=1,nprimal
        temp1(i)=gf(i)
        temp2(i)=acurvn(i)
      enddo
!
      do j=1,ni
        if (xdual(j).lt.1.d-6) xdual(j) = 0.d0 ! = max(xdual(j),0.d0)
      enddo
!
      do ij=1,ns
        j=ncol(ij)
        if (xdual(j).gt.0.d0) then
          xdualj=xdual(j)
          i=nrow(ij)
          gcji=gcvec(ij)
          bcji=bvec(ij)
          temp1(i) = temp1(i) + xdualj*gcji
          temp2(i) = temp2(i) + xdualj*bcji
        endif
      enddo
!
      do i=1,nprimal
        beta = x_h(i)-temp1(i)/max(1.d-6,temp2(i))
        xprimal(i) = min(x_u(i),max(x_l(i),beta))
      enddo

!  evaluate the relaxation variables if needed
      if (relax) then
        sum1 = zero
        do k = 1,nact
          j = iact(k)
          sum1=sum1+xdual(j)
        enddo
        y = (sum1-pen1)/pen2       !   we desire biglam > pen1   !
        y = min(ymax,max(zero,y))
      else
        y = zero
      endif

!  have updated xprimal; now get function value to maximize the dual
      kounts = kounts + 1
      call fun_a (nprimal,xprimal,f_a,x_h,acurv,acurvn,f,gf)
      call conin_a_sprse (nprimal,ni,xprimal,c_a,iact,nact,x_h,
     &                    bcurv,c,ns,ncol,nrow,gcvec,bvec)
      flambda = -f_a
      do k = 1,nact
        j = iact(k)
        flambda = flambda - xdual(j)*c_a(j)
      enddo

!  compute gradient g
      do k = 1,nact
        j = iact(k)
        g(j) = -c_a(j)
      enddo

!  do some corrections for relaxation
      if (relax) then
        flambda = flambda - pen1*y - pen2*y**2/two
        do k = 1,nact
          j = iact(k)
          flambda = flambda + xdual(j)*y
          g(j) = g(j) + y
        enddo
        yhi = y
      endif

!
      ilo = 0
      ihi = 0
      do i=1,nprimal
        if (xprimal(i).eq.x_l(i)) ilo = ilo + 1
        if (xprimal(i).eq.x_h(i)) ihi = ihi + 1
      enddo
!
      glo = 1.d99
      ghi = -1.d99
      nga = 0
      nga2 = 0
      xlo = 1.d99
      xhi = -1.d99
      do k = 1,nact
        j = iact(k)
        if (g(j).gt.ghi) ghi = g(j)
        if (g(j).lt.glo) glo = g(j)
        if (g(j).gt.0.d0) nga = nga + 1
        if (g(j).gt.feaslim) nga2 = nga2 + 1
        if (xdual(j).gt.xhi) xhi = xdual(j)
        if (xdual(j).lt.xlo) xlo = xdual(j)
      enddo
!       write(*,1000)nprimal,nprimal-ilo-ihi,ilo,ihi,nact,nga,nga2,
!      &             nga-nga2,glo,ghi,xlo,xhi,-flambda,f
!
      return
 1000 format (8i5,6es16.4)
      end
!----------------------------------------------------------------------c
      subroutine falk_dq_full (nprimal,ni,xprimal,xdual,x_l,x_u,x_h,
     &                         f,c,gf,gc,f_a,c_a,acurv,bcurv,acurvn,
     &                         bcurvn,g,kounts,iact,nact,flambda,
     &                         ictrl,lctrl,rctrl,cctrl,yhi,xdual_h,
     &                         iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,nprimal,ni,kounts,nact,iact(nact)
      double precision xdual(ni),x_l(nprimal),x_u(nprimal),zero,one,two
      double precision xprimal(nprimal),gf(nprimal),gc(ni,nprimal)
      double precision acurv,bcurv(ni),acurvn (nprimal),gcji,bcji
      double precision bcurvn(ni,nprimal),c(ni),f,f_a,c_a(ni)
      double precision temp1,temp2,beta,flambda,g(ni),xdual_h(ni)
      double precision x_h(nprimal),xdualj,y(ni),sum1,yhi,temp3
      data             zero /0.d0/,one /1.d0/,two /2.d0/
      include          'ctrl_get.inc'

!  construct the bounded dual
      do i=1,nprimal
        temp1 = 0.d0
        temp2 = 0.d0
        do k = 1,nact
          j = iact(k)
          xdualj=xdual(j)
          if (xdualj.gt.0.d0) then
            gcji=gc(j,i)
            if (gcji.ne.0.d0) then
              temp1 = temp1 + xdualj*gcji
            endif
            bcji=bcurvn(j,i)
            if (bcji.gt.0.d0) then 
              temp2 = temp2 + xdualj*bcji ! we do not allow concave functions !
            endif
          endif
        enddo
        temp3 = acurvn(i)+temp2
        beta = x_h(i)-(gf(i) + temp1)/max(1.d-6,temp3)
        xprimal(i) = min(x_u(i),max(x_l(i),beta))
      enddo

 !  evaluate the ni relaxation variables y(j) if needed
      if (relax) then
        do j=1,ni
          y(j) = (xdual(j)-pen1)/pen2
          y(j) = min(ymax,max(zero,y(j)))
        enddo
      endif

!  have updated xprimal; now get function value to maximize the dual
      kounts = kounts + 1
      call fun_a (nprimal,xprimal,f_a,x_h,acurv,acurvn,f,gf)
      call conin_a (nprimal,ni,xprimal,c_a,iact,nact,x_h,bcurv,
     &              bcurvn,c,gc)
      flambda = -f_a
      do k = 1,nact
        j = iact(k)
        flambda = flambda - xdual(j)*c_a(j)
      enddo

!  compute gradient g
      do k = 1,nact
        j = iact(k)
        g(j) = -c_a(j)
      enddo

!  do some corrections for relaxation
      yhi=0.d0
      if (relax) then
        do k = 1,nact
          j = iact(k)
          flambda = flambda - pen1*y(j) - pen2*y(j)**2/two
          flambda = flambda + xdual(j)*y(j)
          g(j) = g(j) + y(j)
          yhi = max(yhi,y(j))
        enddo
      endif
!
      return
      end
!----------------------------------------------------------------------c
      subroutine falk_dqq (nprimal,ni,xprimal,xdual,x_l,x_u,x_h,
     &                     f,c,gf,gc,f_a,c_a,acurv,bcurv,acurvn,
     &                     bcurvn,g,kounts,iact,nact,flambda,
     &                     ictrl,lctrl,rctrl,cctrl,yhi,xdual_h,
     &                     iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,nprimal,ni,kounts,nact,iact(nact)
      double precision xdual(ni),x_l(nprimal),x_u(nprimal),zero,one,two
      double precision xprimal(nprimal),gf(nprimal),gc(ni,nprimal)
      double precision acurv,bcurv(ni),acurvn (nprimal),gcji,bcji
      double precision bcurvn(ni,nprimal),c(ni),f,f_a,c_a(ni)
      double precision temp1,temp2,beta,flambda,g(ni),xdual_h(ni)
      double precision x_h(nprimal),xdualjt,xdualjb,y,sum1,yhi,temp3
      data             zero /0.d0/,one /1.d0/,two /2.d0/
      include          'ctrl_get.inc'

!  construct the bounded dual
      do i=1,nprimal
        temp1 = 0.d0
        temp2 = 0.d0
        do k = 1,nact
          j = iact(k)
          xdualjt=xdual(j)
          xdualjb=xdual_h(j)
          if (xdualjt.ne.0.d0) then
            gcji=gc(j,i)
            if (gcji.ne.0.d0) then
              temp1 = temp1 + xdualjt*gcji
            endif
          endif
          if (xdualjb.ne.0.d0) then
            bcji=bcurvn(j,i)
            !if (bcji.gt.0.d0) then
              temp2 = temp2 + xdualjb*bcji ! we allow concave functions !
            !endif
          endif
        enddo
        temp3 = acurvn(i)+temp2
        beta = x_h(i)-(gf(i) + temp1)/max(1.d-6,temp3)
        xprimal(i) = min(x_u(i),max(x_l(i),beta))
      enddo

!  evaluate the single relaxation variable y if needed
      if (relax) then
        sum1 = zero
        do k = 1,nact
          j = iact(k)
          sum1=sum1+xdual(j)
        enddo
        y = (sum1-pen1)/pen2       !   we desire biglam > pen1   !
        y = min(ymax,max(zero,y))
      else
        y = zero
      endif

!  have updated xprimal; now get function value to maximize the dual
      kounts = kounts + 1
      call fun_nabla2q (nprimal,ni,xprimal,f_a,x_h,acurv,acurvn,f,gf,
     &                  bcurv,bcurvn,xdual_h)
      call conin_nabla2q (nprimal,ni,xprimal,c_a,iact,nact,x_h,c,gc)

      flambda = -f_a
      do k = 1,nact
        j = iact(k)
        flambda = flambda - xdual(j)*c_a(j)
      enddo

!  compute gradient g
      do k = 1,nact
        j = iact(k)
        g(j) = -c_a(j)
      enddo

!  do some corrections for relaxation
      if (relax) then
        flambda = flambda - pen1*y - pen2*y**2/two
        do k = 1,nact
          j = iact(k)
          flambda = flambda + xdual(j)*y
          g(j) = g(j) + y
        enddo
        yhi = y
      endif
!
      return
      end
!----------------------------------------------------------------------c
      subroutine falk_dq_pseudoq (nprimal,ni,xprimal,xdual,x_l,x_u,
     &                           x_h,f,c,gf,f_a,c_a,acurv,bcurv,
     &                           acurvn,g,kounts,iact,nact,flambda,
     &                           ns,ncol,nrow,gcvec,bvec,ictrl,lctrl,
     &                           rctrl,cctrl,yhi,xdual_h,
     &                           iuser,luser,cuser,ruser)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,nprimal,ni,kounts,nact,iact(nact)
      integer          ij,i1
      integer          ns,ncol(ns),nrow(ns),icnt
      double precision gcvec(ns),bvec(ns)
      double precision xdual(ni),x_l(nprimal),x_u(nprimal),zero,one,two
      double precision xprimal(nprimal),gf(nprimal),xdual_h(ni)
      double precision acurv,bcurv(ni),acurvn (nprimal),gcji,bcji
      double precision c(ni),f,f_a,c_a(ni)
      double precision temp1(nprimal),temp2(nprimal)
      double precision beta,flambda,g(ni)
      double precision x_h(nprimal),xdualt,xdualb,y,sum1,yhi
      data             zero /0.d0/,one /1.d0/,two /2.d0/
      include          'ctrl_get.inc'

!  construct the bounded dual,looping over ndual
      do i=1,nprimal
        temp1(i)=gf(i)
        temp2(i)=acurvn(i)
      enddo
!
      do ij=1,ns
        j=ncol(ij)
        xdualt=xdual(j)
        xdualb=xdual_h(j)
        i=nrow(ij)
        gcji=gcvec(ij)
        bcji=bvec(ij)
        if (xdualt.ne.0.d0) temp1(i) = temp1(i) + xdualt*gcji
        if (xdualb.ne.0.d0) temp2(i) = temp2(i) + xdualb*bcji
      enddo
!
      do i=1,nprimal
        beta = x_h(i)-temp1(i)/max(1.d-6,temp2(i))
        xprimal(i) = min(x_u(i),max(x_l(i),beta))
      enddo

!  evaluate the relaxation variables if needed
      if (relax) then
        sum1 = zero
        do k = 1,nact
          j = iact(k)
          sum1=sum1+xdual(j)
        enddo
        y = (sum1-pen1)/pen2       !   we desire biglam > pen1   !
        y = min(ymax,max(zero,y))
      else
        y = zero
      endif

!  have updated xprimal; now get function value to maximize the dual
      kounts = kounts + 1
      call fun_nabla2sq (nprimal,ni,xprimal,f_a,x_h,acurv,acurvn,f,gf,
     &                         xdual_h,ns,ncol,nrow,bvec)
      call conin_nabla2_sprseq (nprimal,ni,xprimal,c_a,iact,nact,x_h,
     &                               c,ns,ncol,nrow,gcvec)
      flambda = -f_a
      do k = 1,nact
        j = iact(k)
        flambda = flambda - xdual(j)*c_a(j)
      enddo

!  compute gradient g
      do k = 1,nact
        j = iact(k)
        g(j) = -c_a(j)
      enddo

!  do some corrections for relaxation
      if (relax) then
        flambda = flambda - pen1*y - pen2*y**2/two
        do k = 1,nact
          j = iact(k)
          flambda = flambda + xdual(j)*y
          g(j) = g(j) + y
        enddo
        yhi = y
      endif

!
      return
      end
!----------------------------------------------------------------------c

! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold,S tellenbosch
! SAOi: The diagonal quadratic SAOI approximations
! SAOi: To-do: blasify everything
! SAOi:


      subroutine fun_a (n,x,f_a,x_h,acurv,acurvn,f,gf)
      implicit         none
      integer          i,n
      double precision x(n),x_h(n),f_a
      double precision acurv,acurvn(n),f,gf(n)
!
      call fun_nabla2 (n,x,f_a,x_h,acurv,acurvn,f,gf)  ! replace with conditional library call
!
      return
      end subroutine fun_a
!-----------------------------------------------------------------------
      subroutine conin_a (n,ni,x,c_a,iact,nact,x_h,bcurv,bcurvn,c,gc)
      implicit         none
      integer          i,n,ni
      integer          nact,iact(nact)
      double precision x(n),c_a(ni),x_h(n),c(ni),gc(ni*n)
      double precision bcurv(ni),bcurvn(ni*n)
!
      call conin_nabla2 (n,ni,x,c_a,iact,nact,x_h,bcurv,bcurvn,c,gc)
!
      return
      end subroutine conin_a
!-----------------------------------------------------------------------
      subroutine gradf_a (n,x,f_a,x_h,acurv,acurvn,f,gf,gf_a)
      implicit         none
      integer          i,n
      double precision x(n),x_h(n),delx(n),f_a
      double precision acurv,acurvn(n),f,gf(n),gf_a(n)
!
      call gradf_nabla2 (n,x,f_a,x_h,acurv,acurvn,f,gf,gf_a)
!
      return
      end subroutine gradf_a
!-----------------------------------------------------------------------
      subroutine gradc_a (n,ni,x,c_a,iact,nact,x_h,bcurv,bcurvn,
     &                    c,gc,gc_a)
      implicit         none
      integer          i,n,ni
      integer          nact,iact(nact)
      double precision x(n),c_a(ni),x_h(n),c(ni),gc(ni*n),delx(n)
      double precision bcurv(ni),bcurvn(ni*n),gc_a(ni*n)
!
      call gradc_nabla2 (n,ni,x,c_a,iact,nact,x_h,bcurv,bcurvn,
     &                   c,gc,gc_a)
!
      return
      end subroutine gradc_a
!-----------------------------------------------------------------------
      subroutine gradh_a (n,nq,x,h_a,x_h,ccurv,ccurvn,h,gh,gh_a)
      implicit         none
      integer          i,n,nq
      double precision x(n),h_a(nq),x_h(n),h(nq),gh(nq,n),delx(n)
      double precision ccurv(n),ccurvn(nq,n),gh_a(nq,n)
!
      call gradh_nabla2 (n,nq,x,h_a,x_h,ccurv,ccurvn,h,gh,gh_a)
!
      return
      end subroutine gradh_a
!----------------------------------------------------------------------c
      subroutine fun_nabla2 (n,x,f_a,x_h,acurv,acurvn,f,gf)
      implicit         none
      integer          i,n
      double precision x(n),delx(n),x_h(n),f,gf(n),f_a,dp,temp,ddot,xax
      double precision acurv,acurvn(n),delx2(n)
!
      do i=1,n
        temp     = x(i)-x_h(i)
        delx(i)  = temp
        delx2(i) = temp**2
      end do
!
      xax = ddot(n,delx2,1,acurvn,1)
      dp  = ddot(n,delx,1,gf,1)
!
      f_a = f+dp+xax/2.d0
!
      return
      end subroutine fun_nabla2
!-----------------------------------------------------------------------
      subroutine conin_nabla2 (n,ni,x,c_a,iact,nact,x_h,bcurv,bcurvn,
     &                        c,gc)
      implicit         none
      integer          i,n,ni
      integer          nact,iact(nact)
      double precision x(n),delx(n),c_a(ni),x_h(n),c(ni),gc(ni*n),temp
      double precision bcurv(ni),bcurvn(ni*n),delx2(n)
!
      do i=1,n
        temp    = x(i)-x_h(i)
        delx(i) = temp
        delx2(i)= temp**2
      end do
!
!  uses BLAS since Version 0.1.5
      call dcopy (ni,c,1,c_a,1)
      call dgemv ('n',ni,n,1.0d0,gc,ni,delx,1,1.d0,c_a,1)
      call dgemv ('n',ni,n,0.5d0,bcurvn,ni,delx2,1,1.d0,c_a,1)
!
      return
      end subroutine conin_nabla2
!-----------------------------------------------------------------------
      subroutine gradf_nabla2 (n,x,f_a,x_h,acurv,acurvn,f,gf,gf_a)
      implicit         none
      integer          i,n
      double precision x(n),x_h(n),delx(n),f_a
      double precision acurv,acurvn(n),f,gf(n),gf_a(n)
!
      do i=1,n
        delx(i)=(x(i)-x_h(i))
      end do
      do i=1,n
        gf_a(i)=gf(i)+acurvn(i)*delx(i)
      end do
!
      return
      end subroutine gradf_nabla2
!-----------------------------------------------------------------------
      subroutine gradc_nabla2 (n,ni,x,c_a,iact,nact,x_h,
     &                         bcurv,bcurvn,c,gc,gc_a)
      implicit         none
      integer          i,j,k,n,ni,nact,iact(nact)
      double precision x(n),c_a(ni),x_h(n),c(ni),gc(ni*n),delx(n)
      double precision bcurv(ni),bcurvn(ni,n),gc_a(ni*n)
!
      do i=1,n
        delx(i)=(x(i)-x_h(i))
      end do
!
      do i=1,n
        do j=1,ni
          k=(i-1)*ni+j
          gc_a(k)=gc(k)+bcurvn(j,i)*delx(i)
        end do
      end do
!
      return
      end subroutine gradc_nabla2
!-----------------------------------------------------------------------
      subroutine conin_a_sprse (n,ni,x,c_a,iact,nact,x_h,bcurv,
     &                          c,ns,ncol,nrow,gcvec,bvec)
      implicit         none
      integer          i,n,ni
      integer          nact,iact(nact)
      integer          ns,ncol(ns),nrow(ns),icnt
      double precision gcvec(ns),bvec(ns)
      double precision x(n),c_a(ni),x_h(n),c(ni)
      double precision bcurv(ni)
!
      call conin_nabla2_sprse (n,ni,x,c_a,iact,nact,x_h,bcurv,
     &                         c,ns,ncol,nrow,gcvec,bvec)
!
      return
      end subroutine conin_a_sprse
!-----------------------------------------------------------------------
      subroutine conin_nabla2_sprse (n,ni,x,c_a,iact,nact,x_h,bcurv,
     &                               c,ns,ncol,nrow,gcvec,bvec)
      implicit         none
      integer          i,n,ni,ij,j
      integer          nact,iact(nact)
      integer          ns,ncol(ns),nrow(ns),icnt
      double precision gcvec(ns),bvec(ns),temp
      double precision x(n),delx(n),c_a(ni),x_h(n),c(ni)
      double precision bcurv(ni)
!
      do i=1,n
        temp    = x(i)-x_h(i)
        delx(i) = temp
      end do
!
      call dcopy (ni,c,1,c_a,1)
      do ij=1,ns
        j=ncol(ij)
        i=nrow(ij)
        c_a(j)=c_a(j)+delx(i)*gcvec(ij)
      enddo
!
      do ij=1,ns
        j=ncol(ij)
        i=nrow(ij)
        c_a(j)=c_a(j)+delx(i)**2*bvec(ij)/2.d0
      enddo
!
      return
      end subroutine conin_nabla2_sprse
!----------------------------------------------------------------------c
      subroutine fun_nabla2q (n,ni,x,f_a,x_h,acurv,acurvn,f,gf,bcurv,
     &                        bcurvn,xlam_h)
      implicit         none
      integer          i,j,n,ni
      double precision x(n),delx(n),x_h(n),f,gf(n),f_a,dp,temp,ddot,xax
      double precision acurv,acurvn(n),delx2(n),acurvnn(n)
      double precision bcurv(ni),bcurvn(ni,n),xlam_h(ni)
!
      do i=1,n
        temp     = x(i)-x_h(i)
        delx(i)  = temp
        delx2(i) = temp**2
      end do
!
      call dcopy(n,acurvn,1,acurvnn,1)
      do i=1,n
        do j=1,ni
          acurvnn(i) = acurvnn(i) + xlam_h(j)*bcurvn(j,i)
        enddo
      enddo
!
      xax = ddot(n,delx2,1,acurvnn,1)
      dp  = ddot(n,delx,1,gf,1)
!
      f_a = f+dp+xax/2.d0
!
      return
      end subroutine fun_nabla2q
!-----------------------------------------------------------------------
      subroutine conin_nabla2q (n,ni,x,c_a,iact,nact,x_h,c,gc)
      implicit         none
      integer          i,n,ni
      integer          nact,iact(nact)
      double precision x(n),delx(n),c_a(ni),x_h(n),c(ni),gc(ni*n),temp
!
      do i=1,n
        temp    = x(i)-x_h(i)
        delx(i) = temp
      end do
!
      call dcopy (ni,c,1,c_a,1)
      call dgemv ('n',ni,n,1.0d0,gc,ni,delx,1,1.d0,c_a,1)
!
      return
      end subroutine conin_nabla2q
!----------------------------------------------------------------------c
      subroutine fun_nabla2sq (n,ni,x,f_a,x_h,acurv,acurvn,f,gf,
     &                         xlam_h,ns,ncol,nrow,bvec)
      implicit         none
      integer          i,j,n,ni,ij
      integer          ns,ncol(ns),nrow(ns),icnt
      double precision bvec(ns)
      double precision x(n),delx(n),x_h(n),f,gf(n),f_a,dp,temp,ddot,xax
      double precision acurv,acurvn(n),delx2(n),acurvnn(n)
      double precision xlam_h(ni)
!
      do i=1,n
        temp     = x(i)-x_h(i)
        delx(i)  = temp
        delx2(i) = temp**2
      end do
!
      call dcopy(n,acurvn,1,acurvnn,1)
!
      do ij=1,ns
        j=ncol(ij)
        i=nrow(ij)
        acurvnn(i)=acurvnn(i) + xlam_h(j)*bvec(ij)
      enddo
!
      xax = ddot(n,delx2,1,acurvnn,1)
      dp  = ddot(n,delx,1,gf,1)
!
      f_a = f+dp+xax/2.d0
!
      return
      end subroutine fun_nabla2sq
!-----------------------------------------------------------------------
      subroutine conin_nabla2_sprseq (n,ni,x,c_a,iact,nact,x_h,
     &                               c,ns,ncol,nrow,gcvec)
      implicit         none
      integer          i,n,ni,ij,j
      integer          nact,iact(nact)
      integer          ns,ncol(ns),nrow(ns),icnt
      double precision gcvec(ns),temp
      double precision x(n),delx(n),c_a(ni),x_h(n),c(ni)
!
      do i=1,n
        temp    = x(i)-x_h(i)
        delx(i) = temp
      end do
!
      call dcopy (ni,c,1,c_a,1)
      do ij=1,ns
        j=ncol(ij)
        i=nrow(ij)
        c_a(j)=c_a(j)+delx(i)*gcvec(ij)
      enddo
!
      return
      end subroutine conin_nabla2_sprseq
!-----------------------------------------------------------------------
      subroutine coneq_a (n,nq,x,h_a,x_h,ccurv,ccurvn,h,gh)
      implicit         none
      integer          n,nq
      double precision x(n),h_a(nq),x_h(n),h(nq),gh(nq,n)
      double precision ccurv(nq),ccurvn(nq,n)
!
      call coneq_nabla2 (n,nq,x,h_a,x_h,ccurv,ccurvn,h,gh)
!
      return
      end subroutine coneq_a
!-----------------------------------------------------------------------
      subroutine coneq_nabla2 (n,nq,x,h_a,x_h,ccurv,ccurvn,h,gh)
      implicit         none
      integer          n,nq,i,j,k
      double precision x(n),h_a(nq),delx(n),x_h(n),h(nq),gh(nq,n)
      double precision ccurv(nq),ccurvn(nq,n),dp,dp1
!
      do i=1,n
        delx(i)=(x(i)-x_h(i))
      end do
!
      do k=1,nq
        dp1=0.d0
        do i=1,n
          dp1=dp1+delx(i)**2*ccurvn(k,i) ! recode with blas
        end do
!
        dp=0.d0
        do i=1,n
          dp=dp+delx(i)*gh(k,i)
        end do
        h_a(k)=h(k)+dp+0.5d0*dp1
      end do
!
      return
      end subroutine coneq_nabla2
!-----------------------------------------------------------------------
      subroutine gradh_nabla2 (n,nq,x,h_a,x_h,ccurv,ccurvn,h,gh,gh_a)
      implicit         none
      integer          n,nq,i,j,k
      double precision x(n),h_a(nq),x_h(n),h(nq),gh(nq,n),delx(n)
      double precision ccurv(n),ccurvn(nq,n),gh_a(nq,n)
!
      do i=1,n
        delx(i)=(x(i)-x_h(i))
      end do
!
      do j=1,n
        do k=1,nq
          gh_a(k,j)=gh(k,j)+ccurvn(k,j)*delx(j)
        end do
      end do
!
      return
      end subroutine gradh_nabla2
!-----------------------------------------------------------------------

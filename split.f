! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Initiate either a single local search, or a global search pattern
! SAOi:


      subroutine SAOi_split (n,m,ni,ne,x,x_l,x_u,timef,timeg,times,
     &                  ictrl,lctrl,rctrl,cctrl,iuser,luser,cuser,
     &                  ruser,nnz,nnzh,eqn,lin,mxloop,time0,ierr,iseed,
     &                  shift,lambda,sloop,z)
      implicit         none
      include          'ctrl.h'
      logical          eqn(*),lin(*)
      logical          warn,feasible
      integer          i,n,m,ni,ne,it,ir,n5,kl,ll,nfe,nge,ierr,iname
      integer          nnz,nnzh,mxloop,nw,iseed,kltt,lltt,sloop,nels
      double precision x(n),x_l(n),x_u(n),xbest(n),f,v,fopt,vopt,p,rdm
      double precision conv,timef,timeg,times,xkkt,xkktopt,lambda(m)
      double precision time0,timei,SAOi_seconds,drandu1,shift(n),z(n)
      external         conv,SAOi_seconds
      include          'ctrl_get.inc'

!  zero the times
      timef         = 0.d0
      timeg         = 0.d0
      times         = 0.d0

!  zero the counters
      nfe           = 0
      nge           = 0

!  initialize
      n5            = min(n,5)
      warn          = .false.
      nw            = m+2*n
      iname         = 0
!     iseed         = 1         ! for repeatability tests !

!  call the secondary driver
      call SAOi_splita (n,m,ni,ne,x,x_l,x_u,timef,timeg,times,f,v,
     &                  xkkt,kl,ll,ictrl,lctrl,rctrl,cctrl,warn,nfe,nge,
     &                  feasible,iuser,luser,cuser,ruser,nnz,nnzh,eqn,
     &                  lin,nw,mxloop,ierr,time0,iseed,shift,lambda,
     &                  iname,sloop,z)
!
      kltt = kl
      lltt = ll

!  return if the problem is assumed to be unimodal
      if (.not.multimodal.or.itglobalmax.le.1) then
        if (warn) write (6,10)
        open  (16,file='test.out',status='old',position='append')
        timei = SAOi_seconds()
        write (16,9616) cname1,n,ni,ne,subsolver,f,v,xkkt,kl,ll,
     &                  timei-time0
        close (16)
        call mesh_case(1)
        call matid_end(0,x,n,nels)
        return
      endif

!  return if there are errors
      if (ierr.lt.0) return

!  the problem is assumed to be multimodal; write a header
      iname = 1
      write(*,6000)
      
      lambda = 0.d0
      z      = 0.d0

!  initialize
      fopt = f
      do i=1,n
        xbest(i)=x(i)
      enddo
      it = 1
      ir = 1
      p = 1.0d0-conv(it,ir)
      xkktopt = xkkt

      call matid_end(it-1,x,n,nels)

!  output the results of the first global iteration
      write(*,1000)it,ir,p,f,v,fopt

      do it = 2,itglobalmax

        write  (8,100) it
        write  (9,100) it
        write (10,100) it
        write (11,100) it
        if (debug) write (12,100) it
        write (13,100) it
        if (warn) write (14,100) it
        write (15,100) it

!  construct a random starting point
        do i=1,nels
          rdm = drandu1(iseed)
          x(i)=(x_u(i)-x_l(i))*(rdm)+x_l(i)
        enddo

!  alternatively, use some clever strategy to prescribe a (hopefully) good starting guess
!
!       if (it.eq.2) then                 ! conditionally initialize x here
!         ...
!       elseif (it.eq.3) then             ! and here
!         ...
!       elseif (it.eq.4) then             ! and here
!         ...
!       else                              ! but not here
!         do i=1,n
!           rdm = drandu1(iseed)
!           x(i)=(x_u(i)-x_l(i))*(rdm)+x_l(i)
!         enddo
!       endif

!  even better, use some algorithm with a global search capability here to calculate a good starting point !
!
!       call DE_PSO (n,m,ni,ne,x,x_l,x_u,timef,timeg,times,f,v,
!      &             xkkt,kl,ll,ictrl,lctrl,rctrl,cctrl,warn,nfe,
!      &             nge,feasible,iuser,luser,cuser,ruser,nnz,nnzh,
!      &             eqn,lin,nw,mxloop,ierr)
!

!  call the secondary driver
        call SAOi_splita (n,m,ni,ne,x,x_l,x_u,timef,timeg,times,f,v,
     &                    xkkt,kl,ll,ictrl,lctrl,rctrl,cctrl,warn,nfe,
     &                    nge,feasible,iuser,luser,cuser,ruser,nnz,nnzh,
     &                    eqn,lin,nw,mxloop,ierr,time0,iseed,shift,
     &                    lambda,iname,sloop,z)
!     
      kltt = kltt + kl
      lltt = lltt + ll

        call matid_end(it-1,x,n)
!  return if there are errors
        if(ierr.lt.0) return

!       if (dabs(fopt-f).ge.tol_bayes) then                    ! absolute
        if (dabs(fopt-f)/dabs(fopt+1.d0).ge.tol_bayes) then    ! relative
          if (f.gt.fopt.or..not.feasible) then
            write(*,1000) it,ir,p,f,v,fopt
          else
            fopt = f
            ir   = 1
            xkktopt = xkkt
            write(*,1000) it,ir,p,f,v,fopt
          endif
        else
          if (f.lt.fopt.and.feasible) fopt=f
          ir = ir+1
          do i=1,n
            xbest(i)=x(i)
          enddo
          xkktopt = xkkt
          p  = 1.0d0-conv(it,ir)
          write(*,1000) it,ir,p,f,v,fopt

!  successful search
          if (p.ge.ptarget) then
            write(*,2000) p
            write(*,*) ' x* = ',(xbest(i),i=1,n5)
            call mesh_case(it-1)
            write(*,*) ' '
            write(*,*) ' Nfe = ',nfe,' Nge = ',nge
            write(*,*) ' '
            if (fapriori.lt.big) then
              write(*,*) ' fapriori - fopt = ',fapriori-fopt,
     &                   ' fopt = ',fopt
            else
              write(*,*) ' fopt = ',fopt
            endif
            write(*,*) ' '
            write(*,7000)
            write(*,*) ' '
!
            timei = SAOi_seconds()
            open  (16,file='test.out',status='old',position='append')
            write (16,9617) cname1,n,ni,ne,subsolver,fopt,v,
     &                      xkktopt,kltt,lltt,timei-time0,p,ptarget
            close (16)!
            return
          endif
        endif
      enddo

!  unsuccessful search ...
      write(*,3000)
      write(*,*) ' The best candidate for x* = ',
     &            (xbest(i),i=1,n5)
      call mesh_case(it-1)
      write(*,*) ' '
      write(*,*) ' Nfe = ',nfe,' Nge = ',nge
      write(*,*) ' '
      if (fapriori.lt.big) then
        write(*,*) ' fapriori - fopt = ',fapriori-fopt,
     &             ' fopt = ',fopt
      else
        write(*,*) ' fopt = ',fopt
      endif
      write(*,*) ' '
      write(*,7000)
      write(*,*) ' '
!
      if (warn) write (6,10)
!
      timei = SAOi_seconds()
      open  (16,file='test.out',status='old',position='append')
      write (16,9618) cname1,n,ni,ne,subsolver,f,v,xkktopt,
     &                nfe,nge,timei-time0,p,ptarget
      close (16)
!
      return
!
 10   format (/,' There were warnings and/or errors; see the file ',
     &        'Warnings.out',/)
 100  format (/,' Initiating global search # ',i6,5x,'<',69('='))
 1000 format (' trials = ',i6,'   successes = ',i3,
     &        '   Prob = ',1f7.4,'   f = ',1es11.4,
     &        '   h = ',1es11.4,'   fopt = ',1es14.7)
 2000 format (/,' Global optimum probably found; the probability is ',
     &           1f8.6,' <======= ',/)
 3000 format (/,' Global optimum possibly not found - maximum ',
     &           ' number of global iterations exceeded',/)
 6000 format (/,36x,'Global optimization using the SAOi algorithm',
     &        //,111('-'))
 7000 format (' ===> See History.out and Variables.out for ',
     &         'additional information ')
 9616 format (a24,1x,4i10,3es15.7,2i8,f13.2)
 9617 format (a24,1x,4i10,3es15.7,2i8,f13.2,4x,1f8.6,' > ',1f8.6)
 9618 format (a24,1x,4i10,3es15.7,2i8,f13.2,4x,1f8.6,' < ',1f8.6,' ! ')
 9626 format (4i8,3es15.7,2i8,1f10.2)
      end subroutine SAOi_split
!-----------------------------------------------------------------------
      subroutine SAOi_splita (n,m,ni,ne,x,x_l,x_u,timef,timeg,times,f,v,
     &                        xkkt,kl,ll,ictrl,lctrl,rctrl,cctrl,warn,
     &                        nfe,nge,feasible,iuser,luser,cuser,ruser,
     &                        nnz,nnzh,eqn,lin,nw,mxloop,ierr,time0,
     &                        iseed,shift,lambda,iname,sloop,z)
      implicit         none
      include          'ctrl.h'
      logical          eqn(*), lin(*)
      integer          n,m,ni,ne,nfe,nge,kloop,llooptot,kl,ll,nnz,nnzh
      integer          outerloop,nw,mxloop,ierr,iseed,iname,sloop
      logical          warn,feasible
      double precision timef,timeg,times,f,v,xkkt,shift(*)
      double precision x(n),x_l(n),x_u(n),time0,lambda(m),z(n)
      include          'ctrl_get.inc'

!  sparse implementation
      if (structure.eq.3) then
        call sao_sparse(n,m,ne,x,x_l,x_u,timef,timeg,times,f,v,
     &                  xkkt,kl,ll,ictrl,lctrl,rctrl,cctrl,warn,nfe,nge,
     &                  feasible,iuser,luser,cuser,ruser,nnz,nnzh,eqn,
     &                  lin,nw,ictrl(15),mxloop,ierr,time0,iseed,
     &                  lambda,iname,sloop,z)
     
!  dense or pseudo-sparse implementation
      elseif (structure.eq.1.or.structure.eq.2) then
        call sao_dense(n,m,ne,x,x_l,x_u,timef,timeg,times,f,v,
     &                 xkkt,kl,ll,ictrl,lctrl,rctrl,cctrl,warn,nfe,nge,
     &                 feasible,iuser,luser,cuser,ruser,shift,nnz,nnzh,
     &                 eqn,lin,nw,ictrl(15),mxloop,ierr,time0,iseed,
     &                 lambda,iname,sloop,z)

!
      else
        stop ' illegal structure in SAOi_splita '
      endif
!
      return
      end subroutine SAOi_splita


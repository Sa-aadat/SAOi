! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Approximate diagonal `Hessian' terms, to be specified by the
! SAOi: user if approx_f.eq.100 or approx_c.eq.100
! SAOi:


      subroutine diaHessUser (n,ni,ne,acurv,bcurv,ccurv,f,c,h,
     &                        x,x_h,x_h2,gf,gc,gh,gf_h,gc_h,gh_h,
     &                        acurvn,bcurvn,ccurvn,ictrl,lctrl,
     &                        rctrl,cctrl,outerloop,iuser,luser,
     &                        cuser,ruser,shift)
      implicit         none
      include          'ctrl.h'
      integer          i,j,k,l,n,ni,ne,outerloop
      double precision acurv,bcurv(ni) ,ccurv(ne),f,c(ni),h(ne)
      double precision acurvn(n),bcurvn(ni,n),ccurvn(ne,n),gf_h(n)
      double precision x(n),x_h(n),x_h2(n),gf(n),gc(ni,n),gh(ne,n)
      double precision gc_h(ni,n),gh_h(ne,n),shift(*)
      include          'ctrl_get.inc'
!=================================================================================================!
!                                                                                                 !
!  This soubroutine allows for specification of exact DIAGONAL                                    !
!    Hessian information, i.e. it is assumed that the problem is                                  !
!    separable. For non-separable problems, see dgalQP for examples.                              !
!                                                                                                 !
!    Different approximations may be used for different iterations,                               !
!    which may be handy when the dual multipliers are known/unknown.                              !
!                                                                                                 !
!    NOTE: all of this is efficiently available in sparse format!                                 !
!                                                                                                 !
!=================================================================================================!
!
!
!
!  n          : number of design variables           : int  scalar
!  ni         : number of inequality constraints     : int  scalar
!  ne         : not used in this routine             : int  scalar
!  x          : design variables                     : dble vector [n]
!  x_h        : design variables at previous iterate : dble vector [n]
!  x_h2       : design variables two iterates back   : dble vector [n]
!  f          : the objective function               : dble scalar
!  c          : the constraint functions             : dble vector [ni]
!  h          : not used in this routine             : dble vector [ne]
!  gf         : gradients of f                       : dble vector [n]
!  gc         : gradients of c                       : dble matrix [ni x n]
!  gh         : not used in this routine             : dble matrix [ne x n]
!  gf_h       : gradients of f at previous iterate   : dble vector [n]
!  gc_h       : gradients of c at previous iterate   : dble matrix [ni x n]
!  gh_h       : not used in this routine             : dble matrix [ne x n]
!  acurv      : not used in this routine             : dble scalar
!  bcurv      : not used in this routine             : dble vector [ni]
!  ccurv      : not used in this routine             : dble vector [ne]
!  acurvn     : the estimated curvatures for f       : dble vector [n]
!  bcurvn     : the estimated curvatures for c       : dble matrix [ni x n]
!  ccurvn     : not used in this routine             : dble matrix [ne x n]
!
!

!
      if (outerloop.eq.0) then
!
        if (approx_f.eq.100) then
          ! initialize the user curvatures for f here when outerloop = 0
        endif
!
        if (approx_c.eq.100) then
          ! initialize the user curvatures for c here when outerloop = 0
        endif
!
      else  ! (outerloop > 0 )
!
        if (approx_f.eq.100) then
          ! construct the user curvatures for f here when outerloop > 0
        endif
!
        if (approx_c.eq.100) then
          ! construct the user curvatures for c here when outerloop > 0
        endif
!
      endif  ! (outerloop conditions)
!
      stop ' Nothing has been coded in diaHessUser' ! remove this when done above
!
      return
      end subroutine diaHessUser
!----------------------------------------------------------------------c


!  optionally determine first and second order derivatives here - for testH2 and debugging only... ! 
        if (GetFinDiffGradHessDense.and.cname1(1:6).eq.'testH2') then
        
        
          write(*,*) 'GetFinDiffGradHessDense',cname1(1:6)
          
        
          call FinDiffGradHessDense1 (n,x,ni,ni-ne,ne,f,c,ictrl,lctrl,
     &                               rctrl,cctrl,do_diff,iuser,luser,
     &                               cuser,ruser,eqn,lin,sloop,acurv,
     &                               nnzh,xlam,Hne,Hval,Hrow,Hcol,Hptr) 
        endif

!  form the exact sparse COO Hessian if A is diagonal  - for testH3 and debugging only... !
        if (cname1(1:6).eq.'testH3') then
        
        
          write(*,*) 'testH3'
          
        
          call FormHessSparseD (n,ni,x,xlam,Hne,Hval,Hrow,Hcol,
     &                          ictrl,lctrl,rctrl,cctrl,iuser,luser,
     &                          cuser,ruser,nnzh,eqn,lin,acurv) 
        endif

!  optionally determine first and second order derivatives here - for testH4 and debugging only... ! 
        if (GetFinDiffGradHessDense.and.cname1(1:6).eq.'testH4') then
      
      
          write(*,*) 'GetFinDiffGradHessDense',cname1(1:6)
      
      
          call FinDiffGradHessDense (n,x,ni,ni-ne,ne,f,c,ictrl,lctrl,
     &                             rctrl,cctrl,do_diff,iuser,luser,
     &                             cuser,ruser,eqn,lin,sloop) 
     
!  call reconstruct sparse/dense here...
          call FinDiffGradHessDense0 (n,x,ni,ni-ne,ne,f,c,ictrl,lctrl,
     &                               rctrl,cctrl,do_diff,iuser,luser,
     &                               cuser,ruser,eqn,lin,sloop,acurv,
     &                               nnzh,xlam,Hne,Hval,Hrow,Hcol,Hptr) 
     
        endif
        
!  form the exact sparse COO Hessian if A is not diagonal - for testH5
        if (cname1(1:6).eq.'testH5') then
        
        
          write(*,*) 'FormHessSparseN','testH5' 
          
        
          call FormHessSparseN (n,ni,x,xlam,Hne,Hval,Hrow,Hcol,
     &                          ictrl,lctrl,rctrl,cctrl,iuser,luser,
     &                          cuser,ruser,nnzh,eqn,lin,acurv) 
        endif


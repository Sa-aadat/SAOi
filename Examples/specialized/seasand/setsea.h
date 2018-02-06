      mult=2          !mesh multiplier
      ground=1        !groundstructure       (1: 2bar, 
                      !                       2: beam, 
                      !                       3: lbrac)
      approx=1        !approximation         (0: noncon. ana. 
                      !                       1: abs. noncon. ana. 
                      !                       2: noncon. sph.
                      !                       3: abs. noncon. sph)
      prob=4         !problem               (1: min. weight st. local stress, 
                      !                       2: min. comp. st. vol (0.5)
                      !                       3: min. comp. st. vol (0.5) and slope  const.
                      !                       4: min. weight st. local stress and slope)
      filt0=25.  ! 1d8       !slope const. paramter ( :'mu')
      solver=1        !solver                (1: cplex, 
                      !                       2: gurobi, not implemented)
      trust=1         !trust-region          ( :yes/no)
      consv=0         !conservatism          ( :yes/no)
      dl0=1.0d-1      !init. move-limit
      dlm=1.0d0       !max. move-limit
      dlc=1d-1        !convergence tolerance
      conv=0          !convergence norm:     (0: Euclid,
                      !                       1: Abs. max)
      exprm=1         !number of experiments in series
      conti=0         !slope continuation     ( :yes/no)
      rs=0            !random start
      warm=0          !warm start
      x0=1.0d0        !init. material

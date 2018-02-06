/*			#####DCPLEX.c#####
=============================================================================

@author: C.G.P. Cleven
@license: Creative Commons 4.0
@date: 24 November 2011
@version: 0.7.11


This file contains C routines which are used by the SAOi algorithm. These 
routines call functions of the CPLEX callable libaries which is an 
optimization toolbox. The following routines are declared:

1. print_error
2. cplex_init_
3. cplex_stop_
4. cplex_solve_sqp_
future: 5. cplex_solve_qp_

Included in this file are a CPLEX header file C{cplex.h} which contains the
CPLEX callable library.

@var env Pointer to the CPLEX optimization environment
@var lp Pointer to optimization problem specified in the CPLEX
        environment.

=============================================================================*/
#include </home/sparker/opt/ibm/ILOG/CPLEX_Studio1271/cplex/include/ilcplex/cplex.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>

CPXENVptr env = NULL;
CPXLPptr  lp  = NULL;

/**
 * This routine writes error messages to the screen when a CPLEX call fails 
 *
 * @param status[in] Status of a CPLEX routine execution (0=ok, otherwise 
 *                   error code).
 * @param text[in] Text message displayed on the screen before the actual
 *                 error message from CPLEX.                 
 */
void print_error(int status, const char *text)
{
    char errmsg[1024]; ///> Variable used to store error message of CPLEX

    CPXgeterrorstring (env, status, errmsg);
    fprintf (stderr, "%s\n%s", text, errmsg);
}

/**
 * This routine starts the CPLEX environment and creates an empty problem in
 * the CPLEX environment used to store the subproblems.
 *
 * @param status[out] Status of execution (0=ok, otherwise error code).
 */
void cplex_init_(int *status)
{
    /* Start the CPLEX environment */
    if (env == NULL) {
        env = CPXopenCPLEX (status);
    }
    
    if (env == NULL) {
        print_error (*status, "Could not open CPLEX environment.");
        return;
    }

    /* Create an empty subproblem in the CPLEX environment which is
       named 'SAOi_Cplex'. */
    if (lp == NULL) {
        lp = CPXcreateprob (env, status, "SAOi_Cplex");
    }
    if (lp == NULL) {
    	print_error(*status,"CPLEX failed to create problem.");
        return;
    }

    return;
}

/**
 * This routine first closes the problem which is pressent in the CPLEX 
 * environment and then closes the CPLEX environment.
 *
 * @param status[out] Status of execution (0=ok, otherwise error code)
 */
void cplex_stop_(int *status)
{
    /* Free space used for the specified problem in the CPLEX environment.*/
    if (lp != NULL) {
        *status = CPXfreeprob(env, &lp);
        if (*status) {
        	print_error(*status, "CPLEX cpxfreeproblem failed.");
            return;
        }
    }

    /* Close the CPLEX environment.*/
    if (env != NULL) {
        *status = CPXcloseCPLEX(&env);

        if (*status) {
        	print_error(*status, "Could not close CPLEX environment.");
            return;
        }
    }

    return;
}
/**
 * This routine puts the optimization problem specified by the input
 * data into the problem C{lp} which is specified in the CPLEX environment.
 * The constructed problem is of the following form:
 * 
 *     xopt =  min gf*x + 0.5x g2Lk x
 *            s.t. gc*x <= -c
 *                 lb <= x <= ub
 * 
 * In this subroutine, first the linear part of the above problem and
 * the upper and lower bound is put into the CPLEX environment. When
 * L{loop} is equal to one (function is called for the first time) this is
 * done by using the function C{CPXcopylp}. In succeeding calls, other CPLEX
 * function are used which only change the necessary parts of the previous
 * problem.
 * 
 * After the Linear part is put into the CPLEX environment the quadratic part
 * is transfered to CPLEX. The quadratic part C{g2Lk} has only non-zero
 * elements on the diagonal which are stored in C{qpdata}.
 * 
 * As a third step the optimization problem is solved by CPLEX and the results
 * are obtained from the CPLEX environment. The results of the optimization
 * are the solution status C{solstat}, the optimal solution C{xopt} and the
 * Lagrange multipliers for the optimal solution C{lambda}.
 * 
 * @param n[in] Number of optimization variables (columns)
 * @param ni[in] Number of inequality constraints (rows)
 * @param c[in] Inequality constraint values in the approximated point
 * @param gf[in] Gradient of the objective function in the approximated point
 * @param gc[in] Gradients of the constraint functions in the approximated 
 *               point. Is used as a Vector where first columns and then rows
 *               are stored as is done in Fortran (C arrays are stored the other
 *               way around)
 * @param lb[in] Lower bounds on the optimization variables
 * @param ub[in] Upper bounds on the optimization variables
 * @param qpdata[in] Diagonal terms of the Q matrix
 * @param loop[in] Current number of the sum of Outer and inner loops of the
 *                 SAOi algorithm. This variable start at the value 1.
 * @param method[in] Method used to solve the QP problem in CPLEX.
 *                        0 = CPLEX chooses method
 *                        1 = Primal simplex
 *                        2 = dual simplex
 *                        3 = network optimizer
 *                        4 = Barrier method
 * @param feasopt[in] Method used to handle infeasibility of the optimization
 *                    problem. The following options are available:
 *                     0 = Minimize sum of all relaxations in first phase
 *                     1 = Same as 1 but in second phase optimize relaxations
 *                     2 = Minimize number of bound and constraint requiring
 *                         relaxation only in first phase. 
 *                     3 = Same as 2 but in second phase optimize relaxations
 *                     4 = Minimize sum of squares of required relaxations
 *                     5 = Same as 4 but in second phase optimize relaxations
 * @param xopt[out] The found optimum by the CPLEX QP solver. step is not made
 *                  from xprimal but from the origin.
 * @param lambda[out] Lagrange multipliers in the optimum.
 * @param status[out] Status of execution (0=ok, otherwise error code)
 */

void cplex_solve_sqp_(int *status, const int *n,const int *ni, const double *c,
                      const double *gf, const double *gc, const double *lb,
                      const double *ub, const double *qpdata,
                      double *xopt, double *lambda, const int *loop,
                      const int *method, const int *feasopt, const int *eqn)
{
    int i, j;
    int nvar = *n;
    int nc   = *ni;
    int steps; ///> maximum number of iterations.
    int feasible = 1; ///> subproblem is feasible (1 = true, 0 = false).
    char buffer_str[510];
    CPXCHARptr buffer_p;
    double crelax[*ni+10]; ///> Vector used to store relaxation variables
    
    char lu[2 * *n]; ///> List specifying if a value in L{bd} is a upper or 
                     ///> lower bound.
    double bd[2 * *n]; ///> List containing bounds which are put into the CPLEX 
                       ///> environment. Type of bounds is stored in L{lu}.
                       
    int vecind2n[2 * *n]; ///> List of length 2*n incrementing from 0 till 
                          ///> 2*n-1.
    int vecindn[*n]; ///> List of length nc incrementing from 0 till n-1
    int vecindnc[*ni]; ///> List of length nc incrementing from 0 till nc-1
    double weight[*ni]; ///> Vector of length nc containing weight values for
                        ///> relaxation.
    
    char sense[*ni]; ///> Vector specifying type of constraints
    double rscon[*ni]; ///> right hand side of the constraints. Value is set 
                       ///> to C{-1 * c}.
                       
    int optstat = 0; ///> store status of parameter setting of options. Used
                     ///> for error handling after all parameters are set
    
    int matbeg[*n]; ///> list used for storing of constraint array in CPLEX. 
                    ///> L{matbeg} is a list of length L{nvar}(number of 
                    ///> columns in L{gc}) which specifies an index related 
                    ///> to L{gc} at which point a certain column starts.
    int matcnt[*n]; ///> Number of non-zero elements in a column of L{gc}. In
                    ///> this case C{gc} contains no non-zero elements therefor
                    ///> each column is of length C{nc} (number of rows in 
                    ///> L{gc}).
    int matind[*n * *ni]; ///> list which specifies for each value in C{gc} in 
                          ///> which row (to which constraint) a value belongs.
    int colind[*n * *ni];
    int solstat; ///> Status of the solution that is obtained by optimizing the
                 ///> problem.
    int i_r;
    int i_c;
    /* construct right hand side of the constraints which is equal to
       the constraint values L{c} multiplied by -1 */
    for (i = 0; i < nc; i++) {
        rscon[i] = -1 * c[i];
    }

//  	printarrayd(gc, nvar*nc);


    /* generate index vector of lenght L{nvar}. */
    for (i = 0; i < nvar; i++) {
        vecindn[i] = i;
    }
    for (i = 0; i < nc; i++) {
      	vecindnc[i] = i;
    }

    i_r = 0;
    i_c = 0;
    for (j = 0; j< nvar*nc; j++){
   	 matind[j] = i_r;
	 colind[j] = i_c;
   	 i_r++;
 	 if (i_r == nc)
	 {
	   i_r = 0;
	   i_c++;
	 }
     }

//    printarrayd(gf, nvar*nc);
//    printarrayd(qpdata, nvar);

    /* If this is the first call to this routine, construct the whole linear
       part of the optimization problem.*/
    if (*loop == 1) {
        /* construct the vector C{sense} of length nc which specifies the kind
           of constraint:

            L = constraint <=
            E = constraint =
            G = constraint >=
            R = ranged constraint
        */
        for (i=0; i < nc; i++) {
		if(eqn[i])
		{
            		sense[i] = 'E';
		}
		else
		{
			sense[i] = 'L';
		}
	}

        /* construct vectors for the storage of the constraint array */ 
	for(j=0; j<nvar; j++){
	matbeg[j] = j*nc;
	matcnt[j] = nc;
	}

	//printarrayd(gc, nvar*nc);
	
        /* Copy the LP part of the problem data into lp */
        *status = CPXcopylp(env, lp, nvar, nc, CPX_MIN, gf, rscon, sense,
                            matbeg, matcnt, matind, gc, lb, ub,
                            NULL);
        if (*status) {
        	print_error(*status, "CPLEX failed to copy LP problem data");
            return;
        }
    }
    /* if this is not the first call to this routine, only update the gradient
       of the objective function L{fg}, the right hand side of the constraints
       L{rscon}, the constraint array L{gc} and the lower (L{lb}) and upper 
       bounds (L{ub}). */
    else {
        /* change linear part of the objective function*/
        *status = CPXchgobj(env, lp, nvar, vecindn, gf);

        if (*status) {
        	print_error(*status, "CPLEX failed to copy LP problem data");
            return;
        }

        /* put right hands side of constraint into problem. */
        *status = CPXchgrhs(env, lp, nc, vecindnc, rscon);

        if (*status) {
        	print_error(*status, "CPLEX failed to copy LP problem data");
            return;
        }

	//*status = CPXchgcoeflist(env,lp,nnz,matind,Acol,val);
        /* change the constraint matrix of the problem. */
        //for (i = 0; i < nvar; i++) {
         //   for (j = 0; j < nc; j++) {
           //     *status = CPXchgcoef(env, lp, j, i,gc[i*nc+j]);
	*status = CPXchgcoeflist(env,lp,nvar*nc,matind,colind,gc);
                if (*status) {
                	print_error(*status, "Failed to update constrained matrix.");
                    return;
                }
           // }
       // }

        /* change the lower and upper bounds of the problem. */
        for (i = 0; i < nvar; i++) {
            lu[i] = 'L';
            bd[i] = lb[i];
            lu[nvar+i] = 'U';
            bd[nvar+i] = ub [i];
            vecind2n[i] = i;
            vecind2n[nvar+i] = i;
        }

        *status = CPXchgbds(env, lp, 2*nvar, vecind2n,lu,bd);

        if (*status) {
        	print_error(*status, "Failed to update bounds.");
            return;
        }
    }

    /* Copy the separable QP part of the problem data into C{lp}. */
    *status = CPXcopyqpsep(env, lp, qpdata);

    if (*status) {
    	print_error(*status, "CPLEX failed to copy QP problem data.");
        return;
    }

    /* Set options for the optimization of the specified problem. */
    *status = CPXsetintparam(env, 1063, *method); // set method used for 
    if (*status) {optstat = *status;}             // optimization of the QP 
                                                  // problem
    *status = CPXsetdblparam(env, 3002, 1e-3); // set convergence tolerance 
    if (*status) {optstat = *status;}          // for LP and QP
   
    *status = CPXsetdblparam(env, 1014, 1e-3); // optimality tolerance 
    if (*status) {optstat = *status;}    
    	
    *status = CPXsetdblparam(env, 1016, 1e-6); // feasibility tolerance 
    if (*status) {optstat = *status;}   
    
    steps = 10000;
    if (nc > steps) {
    	steps = nc;
    }
    *status = CPXsetintparam(env, 3012, steps); // barrier iteration limit 
    if (*status) {optstat = *status;}    
    	
    *status = CPXsetintparam(env, 1020, steps); // simplex iteration limit 
    if (*status) {optstat = *status;}  
                                               
	*status = CPXsetintparam(env, 1084, *feasopt);// Set method used to handle
    if (*status) {optstat = *status;}	          // infeasible subproblems
    	
	/* Check if all parameters are set correctly */
    if (optstat) {
    	print_error(optstat, "CPLEX failed to copy parameter value.");
		return;
	}

        	*status = CPXwriteprob(env, lp, "qp.lp", "LP");
    /* Solve the QP problem with the CPLEX QP solver */
    *status = CPXqpopt(env, lp);

    if (*status) {
    	print_error(*status, "CPLEX failed to optimize the subproblem.");
        return;
    }

    /* Get the status of the found solution and check if the optimization was
       successful. */
    solstat = CPXgetstat(env, lp);

    if (solstat == 0) {
        fprintf (stderr, "Solution status is 0 which indicates an error\n");
        *status = 1;
        return;
    }
    
    /* Check if problem was infeasible. If so, generate infeasible problem by
       relaxing the rhs of contraints. The actual qp model is updated with 
       a purposed relaxation after which the optmization problem is again 
       optmized. */
    if (solstat == CPX_STAT_INFEASIBLE) {

    	feasible = 0;
    	
    	/* construct weight vector containing only ones. */
    	for (i = 0; i < nc; i++) {
    		weight[i] = 1e-6;
    	}

    	while (feasible < 5) {
        	/* Calculate feasible solution */
        	
        	*status = CPXfeasopt(env, lp, weight, NULL, NULL, NULL);
        	
        	if (*status != 0) {
        		print_error(*status, "CPXfeasopt failed\n");
            	*status = 1;
            	return;
        	}
        	
        	/* Determine solution status. */
        	solstat = CPXgetstat(env, lp);

//		fprintf(stderr, "solstat %d \n", solstat);
        	
    
//        	if ((solstat < 14) || ((solstat > 20) && (solstat != 122))) {
//            	fprintf (stderr, "CPXfeasopt gave no solution\n");
//            	buffer_p = CPXgetstatstring(env, solstat, buffer_str);
//            	printf("solution status of feasopt is: %d\n", solstat);
//            	printf("%s\n", buffer_str);
//            	*status = 1;
//            	return;
//            }
            
            /* Get relaxation variables for updating of the rhs of the sub
               problem.*/ 
            *status = CPXgetrowinfeas(env, lp, NULL, crelax, 0, CPXgetnumrows(env,lp)-1 );
            
            if (*status) {
        		print_error(*status, "Could not obtain relaxed variables.");
            	return;
        	}
        	
            *status = CPXgetx(env, lp, xopt, 0, nvar-1);
    
        	if (*status) {
        		print_error(*status, "Failed to obtain infeasible xopt.");
            	return;
        	}

        	/* update right hand side of the linear constraints. */
        	for (i = 0; i < nc; i++) {
        		rscon[i] = rscon[i] + crelax[i];	
            }
    
    		*status = CPXchgrhs(env, lp, nc, vecindnc, rscon);
    
            if (*status) {
            	print_error(*status, "CPLEX failed to update rhs with relaxation.\n");
                return;
            }
            
            /* set starting point for optimization of relaxed problem. */
            
            
            /* optimize relaxed QP problem. */
            *status = CPXqpopt(env, lp);
    
        	if (*status) {
        		print_error(*status, "CPLEX failed to optimize relaxed problem.");
            	return;
        	}

        	
        	/* Get the status of the found solution and check if the optimization was
               successful. */
        	solstat = CPXgetstat(env, lp);
        	
    
        	if (solstat == 0) {
            	fprintf (stderr, "Solution status is 0 which indicates an error\n");
            	*status = 1;
            	return;
        	}
        	
        	/* write problem */
        	*status = CPXwriteprob(env, lp, "fchecke.lp", "LP");
    
        	if (solstat == 3) {
        		printf("CPLEX cannot construct feasible problem\n");
			*status = 1;
			return;
        	}
        	if (solstat == 1) {
        		// printf("foundsolution is feasible, continue");
        		feasible = 10;
        	}
        }
    }

    /* Get the optimization results from the CPLEX environment. First, get the
       minimizer.*/
    *status = CPXgetx(env, lp, xopt, 0, nvar-1);

    if (*status) {
    	print_error(*status, "Failed to obtain xopt from CPLEX.");
        return;
    }

    /* update lagrange multipliers if feasible subproblem is solved. */
    *status = CPXgetpi(env, lp, lambda, 0, nc-1);

    if (*status) {
    	print_error(*status, "Failed to obtain Lagrange multipliers form CPLEX");
        return;
    }

    return;
}

/**
 *
 *
 *
 *
 *
 *
 *
 */
void cplex_solve_qp_(int *status) {
	
	
	
}

void printarrayi(const int *array, int size)
{
	int i;
	for (i = 0; i < size; i++)
	{
		fprintf(stderr, "%d, ", array[i]);
	}
	fprintf(stderr, "\n");
}

void printarrayd(const double *array, int size)
{
	        int i;
		for (i = 0; i < size; i++)	
		{
	                fprintf(stderr, "%f, ", array[i]);	
		}
		fprintf(stderr, "\n");
}

/**
 * 
 * Solves separable qp problems with sparse arithmetic
 * 
 */
void cplex_solve_sqps_sep_(int *status, const int *n,const int *ni,
                           const double *c,const double *gf, const double *gc,
                      const double *lb,const double *ub, const double *qpdata,
                      double *xopt, double *lambda, const int *loop,
                      const int *method, const int *feasopt, const int *nnz_c,
                      const int *Acol_c, const int *Aptr_c,
                       int *storind, double *delta_t, const int *eqn)
{
    int i, j;
    int nnz = *nnz_c;
    int nvar = *n;
    int nc   = *ni;
    int steps; ///> maximum number of iterations.
    int feasible = 1; ///> subproblem is feasible (1 = true, 0 = false).
    char buffer_str[510];
    CPXCHARptr buffer_p;
    double crelax[*ni+10]; ///> Vector used to store relaxation variables
    
    char lu[2 * *n]; ///> List specifying if a value in L{bd} is a upper or 
                     ///> lower bound.
    double bd[2 * *n]; ///> List containing bounds which are put into the CPLEX 
                       ///> environment. Type of bounds is stored in L{lu}.
                       
    int vecind2n[2 * *n]; ///> List of length 2*n incrementing from 0 till 
                          ///> 2*n-1.
    int vecindn[*n]; ///> List of length nc incrementing from 0 till n-1
    int vecindnc[*ni]; ///> List of length nc incrementing from 0 till nc-1
    double weight[*ni]; ///> Vector of length nc containing weight values for
                        ///> relaxation.
    
    char sense[*ni]; ///> Vector specifying type of constraints
    double rscon[*ni]; ///> right hand side of the constraints. Value is set 
                       ///> to C{-1 * c}.
                       
    int optstat = 0; ///> store status of parameter setting of options. Used
                     ///> for error handling after all parameters are set
                     
                     //Maybe this should also run only the first time
    int Aptr[nc+ 1]; 	//memcpy(Aptr, Aptr_c,sizeof(Aptr_c));
    int Acol[nnz]; 	//memcpy(Acol, Acol_c,sizeof(Acol_c));

    double val[nnz];

    int matbeg[nvar]; ///> list used for storing of constraint array in CPLEX. 
                    ///> L{matbeg} is a list of length L{nvar}(number of 
                    ///> columns in L{gc}) which specifies an index related 
                    ///> to L{gc} at which point a certain column starts.
    int matcnt[nvar]; ///> Number of non-zero elements in a column of L{gc}. In
                    ///> this case C{gc} contains no non-zero elements therefor
                    ///> each column is of length C{nc} (number of rows in 
                    ///> L{gc}).
    int matind[nnz]; ///> list which specifies for each value in C{gc} in 
                          ///> which row (to which constraint) a value belongs.
    int solstat; ///> Status of the solution that is obtained by optimizing the
                 ///> problem.
    int Acol_tempvec[nnz];
    int matind_tempvec[nnz];
    double val_temp;
    double t_1, t_b;
    double t_2, t_e;
    int cnt = 0;
    int col;
    int prog[nvar];
    int Cptr[nvar+1];


   for (i = 0; i < nnz; i++)
    {
	Acol_tempvec[i] = Acol_c[i]-1;	
    } 
   for (i = 0; i < nc+ 1; i ++)
    {
	Aptr[i] = Aptr_c[i]-1;
    }
     /* This section determines the row index of each entry in the constraint
      * jacobian.
      * Function call: fmatind(n,ptr,matind);
      */
      t_1 = (double)clock();
      j = 0;
	  for (i = 0; i < nc+1; i++)
	  {
	    matind_tempvec[j] = 0;  
	    for (j; j < Aptr[i+1]; j++)
	    {
	      matind_tempvec[j] = i;   
	    }
	  }

    t_2 = (double)clock();
    *delta_t = *delta_t + (t_2 - t_1)/CLOCKS_PER_SEC;

    /* construct right hand side of the constraints which is equal to
       the constraint values L{c} multiplied by -1 --> "constant term in taylor series expansion"*/
    for (i = 0; i < nc; i++) {
       rscon[i] = -1 * c[i];
    }
    
    /* generate index vector of lenght L{nvar}. */
    for (i = 0; i < nvar; i++) {
        vecindn[i] = i;
    }
    for (i = 0; i < nc; i++) {
      	vecindnc[i] = i;
    }
    /* If this is the first call to this routine, construct the whole linear
       part of the optimization problem.*/

    
    if (*loop == 1) {
	
	/*Copy values for use*/
        /* construct the vector C{sense} of length nc which specifies the kind
           of constraint:

            L = constraint <=
            E = constraint =
            G = constraint >=
            R = ranged constraint
        */

        for (i=0; i < nc; i++){
            if(eqn[i])
	    {
		sense[i] = 'E';
	    }
	    else
	    {
            	sense[i] = 'L';
	    }
        }
        t_1 = (double)clock();
        /* construct vectors for the storage of the constraint array */

	/*Insert first time stuff*/
	/*Construct and build matcnt---> number of elements in each column*/	
	 for (j = 0; j < nvar; j++)
	 {
		 matcnt[j] = 0;

	 }
	 for (j = 0; j < nnz; j++)
	 {
		matcnt[Acol_tempvec[j]] = matcnt[Acol_tempvec[j]] + 1;
		
	 }
	 
	 /*Construct matbeg-->index in constraint array of first entry in each collumn (stored columnwise)
	  * prog--> counter
	  * Cptr--> see Aptr, for columns*/
	 Cptr[0] = 0;
	 for (j = 0; j < nvar; j++)
	 {
		 matbeg[j] = 0;
		 prog[j] = 0;
		 Cptr[j+1] = matcnt[j] + Cptr[j];
		 matbeg[j] = Cptr[j];
	 }

	 /*Sort gc rowise --> columnwise. Construct index vector to remember mapping and sort matind.*/
	 for(i = 0; i < nnz; i++)
	 {
		col = Acol_tempvec[i];
		val[Cptr[col] + prog[col] ] = gc[i];
		storind[Cptr[col]+prog[col]] = i;

		matind[Cptr[col] + prog[col]] = matind_tempvec[i];
		prog[col] = prog[col] + 1;

	 }

	t_e = (double)clock();
        t_2 = (double)clock();
	*delta_t = *delta_t + (t_2 - t_1)/CLOCKS_PER_SEC;

	
        /* Copy the LP part of the problem data into lp */
        *status = CPXcopylp(env, lp, nvar, nc, CPX_MIN, gf, rscon, sense,
                            matbeg, matcnt, matind, val, lb, ub,
                            NULL);
        if (*status) {
        	print_error(*status, "CPLEX failed to copy LP problem data");
            return;
        }

        /* Copy the LP part of the problem data into lp */
	    

    }
    /* if this is not the first call to this routine, only update the gradient
       of the objective function L{fg}, the right hand side of the constraints
       L{rscon}, the constraint array L{gc} and the lower (L{lb}) and upper 
       bounds (L{ub}). */
    else {
	    
        /* change linear part of the objective function*/
        *status = CPXchgobj(env, lp, nvar, vecindn, gf);

        if (*status) {
        	print_error(*status, "CPLEX failed to copy LP problem data");
            return;
        }

        /* put right hands side of constraint into problem. */
        *status = CPXchgrhs(env, lp, nc, vecindnc, rscon);
        if (*status) {
        	print_error(*status, "CPLEX failed to copy LP problem data");
            return;
        }
               
        /*UPDATE*/
        /* change the constraint matrix of the problem. Check how CPXchgcoef works....,.
	 env --> pointer to CPLEX environment
	 lp--> pointer to problem object
	 j--> row where coef is located
	 i--> columnn where coef is located
	 gc--> new value*/
        t_1 = (double)clock();	
	for(i = 0; i < nnz; i++)
        {
            val[i] = gc[storind[i]];
	    Acol[i] = Acol_tempvec[storind[i]];
            matind[i] = matind_tempvec[storind[i]];
        }
	t_2 = (double)clock();
	*delta_t = *delta_t + (t_2 - t_1)/CLOCKS_PER_SEC;

                   *status = CPXchgcoeflist(env,lp,nnz,matind,Acol,val);
                   if (*status) {
                	print_error(*status, "Failed to update constrained matrix.");
                   return;
               }
	
        /* change the lower and upper bounds of the problem. */
        for (i = 0; i < nvar; i++) {
            lu[i] = 'L';
            bd[i] = lb[i];
            lu[nvar+i] = 'U';
            bd[nvar+i] = ub [i];
            vecind2n[i] = i;
            vecind2n[nvar+i] = i;
        }

        *status = CPXchgbds(env, lp, 2*nvar, vecind2n,lu,bd);

        if (*status) {
        	print_error(*status, "Failed to update bounds.");
            return;
        }


    }

    /* Copy the separable QP part of the problem data into C{lp}. */
    *status = CPXcopyqpsep(env, lp, qpdata);
    
    if (*status) {
    	print_error(*status, "CPLEX failed to copy QP problem data.");
        return;
    }
    /* Set options for the optimization of the specified problem. */
    *status = CPXsetintparam(env, 1063, *method); // set method used for 
    if (*status) {optstat = *status;}             // optimization of the QP 
                                                  // problem
    *status = CPXsetdblparam(env, 3002, 1e-6); // set convergence tolerance 
    if (*status) {optstat = *status;}          // for LP and QP
   
    *status = CPXsetdblparam(env, 1014, 1e-6); // optimality tolerance 
    if (*status) {optstat = *status;}    
    	
    *status = CPXsetdblparam(env, 1016, 1e-6); // feasibility tolerance 
    if (*status) {optstat = *status;}   
    
    steps = 10000;
    if (nc > steps) {
    	steps = nc;
    }
    *status = CPXsetintparam(env, 3012, steps); // barrier iteration limit 
    if (*status) {optstat = *status;}    
    	
    *status = CPXsetintparam(env, 1020, steps); // simplex iteration limit 
    if (*status) {optstat = *status;}  

                                               
	*status = CPXsetintparam(env, 1084, *feasopt);// Set method used to handle
    if (*status) {optstat = *status;}	          // infeasible subproblems
    	
	/* Check if all parameters are set correctly */
    if (optstat) {
    	print_error(optstat, "CPLEX failed to copy parameter value.");
		return;
	}
    
    /* Solve the QP problem with the CPLEX QP solver */
    *status = CPXqpopt(env, lp);
    if (*status) {
    	print_error(*status, "CPLEX failed to optimize the subproblem.");
        return;
    }

    /* Get the status of the found solution and check if the optimization was
       successful. */
    solstat = CPXgetstat(env, lp);
    if (solstat == 0) {
        fprintf (stderr, "Solution status is 0 which indicates an error\n");
        *status = 1;
        return;
    }

    
    /* Check if problem was infeasible. If so, generate infeasible problem by
       relaxing the rhs of contraints. The actual qp model is updated with 
       a purposed relaxation after which the optmization problem is again 
       optmized. */


    if (solstat == CPX_STAT_INFEASIBLE) {

    	feasible = 0;
    	/* construct weight vector containing only ones. */
    	for (i = 0; i < nc; i++) {
    		weight[i] = 1e-6;
    	}

    	while (feasible < 5) {
        	*status = CPXwriteprob(env, lp, "prerelax.lp", "LP");
        	/* Calculate feasible solution */
        	
        	*status = CPXfeasopt(env, lp, weight, NULL, NULL, NULL);
        	if (*status != 0) {
        		print_error(*status, "CPXfeasopt failed\n");
            	*status = 1;
            	return;
        	}
        	
        	/* Determine solution status. */
        	solstat = CPXgetstat(env, lp);
        	
    
//        	if ((solstat < 14) || ((solstat > 20) && (solstat != 122))) {
 //           	fprintf (stderr, "CPXfeasopt gave no solution\n");
 //           	buffer_p = CPXgetstatstring(env, solstat, buffer_str);
 //           	printf("solution status of feasopt is: %d\n", solstat);
 //           	printf("%s\n", buffer_str);
  //          	*status = 1;
  //          	return;
  //          }
            
            /* Get relaxation variables for updating of the rhs of the sub
               problem. Also get x for a good estimate of the starting point. */ 
            *status = CPXgetrowinfeas(env, lp, NULL, crelax, 0, CPXgetnumrows(env,lp)-1 );
            
            if (*status) {
       		print_error(*status, "Could not obtain relaxed variables.");
           	return;
        	}
        	
            *status = CPXgetx(env, lp, xopt, 0, nvar-1);
        	if (*status) {
        		print_error(*status, "Failed to obtain infeasible xopt.");
            	return;
        	}

        	/* update right hand side of the linear constraints. */
        	for (i = 0; i < nc; i++) {
        		rscon[i] = rscon[i] + crelax[i];
            }
    
    		*status = CPXchgrhs(env, lp, nc, vecindnc, rscon);
    
            if (*status) {
            	print_error(*status, "CPLEX failed to update rhs with relaxation.\n");
                return;
            }
            
            /* set starting point for optimization of relaxed problem. */
            
            
            /* optimize relaxed QP problem. */
            *status = CPXqpopt(env, lp);
    
        	if (*status) {
        		print_error(*status, "CPLEX failed to optimize relaxed problem.");
            	return;
        	}
        	
        	/* Get the status of the found solution and check if the optimization was
               successful. */
        	solstat = CPXgetstat(env, lp);
        	
    
        	if (solstat == 0) {
            	fprintf (stderr, "Solution status is 0 which indicates an error\n");
            	*status = 1;
            	return;
        	}
        	
        	/* write problem */
	
        	*status = CPXwriteprob(env, lp, "postrelax.lp", "LP");
    
        	if (solstat == 3) {
        		printf("CPLEX cannot construct feasible problem\n");
        		printf("Set cplex_method to 1,2 or 3.\n");
        		printf("Press Enter:\n");
        		getchar();
        		*status = 1;
        		return;
        	}
		
        	if (solstat == 1) {
        		// printf("foundsolution is feasible, continue");
        		feasible = 10;
        	}
        }
    }

    /* Get the optimization results from the CPLEX environment. First, get the
       minimizer.*/
    *status = CPXgetx(env, lp, xopt, 0, nvar-1);

    if (*status) {
    	print_error(*status, "Failed to obtain xopt from CPLEX.");
        return;
    }

    /* update lagrange multipliers if feasible subproblem is solved. */
    *status = CPXgetpi(env, lp, lambda, 0, nc-1);

    if (*status) {
    	print_error(*status, "Failed to obtain Lagrange multipliers form CPLEX");
        return;
    }
    //fprintf(stderr, "Time in solve_qps: %f", (t_e - t_b)/CLOCKS_PER_SEC );
    return;
}

/**
 * Sa-aadat Parker
 * Wed 11 Oct 2017 07:51:45 SAST
 * Solves qp problems with sparse arithmetic
 */
void cplex_solve_sqps_(int *status, const int *n,const int *ni,
                       const double *c, const double *gf, const double *gc,
                       const double *lb, const double *ub,
                      double *xopt, double *lambda, const int *loop,
                      const int *method, const int *feasopt, const int *nnz_c,
                      const int *Acol_c, const int *Aptr_c, const int *nnz_h,
                      const int *Acol_h, const int *Aptr_h,const double *hval,
                       int *storind, double *delta_t, const int *eqn)
{
    int i;
    int j;
    int k;
    int nnz = *nnz_c; /* nnz equals value pointed at by nnz_c*/
    int nnzh = *nnz_h;
    int nvar = *n;
    int nc   = *ni;
    int steps; ///> maximum number of iterations.
    int feasible = 1; ///> subproblem is feasible (1 = true, 0 = false).
    char buffer_str[510];
    CPXCHARptr buffer_p;
    double crelax[*ni+10]; ///> Vector used to store relaxation variables

    char lu[2 * *n]; ///> List specifying if a value in L{bd} is a upper or
                     ///> lower bound.
    double bd[2 * *n]; ///> List containing bounds which are put into the CPLEX
                       ///> environment. Type of bounds is stored in L{lu}.

    int vecind2n[2 * *n]; ///> List of length 2*n incrementing from 0 till
                          ///> 2*n-1.
    int vecindn[*n]; ///> List of length nc incrementing from 0 till n-1
    int vecindnc[*ni]; ///> List of length nc incrementing from 0 till nc-1
    double weight[*ni]; ///> Vector of length nc containing weight values for
                        ///> relaxation.
    double rscon[*ni]; ///> right hand side of the constraints. Value is set
                       ///> to C{-1 * c}.
    char sense[*ni]; ///> Vector specifying type of constraints

    int optstat = 0; ///> store status of parameter setting of options. Used
                     ///> for error handling after all parameters are set
                     //Maybe this should also run only the first time
    int Aptr[nc+1]; //memcpy(Aptr, Aptr_c,sizeof(Aptr_c));
    int Acol[nnz]; 	//memcpy(Acol, Acol_c,sizeof(Acol_c));

    double val[nnz];

    int matbeg[nvar]; ///> list used for storing of constraint array in CPLEX.
                    ///> L{matbeg} is a list of length L{nvar}(number of
                    ///> columns in L{gc}) which specifies an index related
                    ///> to L{gc} at which point a certain column starts.
    int matcnt[nvar]; ///> Number of non-zero elements in a column of L{gc}. In
                    ///> this case C{gc} contains no non-zero elements therefor
                    ///> each column is of length C{nc} (number of rows in
                    ///> L{gc}).
    int matind[nnz]; ///> list which specifies for each value in C{gc} in
                     ///> which row (to which constraint) a value belongs.
    int solstat; ///> Status of the solution that is obtained by optimizing the
                 ///> problem.
    int Acol_tempvec[nnz];
    int matind_tempvec[nnz];
    double val_temp;
    double t_1, t_b;
    double t_2, t_e;
    int cnt = 0;
    int col;
    int prog[nvar];
    int Cptr[nvar+1];

//     Variables initialised below for quadratic matrix in SQP
    int qmatAptr[nvar+1];/*int Aptr[nc+ 1]*/
    int qmatAcol[nnzh];
    int qmatArow[nnzh];
    int qmatbeg[nvar];    /*int nvar = *n;*/
    int qmatcnt[nvar]; /*int nc = *ni;*/
    int qmatind[nnzh];
//     double hval[nnzh]; //equivalent to gc
    double qmatval[nnzh]; //equivalent to gc
    int qmatAcol_tempvec[nnzh];
    int qmatind_tempvec[nnzh];

   for (i = 0; i < nnz; i++) {/*LP loop*/
    Acol_tempvec[i] = Acol_c[i]-1;
   }

   for (i = 0; i < nc+1; i++) {
    Aptr[i] = Aptr_c[i]-1;
   }

   for (i = 0; i < nnzh; i++) {/*hessian loop*/
    qmatAcol_tempvec[i] = Acol_h[i]-1;
   }

   for (i = 0; i < nvar+ 1; i++) {
    qmatAptr[i] = Aptr_h[i]-1;
   }

     /* This section determines the row index of each entry in the constraint
      * jacobian. Function call: fmatind(n,ptr,matind);
      */
    t_1 = (double)clock();

    j = 0;
	  for (i = 0; i < nc+1; i++){ /*LP loop*/
	    matind_tempvec[j] = 0;
	    for (j; j < Aptr[i+1]; j++)  {
	      matind_tempvec[j] = i;
	    }
	  }

    j = 0; /*hessian loop*/
	  for (i = 0; i < nvar+1; i++) {
	    qmatind_tempvec[j] = 0;
	    for (j; j < qmatAptr[i+1]; j++)  {
	      qmatind_tempvec[j] = i;
	    }
	  }

    t_2 = (double)clock();
    *delta_t = *delta_t + (t_2 - t_1)/CLOCKS_PER_SEC;

    /* construct right hand side of the constraints which is equal to
       the constraint values L{c} multiplied by -1
       --> "constant term in taylor series expansion"*/
    for (i = 0; i < nc; i++) {
       rscon[i] = -1 * c[i];
    }

    /* generate index vector of length L{nvar}. */
    for (i = 0; i < nvar; i++) {
        vecindn[i] = i;
    }
    for (i = 0; i < nc; i++) {
      	vecindnc[i] = i;
    }
    /* If this is the first call to this routine, construct the whole linear
       part of the optimization problem.*/

    if (*loop == 1) {

      /*Copy values for use*/
      /* construct the vector C{sense} of length nc which specifies the kind
      of constraint:

      L = constraint <=
      E = constraint =
      G = constraint >=
      R = ranged constraint */

      for (i=0; i < nc; i++){
        if(eqn[i])  {
          sense[i] = 'E';/*set up for an equality constraint */
          }
        else {
          sense[i] = 'L';/*set up for an leq inequality */
        }
      }
      t_1 = (double)clock();
        /* construct vectors for the storage of the constraint array */

      /*Insert first time stuff*/
      /*Construct and build matcnt---> number of elements in each column*/
      for (j = 0; j < nvar; j++) {
        matcnt[j] = 0;
      }
      for (j = 0; j < nnz; j++) {
        matcnt[Acol_tempvec[j]] = matcnt[Acol_tempvec[j]] + 1;
      }

      /*Construct matbeg-->index in constraint array of first entry in each column
      * (stored columnwise)
      * prog--> counter
      * Cptr--> see Aptr, for columns*/
      Cptr[0] = 0;
      for (j = 0; j < nvar; j++) {
        matbeg[j] = 0;
        prog[j] = 0;
        Cptr[j+1] = matcnt[j] + Cptr[j];
        matbeg[j] = Cptr[j];
      }

      /*Sort gc rowise --> columnwise.
      *Construct index vector to remember mapping and sort matind.*/
      for(i = 0; i < nnz; i++)	 {
        col = Acol_tempvec[i];
        val[Cptr[col] + prog[col] ] = gc[i];
//         storind[Cptr[col]+prog[col]] = i;

        matind[Cptr[col] + prog[col]] = matind_tempvec[i];
        prog[col] = prog[col] + 1;
      }

      t_e = (double)clock();
      t_2 = (double)clock();
      *delta_t = *delta_t + (t_2 - t_1)/CLOCKS_PER_SEC;

      /* Copy the LP part of the problem data into lp */
      *status = CPXcopylp(env, lp, nvar, nc, CPX_MIN, gf, rscon, sense,
                      matbeg, matcnt, matind, val, lb, ub,NULL);
      if (*status) {
        print_error(*status, "CPLEX failed to copy LP problem data");
      return;
      }
      /* Copy the LP part of the problem data into lp */

      /*Hessian construction */
      /*Construct and build qmatcnt---> number of elements in each column*/
      for (j = 0; j < nvar; j++) {
        qmatcnt[j] = 0;
      }
      for (j = 0; j < nnzh; j++)	 {
        qmatcnt[qmatAcol_tempvec[j]] = qmatcnt[qmatAcol_tempvec[j]] + 1;
      }

      /*Construct qmatbeg-->index in constraint array
       * of first entry in each column (stored columnwise)
      * prog--> counter
      * Cptr--> see Aptr, for columns*/
      Cptr[0] = 0;
      for (j = 0; j < nvar; j++) {
        qmatbeg[j] = 0;
        prog[j] = 0;
        Cptr[j+1] = qmatcnt[j] + Cptr[j];
        qmatbeg[j] = Cptr[j];
      }

      /*Sort gc row wise --> columnwise.
      *Construct index vector to remember mapping and sort qmatind.*/
      for(i = 0; i < nnzh; i++)  {
        col = qmatAcol_tempvec[i];
        qmatval[Cptr[col] + prog[col] ] = hval[i];/*0*/
        storind[Cptr[col]+prog[col]] = i;

        qmatind[Cptr[col] + prog[col]] = qmatind_tempvec[i];
        prog[col] = prog[col] + 1;
      }

      /* Copy the QP part of the problem data into C{lp}. */
      *status = CPXcopyquad(env,lp,qmatbeg,qmatcnt,qmatind,qmatval);
      /* CPXcopyqpsep(env, lp, qpdata);*/

      if (*status) {
        print_error(*status, "CPLEX failed to copy QP problem data.");
          return;
      }
    }
    /* if this is not the first call to this routine, only update the gradient
       of the objective function L{fg}, the right hand side of the constraints
       L{rscon}, the constraint array L{gc} and the lower (L{lb}) and upper
       bounds (L{ub}). */
    else {

        /* change linear part of the objective function*/
        *status = CPXchgobj(env, lp, nvar, vecindn, gf);

        if (*status) {
        	print_error(*status, "CPLEX failed to copy LP problem data");
            return;
        }

        /* put right hands side of constraint into problem. */
        *status = CPXchgrhs(env, lp, nc, vecindnc, rscon);
        if (*status) {
        	print_error(*status, "CPLEX failed to copy LP problem data");
            return;
        }

       /*UPDATE*/
       /* change the constraint matrix of the problem.
       env --> pointer to CPLEX environment
       lp--> pointer to problem object
       j--> row where coef is located
       i--> column where coef is located
       gc--> new value */
        t_1 = (double)clock();
        for(i = 0; i < nnz; i++)  {
          val[i] = gc[storind[i]];
          Acol[i] = Acol_tempvec[storind[i]];
          matind[i] = matind_tempvec[storind[i]];
        }

        t_2 = (double)clock();
        *delta_t = *delta_t + (t_2 - t_1)/CLOCKS_PER_SEC;

        *status = CPXchgcoeflist(env,lp,nnz,matind,Acol,val);
        if (*status) {
          print_error(*status, "Failed to update constraint matrix.");
        return;
        }

        /* change the lower and upper bounds of the problem. */
        for (i = 0; i < nvar; i++) {
            lu[i] = 'L';
            bd[i] = lb[i];
            lu[nvar+i] = 'U';
            bd[nvar+i] = ub [i];
            vecind2n[i] = i;
            vecind2n[nvar+i] = i;
        }

        *status = CPXchgbds(env, lp, 2*nvar, vecind2n,lu,bd);

        if (*status) {
        	print_error(*status, "Failed to update bounds.");
            return;
        }

//         Hessian row vector determination
        k=0;
        for (i = 0; i < nvar+1; i++) {
          for (j = 0; j < qmatAptr[i+1]-qmatAptr[i]; j++) {
          qmatArow[k]=i;
          k=k+1;
          }
        }

//         Put in this code for Hessian
        for (i = 0; i < nnzh; i++) {
          qmatval[i]=hval[i];
          qmatAcol[i] = qmatAcol_tempvec[i];
          *status = CPXchgqpcoef (env, lp, qmatArow[i],
                                 qmatAcol[i], qmatval[i]);
          if (*status) {
            print_error(*status, "Failed to update QP data.");
          return;
          }
        }
    }


    /* Set options for the optimization of the specified problem. */
    *status = CPXsetintparam(env, 1063, *method); // set method used for
    if (*status) {optstat = *status;}             // optimization of the QP
                                                  // problem
    *status = CPXsetdblparam(env, 3002, 1e-6); // set convergence tolerance
    if (*status) {optstat = *status;}          // for LP and QP

    *status = CPXsetdblparam(env, 1014, 1e-6); // optimality tolerance
    if (*status) {optstat = *status;}

    *status = CPXsetdblparam(env, 1016, 1e-6); // feasibility tolerance
    if (*status) {optstat = *status;}

    steps = 10000;
    if (nc > steps) {
    	steps = nc;
    }
    *status = CPXsetintparam(env, 3012, steps); // barrier iteration limit
    if (*status) {optstat = *status;}

    *status = CPXsetintparam(env, 1020, steps); // simplex iteration limit
    if (*status) {optstat = *status;}



    *status = CPXsetintparam(env, 1084, *feasopt);// Set method used to handle
    if (*status) {optstat = *status;}	          // infeasible subproblems

	/* Check if all parameters are set correctly */
    if (optstat) {
    	print_error(optstat, "CPLEX failed to copy parameter value.");
		return;
	}

    /* Solve the QP problem with the CPLEX QP solver */
    *status = CPXqpopt(env, lp);
    if (*status) {
    	print_error(*status, "CPLEX failed to optimize the subproblem, sorry");
        return;
    }

    /* Get the status of the found solution and check if the optimization was
       successful. */
    solstat = CPXgetstat(env, lp);
    if (solstat == 0) {
        fprintf (stderr, "Solution status is 0 which indicates an error\n");
        *status = 1;
        return;
    }


    /* Check if problem was infeasible. If so, generate infeasible problem by
       relaxing the rhs of contraints. The actual qp model is updated with
       a purposed relaxation after which the optmization problem is again
       optmized. */


    if (solstat == CPX_STAT_INFEASIBLE) {

    	feasible = 0;
    	/* construct weight vector containing only ones. */
    	for (i = 0; i < nc; i++) {
    		weight[i] = 1e-6;
    	}

    	while (feasible < 5) {
        	*status = CPXwriteprob(env, lp, "prerelax.lp", "LP");
        	/* Calculate feasible solution */

        	*status = CPXfeasopt(env, lp, weight, NULL, NULL, NULL);
        	if (*status != 0) {
        		print_error(*status, "CPXfeasopt failed\n");
            	*status = 1;
            	return;
        	}

        	/* Determine solution status. */
        	solstat = CPXgetstat(env, lp);


//        	if ((solstat < 14) || ((solstat > 20) && (solstat != 122))) {
 //           	fprintf (stderr, "CPXfeasopt gave no solution\n");
 //           	buffer_p = CPXgetstatstring(env, solstat, buffer_str);
 //           	printf("solution status of feasopt is: %d\n", solstat);
 //           	printf("%s\n", buffer_str);
  //          	*status = 1;
  //          	return;
  //          }

            /* Get relaxation variables for updating of the rhs of the sub
               problem. Also get x for a good estimate of the starting point. */
            *status = CPXgetrowinfeas(env, lp, NULL, crelax, 0, CPXgetnumrows(env,lp)-1 );

            if (*status) {
       		print_error(*status, "Could not obtain relaxed variables.");
           	return;
        	}

            *status = CPXgetx(env, lp, xopt, 0, nvar-1);
        	if (*status) {
        		print_error(*status, "Failed to obtain infeasible xopt.");
            	return;
        	}

        	/* update right hand side of the linear constraints. */
        	for (i = 0; i < nc; i++) {
        		rscon[i] = rscon[i] + crelax[i];
            }

    		*status = CPXchgrhs(env, lp, nc, vecindnc, rscon);

            if (*status) {
            	print_error(*status, "CPLEX failed to update rhs with relaxation.\n");
                return;
            }

            /* set starting point for optimization of relaxed problem. */


            /* optimize relaxed QP problem. */
            *status = CPXqpopt(env, lp);

        	if (*status) {
        		print_error(*status, "CPLEX failed to optimize relaxed problem.");
            	return;
        	}

        	/* Get the status of the found solution and check if the optimization was
               successful. */
        	solstat = CPXgetstat(env, lp);


        	if (solstat == 0) {
            	fprintf (stderr, "Solution status is 0 which indicates an error\n");
            	*status = 1;
            	return;
        	}

        	/* write problem */

        	*status = CPXwriteprob(env, lp, "postrelax.lp", "LP");

        	if (solstat == 3) {
        		printf("CPLEX cannot construct feasible problem\n");
        		printf("Set cplex_method to 1,2 or 3.\n");
        		printf("Press Enter:\n");
        		getchar();
        		*status = 1;
        		return;
        	}

        	if (solstat == 1) {
        		// printf("foundsolution is feasible, continue");
        		feasible = 10;
        	}
        }
    }

    /* Get the optimization results from the CPLEX environment. First, get the
       minimizer.*/
    *status = CPXgetx(env, lp, xopt, 0, nvar-1);

    if (*status) {
    	print_error(*status, "Failed to obtain xopt from CPLEX.");
        return;
    }

    /* update lagrange multipliers if feasible subproblem is solved. */
    *status = CPXgetpi(env, lp, lambda, 0, nc-1);

    if (*status) {
    	print_error(*status, "Failed to obtain Lagrange multipliers form CPLEX");
        return;
    }
    //fprintf(stderr, "Time in solve_qps: %f", (t_e - t_b)/CLOCKS_PER_SEC );
    return;
}


/**
 *
 *
 *
 *
 *
 *
 *
 */
void cplex_solve_qps_(int *status) {
	
	
	
}

void cplex_mod_qc(int *status, const double new_val, const int prop,const int con, const  int n, const int ni)
{
	int linnnz, quadnnz, linspace, quadspace;
	double rhs;
	char sense;

	int linsurp;
        int quadsurp;

	int linind[n], quadrow[n], quadcol[n];
	double linval[n], quadval[n];


	//Remember constraint spaces....
	linspace = 1;
	quadspace = 1;

	 linind[linspace];
	 quadrow[quadspace];
         quadcol[quadspace];

	 int num = CPXgetnumqconstrs(env, lp);
	 fprintf(stderr, "num: %d\n", num);
	*status = CPXgetqconstr (env, lp, &linnnz, &quadnnz, &rhs, &sense, linind, linval , linspace, &linsurp, quadrow, quadcol, quadval, quadspace, &quadsurp, con);
        if (*status){ fprintf(stderr, "error 1 %d", *status);}


	*status =  CPXdelqconstrs(env, lp, con, con);

	fprintf(stderr, "rhs: %f\n", rhs);
	fprintf(stderr, "newval: %f\n", new_val);
	if (prop == 1) //change rhs
	{
        	*status = CPXaddqconstr (env, lp, linnnz, quadnnz, new_val, sense, linind, linval, quadrow, quadcol, quadval, NULL);
       		 if (*status){ fprintf(stderr, "error 2");}

		 //Remember to switch "remember arrays....." put last move up
	}

}

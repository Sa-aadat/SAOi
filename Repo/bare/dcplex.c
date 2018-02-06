/*=============================================================================

@author: C.G.P. Cleven
@license: ??
@date: 24 November 2011
@version: 0.7.11


This file contains dummy C routines of CPLEX such that the SAOi algorithm can
be used whithout cplex.

1. print_error
2. cplex_init_
3. cplex_stop_
4. cplex_solve_sqp_

=============================================================================*/
#include <stdio.h>

void CPXgeterrorstring()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXopenCPLEX()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXcreateprob()
{
	int *stat;
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXcloseCPLEX()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXcopylp()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXchgobj()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXchgrhs()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXchgcoef()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXchgbds()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXcopyqpsep()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXsetintparam()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXsetdblparam()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXfeeprob()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXqpopt()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXgetstat()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXfeasopt()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXgetnumrows()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXgetrowinfeas()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXgetx()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXwriteprob()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXfreeprob()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

void CPXgetpi()
{
    fprintf (stderr, "Terminal: This is a CPLEX dummy routine....\n");
}

/******************************************************************
 * File          : pvbbdpre.h                                     *
 * Programmers   : Michael Wittman and Alan C. Hindmarsh @ LLNL   *
 * Last Modified : 3 August 1998                                  *
 *----------------------------------------------------------------*
 * This is the header file for the PVBBDPRE module, for a         *
 * band-block-diagonal preconditioner, i.e. a block-diagonal      *
 * matrix with banded blocks, for use with PVODE and CVSpgmr.     *
 *                                                                *
 * Summary:                                                       *
 *                                                                *
 * These routines provide a preconditioner matrix for PVODE that  *
 * is block-diagonal with banded blocks.  The blocking corresponds*
 * to the distribution of the dependent variable vector y among   *
 * the processors.  Each preconditioner block is generated from   *
 * the Jacobian of the local part (on the current processor) of a *
 * given function g(t,y) approximating f(t,y).  The blocks are    *
 * generated by a difference quotient scheme on each processor    *
 * independently.  This scheme utilizes an assumed banded         *
 * structure with given half-bandwidths, mudq and mldq.           *
 * However, the banded Jacobian block kept by the scheme has      *
 * half-bandwiths mukeep and mlkeep, which may be smaller.        *
 *                                                                *
 * The user's calling program should have the following form:     *
 *                                                                *
 *   #include "pvbbdpre.h"                                        *
 *   ...                                                          *
 *   PVBBDData p_data;                                            *
 *   ...                                                          *
 *   machEnv = PVecInitMPI(...);                                  *
 *   ...                                                          *
 *   cvode_mem = CVodeMalloc(...);                                *
 *   ...                                                          *
 *   p_data = PVBBDAlloc(Nlocal, mudq ,mldq, mukeep, mlkeep,      *
 *              dqrely, gloc,cfn, f_data);                        *
 *   ...                                                          *
 *   CVSpgmr(cvode_mem, pretype, gstype, maxl, delt, PVBBDPrecon, *
 *           PVBBDPSol, p_data);                                  *
 *   ...                                                          *
 *   ier = CVode(...);                                            *
 *   ...                                                          *
 *   PVBBDFree(p_data);                                           *
 *   ...                                                          *
 *   CVodeFree(...);                                              *
 *                                                                *
 *   PVecFreeMPI(machEnv);                                        *
 *                                                                *
 *                                                                *
 * The user-supplied routines required are:                       *
 *                                                                *
 *   f      = function defining the ODE right-hand side f(t,y).   *
 *                                                                *
 *   gloc = function defining the approximation g(t,y).           *
 *                                                                *
 *   cfn = function to perform communication need for gloc.       *
 *                                                                *
 *                                                                *
 * Notes:                                                         *
 *                                                                *
 * 1) This header file is included by the user for the definition *
 *    of the PVBBDData type and for needed function prototypes.   *
 *                                                                *
 * 2) The PVBBDAlloc call includes half-bandwiths mudq and mldq   *
 *    to be used in the difference-quotient calculation of the    *
 *    approximate Jacobian.  They need not be the true            *
 *    half-bandwidths of the Jacobian of the local block of g,    *
 *    when smaller values may provide a greater efficiency.       *
 *    Also, the half-bandwidths mukeep and mlkeep of the retained *
 *    banded approximate Jacobian block may be even smaller,      *
 *    to reduce storage and computation costs further.            *
 *    For all four half-bandwidths, the values need not be the    *
 *    same on every processor.                                    *
 *                                                                *
 * 3) The actual name of the user's f function is passed to       *
 *    CVodeMalloc, and the names of the user's gloc and cfn       *
 *    functions are passed to PVBBDAlloc.                         *
 *                                                                *
 * 4) The pointer to the user-defined data block f_data, which is *
 *    passed to CVodeMalloc, is also passed to PVBBDAlloc, and    *
 *    is available to the user in gloc and cfn.                   *
 *                                                                *
 * 5) For the CVSpgmr call, the preconditioning and Gram-Schmidt  *
 *    types, pretype and gstype, are left to the user to specify. *
 *                                                                *
 * 6) Functions PVBBDPrecon and PVBBDPSol are never called by the *
 *    user explicitly; their names are simply passed to CVSpgmr.  *
 *                                                                *
 * 7) Optional outputs specific to this module are available by   *
 *    way of macros listed below.  These include work space sizes *
 *    and the cumulative number of gloc calls.  The costs         *
 *    associated with this module also include nsetups banded LU  *
 *    factorizations, nsetups cfn calls, and nps banded           *
 *    backsolve calls, where nsetups and nps are CVODE optional   *
 *    outputs.                                                    *
 ******************************************************************/

#ifndef _pvbbdpre_h
#define _pvbbdpre_h

#include "cvode.h"
#include "llnltyps.h"
#include "nvector.h"
#include "band.h"

namespace pvode {


/******************************************************************
 * Type : PVLocalFn                                               *
 *----------------------------------------------------------------*        
 * The user must supply a function g(t,y) which approximates the  *
 * right-hand side function f for the system y'=f(t,y), and which *
 * is computed locally (without inter-processor communication).   *
 * (The case where g is mathematically identical to f is allowed.)*
 * The implementation of this function must have type PVLocalFn.  *
 *                                                                *
 * This function takes as input the local vector size Nlocal, the *
 * independent variable value t, the local real dependent         *
 * variable array ylocal, and a pointer to the user-defined data  *
 * block f_data.  It is to compute the local part of g(t,y) and   *
 * store this in the real array glocal.                           *
 * (Allocation of memory for ylocal and glocal is handled within  *
 * the preconditioner module.)                                    *
 * The f_data parameter is the same as that passed by the user    *
 * to the CVodeMalloc routine.                                    *
 * A PVLocalFn gloc does not have a return value.                 *
 ******************************************************************/

typedef void (*PVLocalFn)(integer Nlocal, real t, real *ylocal,
                          real *glocal, void *f_data);


/******************************************************************
 * Type : PVCommFn                                                *
 *----------------------------------------------------------------*        
 * The user must supply a function of type PVCommFn which performs*
 * all inter-processor communication necessary to evaluate the    *
 * approximate right-hand side function described above.          *
 *                                                                *
 * This function takes as input the local vector size Nlocal,     *
 * the independent variable value t, the dependent variable       *
 * vector y, and a pointer to the user-defined data block f_data. *
 * The f_data parameter is the same as that passed by the user to *
 * the CVodeMalloc routine.  The PVCommFn cfn is expected to save *
 * communicated data in space defined with the structure *f_data. *
 * A PVCommFn cfn does not have a return value.                   *
 *                                                                *
 * Each call to the PVCommFn cfn is preceded by a call to the     *
 * RhsFn f with the same (t,y) arguments.  Thus cfn can omit any  *
 * communications done by f if relevant to the evaluation of g.   *
 ******************************************************************/

typedef void (*PVCommFn)(integer Nlocal, real t, N_Vector y,
                         void *f_data);

 
/*********************** Definition of PVBBDData *****************/

typedef struct pvbbddata_type {

  /* passed by user to PVBBDAlloc, used by Precond/Psolve functions: */
  void *f_data;
  integer mudq, mldq, mukeep, mlkeep;
  real dqrely;
  PVLocalFn gloc;
  PVCommFn cfn;

  /* set by PVBBDPrecon and used by PVBBDPSol: */
  BandMat savedJ;
  BandMat savedP;
  integer *pivots;

  /* available for optional output: */
  integer rpwsize;
  integer ipwsize;
  integer nge;

} *PVBBDData;


/*************** Macros for optional outputs **********************
 *                                                                *
 * PVBBD_RPWSIZE(pdata) returns the size of the real work space,  *
 * in real words, used by this preconditioner module.             *
 * This size is local to the current processor.                   *
 *                                                                *
 * PVBBD_IPWSIZE(pdata) returns the size of the integer work      *
 * space, in integer words, used by this preconditioner module.   *
 * This size is local to the current processor.                   *
 *                                                                *
 * PVBBD_NGE(pdata) returns the number of g(t,y) evaluations,     *
 * i.e. the number of calls to the gloc function, so far.         *
 ******************************************************************/

#define PVBBD_RPWSIZE(pdata) (pdata->rpwsize)
#define PVBBD_IPWSIZE(pdata) (pdata->ipwsize)
#define PVBBD_NGE(pdata) (pdata->nge)


/******************************************************************
 * Function : PVBBDAlloc                                          *
 *----------------------------------------------------------------*
 * PVBBDAlloc allocates and initializes a PVBBDData structure     *
 * to be passed to CVSpgmr (and subsequently used by PVBBDPrecon  *
 * and PVBBDPSol.                                                 *
 *                                                                *
 * The parameters of PVBBDAlloc are as follows:                   *
 *                                                                *
 * Nlocal  is the length of the local block of the vectors y etc. *
 *         on the current processor.                              *
 *                                                                *
 * mudq, mldq  are the upper and lower half-bandwidths to be used *
 *         in the difference-quotient computation of the local    *
 *         Jacobian block.                                        *
 *                                                                *
 * mukeep, mlkeep  are the upper and lower half-bandwidths of the *
 *         retained banded approximation to the local Jacobian    * 
 *         block.                                                 *
 *                                                                *
 * dqrely  is an optional input.  It is the relative increment    *
 *         in components of y used in the difference quotient     *
 *         approximations.  To specify the default, pass 0.       *
 *         The default is dqrely = sqrt(unit roundoff).           *
 *                                                                *
 * gloc    is the name of the user-supplied function g(t,y) that  *
 *         approximates f and whose local Jacobian blocks are     *
 *         to form the preconditioner.                            *
 *                                                                *
 * cfn     is the name of the user-defined function that performs *
 *         necessary inter-processor communication for the        *
 *         execution of gloc.                                     *
 *                                                                *
 * f_data  is a pointer to the optional user data block.          *
 *                                                                *
 * PVBBDAlloc returns the storage allocated (type PVBBDData),     *
 * or NULL if the request for storage cannot be satisfied.        *
 ******************************************************************/

PVBBDData PVBBDAlloc(integer Nlocal, integer mudq, integer mldq,
            integer mukeep, integer mlkeep, real dqrely,
            PVLocalFn gloc, PVCommFn cfn, void *f_data );


/******************************************************************
 * Function : PVBBDFree                                           *
 *----------------------------------------------------------------*
 * PVBBDFree frees the memory block p_data allocated by the call  *
 * to PVBBDAlloc.                                                 *
 ******************************************************************/

void PVBBDFree(PVBBDData pdata);


/****** Prototypes of functions PVBBDPrecon and PVBBDPSol *********/
  
int PVBBDPrecon(integer N, real t, N_Vector y, N_Vector fy, boole jok,
	      boole *jcurPtr, real gamma, N_Vector ewt, real h, real uround,
              long int *nfePtr, void *P_data, N_Vector vtemp1,
	      N_Vector vtemp2, N_Vector vtemp3);

int PVBBDPSol(integer N, real t, N_Vector y, N_Vector fy, N_Vector vtemp,
	    real gamma, N_Vector ewt, real delta, long int *nfePtr,
	    N_Vector r, int lr, void *P_data, N_Vector z);

}
#endif

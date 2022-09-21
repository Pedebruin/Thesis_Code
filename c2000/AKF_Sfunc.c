/*
 * Code for block defenition of the AKF block.
 *
 * Autor: Pim de Bruin
 * Date: 10/05/2022
 * Delft University of Technology
 */

#define S_FUNCTION_NAME  AKF_Sfunc
#define S_FUNCTION_LEVEL 2

/*
 * Need to include simstruc.h for the definition of the SimStruct and
 * its associated macro definitions.
 */
#include "simstruc.h"
#include "MatLib\MatLib.c"

#include <stdio.h>
#include <stdlib.h>

#define y(element) (*yPtrs[element])  /* Pointer to Input Port0 */

#define A_IDX 0
#define A_PARAM(S) ssGetSFcnParam(S,A_IDX)

#define C_IDX 1
#define C_PARAM(S) ssGetSFcnParam(S,C_IDX)

#define Q_IDX 2
#define Q_PARAM(S) ssGetSFcnParam(S,Q_IDX)

#define R_IDX 3
#define R_PARAM(S) ssGetSFcnParam(S,R_IDX)

#define N_IDX 4
#define N_PARAM(S) ssGetSFcnParam(S,N_IDX)

#define P0_IDX 5
#define P0_PARAM(S) ssGetSFcnParam(S,P0_IDX)

#define Bw_IDX 6
#define Bw_PARAM(S) ssGetSFcnParam(S,Bw_IDX)

unsigned int nu = 1;
unsigned int ny = 1;
unsigned int nq = 5; 

Mat* A;
Mat* C; 
Mat* Q;
Mat* R; 
Mat* N; 
Mat* Bw;

/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 10);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, nq+nq*nq);
    /* The state vector for this block contains of [q; vec(P)], where q = [x;u]
    *, and vec(P) is the vectorised version of the P matrix. 
    */

    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, ny);
//     ssSetInputPortRequiredContiguous(S, 0, 1); /*direct input signal access*/
    ssSetInputPortDirectFeedThrough(S, 0, 1); // input 0 has direct feedthrough

    if (!ssSetNumOutputPorts(S, 3)) return;
    ssSetOutputPortWidth(S, 0, ny);     // yhat
    ssSetOutputPortWidth(S, 1, nq-nu);  // xhat
    ssSetOutputPortWidth(S, 2, nu);     // uhat
    
    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);
    ssSetOperatingPointCompliance(S, USE_DEFAULT_OPERATING_POINT);

    ssSetRuntimeThreadSafetyCompliance(S, RUNTIME_THREAD_SAFETY_COMPLIANCE_TRUE);
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

// Specify sample time
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

// Initialise states every time states are reset
#define MDL_INITIALIZE_CONDITIONS 
#if defined(MDL_INITIALIZE_CONDITIONS)
  static void mdlInitializeConditions(SimStruct *S)
  {    
      // Get the P0 parameter
      real_T    *P0_PARAM_ARRAY      = mxGetPr(P0_PARAM(S));
    
      // Get pointer to state vector
      real_T *states0 = ssGetRealDiscStates(S);
        
      // Initialise states to zero
      int_T lp;
      for (lp = 0;lp<nq-nu;lp++){
          states0[lp] = 0;
      }

      // Initialise Input to zero
      for (lp = nq-nu; lp<nq; lp++){
          states0[lp] = 0;
      }
        
      // Initialise the rest at the P0 values (states = [x;u;P], q = [x;u])
      lp = nq;
      for (lp = nq;lp<nq+nq*nq;lp++){
            states0[lp] = (real_T) P0_PARAM_ARRAY[lp-nq];
      }
 }
#endif /* MDL_INITIALIZE_CONDITIONS */

// Initialise once at start of simulation
#define MDL_START
#if defined(MDL_START) 
  static void mdlStart(SimStruct *S)
  {
      real_T    *A_PARAM_ARRAY      = mxGetPr(A_PARAM(S));
      real_T    *C_PARAM_ARRAY      = mxGetPr(C_PARAM(S));
      real_T    *Q_PARAM_ARRAY      = mxGetPr(Q_PARAM(S));
      real_T    *R_PARAM_ARRAY      = mxGetPr(R_PARAM(S));
      real_T    *N_PARAM_ARRAY      = mxGetPr(N_PARAM(S));
      real_T    *Bw_PARAM_ARRAY      = mxGetPr(Bw_PARAM(S));

      A = from_array(nq,nq,A_PARAM_ARRAY);
      C = from_array(ny,nq,C_PARAM_ARRAY);
      Q = from_array(nq,nq,Q_PARAM_ARRAY);
      R = from_array(ny,ny,R_PARAM_ARRAY);
      N = from_array(nq,ny,N_PARAM_ARRAY);
      Bw = from_array(nq,nq,Bw_PARAM_ARRAY);
}
#endif /*  MDL_START */

// Compute outputs of the block
static void mdlOutputs(SimStruct *S, int_T tid)
{
        UNUSED_ARG(tid);

        // Get pointer to block output vectors (need to put values in here later)
        real_T              *yhat_array = ssGetOutputPortRealSignal(S,0);
        real_T              *xhat_array = ssGetOutputPortRealSignal(S,1);
        real_T              *uhat_array = ssGetOutputPortRealSignal(S,2);
    
        // Get states and system output
        real_T              *states_array = ssGetRealDiscStates(S);         // Get states vector
        
        // Split states_array into q_array and u_array
        real_T              *q_array = states_array;
        real_T              *u_array = states_array + nq-nu; 
               
        // Get state and input estimate in matrix
        Mat *xhat = from_array(nq-nu,   1,  q_array);
        Mat *uhat = from_array(nu,      1,  u_array);
        
        // compute the output in matrix
//         Mat *yhat = multiply(C,xhat);
        Mat *yhat = ones(ny,ny);
        
        // Write content of matrices to appropriate output arrays 
        to_array(xhat, xhat_array);
        to_array(yhat, yhat_array);
        to_array(uhat, uhat_array);

//         freemat(xhat); 
//         freemat(yhat);
//         freemat(uhat);
}

// Update discrete states
#define MDL_UPDATE
#if defined(MDL_UPDATE)
  static void mdlUpdate(SimStruct *S, int_T tid)
  {
        UNUSED_ARG(tid);
        InputRealPtrsType   yPtrs = ssGetInputPortRealSignalPtrs(S,0); // Get sys output
        real_T              *states = ssGetRealDiscStates(S);          // Get states vector
        
        // Separate q and P, and make matrix out of P
        real_T *q_array = states;
        real_T *P_array = states + nq;
        
        // Get q1hat, P1 and Y matrix (q(k-1|k-1), P(k-1|k-1))
        Mat *q1hat   =   from_array(nq, 1, q_array);
        Mat *P1      =   from_array(nq, nq, P_array);
        
        // Get Y matrix (y(k))
        double ytemp[1] = {y(0)};
        Mat *Y       =   from_array(ny, 1,ytemp);  

        /* Actual algorithm =============================================*/
        // Often used matrices
        Mat *CT = transpose(C);
                   
        // Tim update
        Mat *qhatkk1 = multiply(A,q1hat);            // q(k|k-1) = Aq(k-1)
        Mat *Pkk1 = sum(   	                    // P(k|k-1) = AP(k-1|k-1)A' + BwQBw'
                                multiply(multiply(A,P1),transpose(A)),
                                multiply(multiply(Bw,Q),transpose(Bw)));  
        
        // Measurement update
        Mat *K = multiply(   multiply(Pkk1,CT),                             // K(k) = P(k|k-1)CT/(CP(k|k-1)CT+R)
                                    inverse(sum(multiply(C,multiply(Pkk1,CT)),R)));
    
        Mat *qhatkk = sum(qhatkk1,              // qhat(k|k) = qhat(k|k-1) + K(k)*(y(k) - C*qhat(k|k-1))
                                      multiply(K,
                                                  minus(Y,multiply(C,qhatkk1))));
        Mat *Pkk = minus(Pkk1,                    // P(k|k-1) - K(k)CP(k|k-1)
                                      multiply(K,
                                                  multiply(C,Pkk1)));
                
        /* ==============================================================*/ 
        // Write states into the simulink state vector
        to_array(qhatkk,  q_array);
        to_array(Pkk,     P_array);

        //  Clean
        freemat(q1hat);
        freemat(P1);
        freemat(Y);
        freemat(CT);
        freemat(qhatkk1);
        freemat(Pkk1);
        freemat(K);
        freemat(qhatkk);
        freemat(Pkk);
  }

#endif /* MDL_UPDATE */

// Clean up
static void mdlTerminate(SimStruct *S)
{
    freemat(A);
    freemat(C);
    freemat(Q);
    freemat(R);
    freemat(N);
    freemat(Bw);
}


/*=============================*
 * Required S-function trailer *
 *=============================*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

// NUMERICAL SIMULATION OF MSO NEURON
// MICHIEL REMME
// OCTOBER 2018

#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double* Make1DDoubleArray(int arraySizeX) {
    double* theArray;
    int i;
    theArray = (double*) mxCalloc(arraySizeX,sizeof(double));
    if (theArray == NULL) mexErrMsgTxt("Cannot allocate temporary variables\n");
    return theArray;
}

void tri( int matrix_size, double *A, double *D, double *C, double *B, double *X) {
    int i;
    double xmult;
    
    for (i = 1; i < matrix_size; i++) {
        xmult = A[i-1]/D[i-1];
        D[i] -= xmult * C[i-1];
        B[i] -= xmult * B[i-1];
    }
    X[matrix_size-1] = B[matrix_size-1]/D[matrix_size-1];
    for (i = matrix_size -2; i>=0; i--)   X[i] = (B[i]- C[i]* X[i+1]) / D[i];
}

double winf(double v) {
    double whalf   = -57.34;
    double wk      = -11.7;
    return 1/(1+exp((v-whalf)/wk));
}

double wtau(double v, double tw_fac) {
    return (100 /(6*exp((v+60)/7) + 24*exp((v+60)/-51)) + 1.59) * 0.22 * tw_fac;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    ////////// INPUT PARAMETERS FROM MATLAB
	double *simparams, *cellparams, *locvec, *gvec;
    int		i;

    simparams       = mxGetPr(prhs[0]);
    cellparams      = mxGetPr(prhs[1]);
    locvec          = mxGetPr(prhs[2]); // vector of input locations
    gvec            = mxGetPr(prhs[3]); // input conductance matrix

    int     tend    = simparams[0];    // (ms)
    double  dt      = simparams[1];    // (ms)
    double  dx      = simparams[2];    // (elect size) spatial discretization
    double  L       = cellparams[0];
    double  rho     = cellparams[1];
    double  tm      = cellparams[2];
    double  gw      = cellparams[3]; // ratio peak density to gleak
    double  tw_fac  = cellparams[4]; // multiplication factor for tauw
    double  Vrest   = cellparams[5];
    double  Ena     = cellparams[6]; // sodium reversal potential
    double  Ek      = cellparams[7]; // potassium reversal potential
    double  Es      = cellparams[8]; // reversal potential of synaptic input
    double  RNinf   = cellparams[9]; // RNinf semi-infinite cable
    double  sdAreaRatio = cellparams[10]; // how much larger is soma area compared to single compartment? 
    double  gwzinf  = gw*((1-0.27)/(1+exp((Vrest+67)/6.16))+0.27);
    double  El      = gwzinf*pow(winf(Vrest),4)*(Vrest-Ek) + Vrest;  // leak reversal potential
    int     Ninputs = mxGetNumberOfElements(prhs[2]);  // number of inputs
    int     iloc;        
    int     nseg    = (int)floor(L/dx);     // number of compartments per cable
    int     Nc      = 2*nseg + 1;           // total number of compartments
    int     Nt      = (int)floor(tend/dt);  // number of iterations
    
    for (i=0;i<Ninputs;i++) {
        if (locvec[i]<0||locvec[i]>=Nc) mexErrMsgTxt("Trying to address non-existent cable location\n");
    }
    if (nseg<3) mexErrMsgTxt("Too few compartments\n"); 
    
    ////////// DECLARE PARAMETERS AND VARIABLES    
    int     t;
    double  s  = dx*dx*tm/dt;
    double  r  = 2*(1 + s + dx*dx);
    double  ss = dx*tm/(dt*rho);
    double  rs = 2 + ss + dx/rho;
            
    double*	A  = Make1DDoubleArray(Nc);   // below diagonal
    double*	C  = Make1DDoubleArray(Nc);   // above diagonal
    double*	D  = Make1DDoubleArray(Nc);   // diagonal
    double*	U  = Make1DDoubleArray(Nc);   // voltage new
    double*	V  = Make1DDoubleArray(Nc);   // voltage old
    
    double*	W  = Make1DDoubleArray(Nc);   // gating variable    
    double*	gklt = Make1DDoubleArray(Nc); // klt conductance divided by gleak
    
	////////// OUTPUT VARIABLES TO MATLAB
    double *vs_result, *vd_result, *ina_syn_result, *ik_syn_result, *ina_leak_result, *ik_leak_result, *ik_klt_result;
    for (i=0;i<7;i++) {
        plhs[i]  = mxCreateDoubleMatrix(1, Nt, mxREAL);
        if (plhs[i] == NULL) mexErrMsgTxt("Cannot allocate output arrays\n");
    }        
    vs_result           = mxGetPr(plhs[0]);    
    vd_result           = mxGetPr(plhs[1]);    
    ina_syn_result      = mxGetPr(plhs[2]);    
    ik_syn_result       = mxGetPr(plhs[3]);   
    ina_leak_result     = mxGetPr(plhs[4]);    
    ik_leak_result      = mxGetPr(plhs[5]);    
    ik_klt_result       = mxGetPr(plhs[6]);
    
    ////////// INITIALIZE    
    for (i = 0; i<Nc-1; i++) {
        A[i] = -1.0;
        D[i] = r;
        C[i] = -1.0;
    }
    D[0]     = r-1; // BC at x=0
    D[Nc-1]  = r-1; // BC at x=end
    
    D[nseg-1]= r+1; // BC at transition cable1-soma
    C[nseg-1]= -2.0;

    D[nseg]  = rs;  // BC at soma

    A[nseg]  = -2.0;
    D[nseg+1]= r+1; // BC at transition soma-cable2

    for (i=0; i<Nc; i++) {
        U[i] = Vrest;
        W[i] = winf(Vrest);
    }
        
    ////////// SIMULATE
    for (t = 0; t < Nt; t++) {
        // compute gklt/gpas in all compartments
        for (i = 0; i<Nc; i++) {
            gklt[i] =  gwzinf*pow(W[i],4);
        }
        
        // UPDATE MATRIX
        D[0]        = r-1 + 2*dx*dx*gklt[0];
        V[0]        = (2*s-1)*U[0] + U[1] + 2*dx*dx*(El + gklt[0]*Ek);
        for (i = 1; i<Nc-1; i++) {
            D[i]    = r + 2*dx*dx*gklt[i];
            V[i]    = U[i-1] + 2*(s-1)*U[i] + U[i+1] + 2*dx*dx*(El + gklt[i]*Ek);
        }
        D[Nc-1]     = r-1 + 2*dx*dx*gklt[Nc-1];
        V[Nc-1]     = U[Nc-2] + (2*s-1)*U[Nc-1] + 2*dx*dx*(El + gklt[Nc-1]*Ek);
        
        // BC at transition cable1-soma
        D[nseg-1]   = r+1 + 2*dx*dx*gklt[nseg-1];
        V[nseg-1]   = U[nseg-2] +(2*s-3)*U[nseg-1]  + 2*U[nseg] + 2*dx*dx*(El + gklt[nseg-1]*Ek);
        // BC at soma
        D[nseg]     = rs + dx/rho*gklt[nseg];
        V[nseg]     = U[nseg-1] +(ss-2)*U[nseg]     + U[nseg+1] + dx/rho*(El + gklt[nseg]*Ek);
        // BC at transition soma-cable2
        D[nseg+1]   = r+1 + 2*dx*dx*gklt[nseg+1];
        V[nseg+1]   = 2*U[nseg] +(2*s-3)*U[nseg+1]  + U[nseg+2] + 2*dx*dx*(El + gklt[nseg+1]*Ek);
        
        // EXTERNAL INPUT
        for (i = 0; i<Ninputs; i++) {
            iloc = (int)locvec[i];
            // conductance input
            D[iloc]    += RNinf*2*dx*gvec[Ninputs*t+i];
            V[iloc]    += RNinf*2*dx*gvec[Ninputs*t+i]*Es; // current is point process
        }
        
        tri(Nc, A, D, C, V, U); // compute new voltage U[i]
        
        for (i = 0; i<Nc; i++) { // update gating variable
            W[i] += (1 - exp(-dt/wtau(U[i],tw_fac)))*(winf(U[i])-W[i]);
        }
        
   		// store variables in the result array
        vs_result[t] = U[nseg]; // somatic voltage
        vd_result[t] = U[1]; // dendritic voltage
        for (i = 0; i<Ninputs; i++) {
            iloc = (int)locvec[i];
            ina_syn_result[t] += 2.0/3*gvec[Ninputs*t+i]*(U[iloc]-Ena);
            ik_syn_result[t]  += 1.0/3*gvec[Ninputs*t+i]*(U[iloc]-Ek);
        }
        for (i=0; i<Nc; i++) { // the following currents still have to be multiplied by gpas
            ina_leak_result[t]+= (U[i]-Ena);    
            ik_leak_result[t] += (U[i]-Ek);    
            ik_klt_result[t]  += gklt[i]*(U[i]-Ek);
        }
        // correct for greater soma area than dend compartments
        ina_leak_result[t]+= (sdAreaRatio-1)*(U[nseg]-Ena);    
        ik_leak_result[t] += (sdAreaRatio-1)*(U[nseg]-Ek);    
        ik_klt_result[t]  += (sdAreaRatio-1)*gklt[nseg]*(U[nseg]-Ek);        
    }
    
    mxFree(A);
    mxFree(C);
    mxFree(D);
    mxFree(U);
    mxFree(V);
    mxFree(W);
    mxFree(gklt);

	return;
}

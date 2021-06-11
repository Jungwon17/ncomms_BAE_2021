#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nbins;
    double jitter;
    
    double *C;
    
    jitter = mxGetScalar(prhs[0]);
    nbins = (int)mxGetScalar(prhs[1]);
        
    plhs[0] = mxCreateDoubleMatrix(nbins, 1, mxREAL);
    C = mxGetPr(plhs[0]);
    
    for (int i = 0; i < nbins; i++) {
        C[i] = jitter*((double)rand()/32767*2 - 1);
    }
}

#include "mex.h"

// input: nx, ny
// output: nx * ny matrix

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nx, ny;
    double *C;
    int k = 0;
    
    nx = (int)mxGetScalar(prhs[0]);
    ny = (int)mxGetScalar(prhs[1]);
        
    plhs[0] = mxCreateDoubleMatrix(nx, ny, mxREAL);
    C = mxGetPr(plhs[0]);
    
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            k+=1;
            C[ny*ix + iy] = k;
        }
    }
}
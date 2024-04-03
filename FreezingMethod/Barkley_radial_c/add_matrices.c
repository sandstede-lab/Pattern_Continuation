// add_matrices.c
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 
/* void mexFunction(int nlhs, mxArray *plhs[], mxArray *prhs[]) { */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("add_matrices:nrhs", "Two input matrices are required.");
    }
    if (nlhs != 1) {
         mexErrMsgIdAndTxt("add_matrices:nlhs", "One output matrix is required.");
    }

    double *matrixA = mxGetPr(prhs[0]);
    double *matrixB = mxGetPr(prhs[1]);
    mwSize numRows = mxGetM(prhs[0]);
    mwSize numCols = mxGetN(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix(numRows, numCols, mxREAL);
    double *resultMatrix = mxGetPr(plhs[0]);

    for (mwSize i = 0; i < numRows * numCols; ++i) {
        resultMatrix[i] = matrixA[i] + matrixB[i];
    }
}

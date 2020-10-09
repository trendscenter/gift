/************************************************************************************************************************
		Function to compute all eigen values and vectors of real symmetric matrix using packed storage method in LAPACK
*************************************************************************************************************************/

#include "mex.h"
#include "string.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Declare Variables */
	char *ch1 = "V", *ch2 = "L";	
	void *A, *W, *Z, *work;
	
#ifdef MWSIZE
	mwSize N, numElem, *iwork, lwork, liwork, info = 0, optional, isDouble = 1;
#else
	int N, numElem, *iwork, lwork, liwork, info = 0, optional, isDouble = 1;
#endif

	/****** Do Error check *************/
	if (nrhs != 3) 
	{
		mexErrMsgTxt("icatb_eig_symm_all function accepts three input arguments.\nThe usage of the function is [eigen_vectors, eigen_values] = icatb_eig_symm_all(A, N, opt);\nwhere A is the lower triangular matrix without zeros entered as a vector, N is the original dimension and opt is parameter to operate on copy of the matrix or the original array.");
	}
	
	/* Inputs Matrix (A), order (N)  and parameter (opt) to operate on copy or original matrix*/
#ifdef MWSIZE
	numElem = (mwSize) mxGetNumberOfElements(prhs[0]);   
	N = (mwSize) mxGetScalar(prhs[1]);
	/* Parameter for using the copy or original array to do operations*/
	optional = (mwSize) mxGetScalar(prhs[2]);
#else
	numElem = (int) mxGetNumberOfElements(prhs[0]);
	N = (int) mxGetScalar(prhs[1]);
	optional = (int) mxGetScalar(prhs[2]);	
#endif
		
	
	if (numElem != (N*(N + 1)/2))
	{
		mexErrMsgTxt("Number of elements of input matrix (A) must be equal to N*(N+1)/2");
	}
	
	if (nlhs > 2)
	{
		mexErrMsgTxt("Maximum allowed output arguments is 2.");
	}
	/********* End for doing error check *********/
	
	/* Work array dimensions */
	lwork = 1 + 6*N + N*N;
	liwork = 3 + 5*N;
	
	
	/* Check isdouble */
	if (mxIsDouble(prhs[0]))
	{		
		isDouble = 1;
		
		if (optional == 1) 
		{			
			A = (void *) mxCalloc(numElem, sizeof(double));
			memcpy(A, mxGetData(prhs[0]), numElem*sizeof(double));   			
		}
		else
		{		
			A = mxGetData(prhs[0]);
		}		

		A = (double *) A;
		
		/* Pointer to work */
		work = (double *) mxCalloc(lwork, sizeof(double));		
		
		/* Pointer to eigen vectors */
		plhs[0] = mxCreateDoubleMatrix(N, N, mxREAL);
		Z = mxGetPr(plhs[0]);			
		
		/* Pointer to eigen values */
		plhs[1] = mxCreateDoubleMatrix(1, N, mxREAL);
		W = mxGetPr(plhs[1]);		
	
	}

	else
	{
		isDouble = 0;
		
		if (optional == 1) 
		{			
			A = (void *) mxCalloc(numElem, sizeof(float));
			memcpy(A, mxGetData(prhs[0]), numElem*sizeof(float));   
		}
		else
		{
			A = mxGetData(prhs[0]);
		}

		A = (float *) A;
		
		/* Pointer to work */
		work = (float *) mxCalloc(lwork, sizeof(float));
		
		/* Pointer to eigen vectors */
		plhs[0] = mxCreateNumericMatrix(N, N, mxSINGLE_CLASS, mxREAL);	
		Z = (float *) mxGetData(plhs[0]);			
		
		/* Pointer to eigen values */
		plhs[1] = mxCreateNumericMatrix(1, N, mxSINGLE_CLASS, mxREAL);	
		W = (float *) mxGetData(plhs[1]);		
	}		
	/* End for check isdouble */
	
	/* Pointer to iwork */
#ifdef MWSIZE
	iwork = (mwSize *) mxCalloc(liwork, sizeof(mwSize));
#else
	iwork = (int *) mxCalloc(liwork, sizeof(int));
#endif
	
	/* Call subroutine to find eigen values and vectors */            
#ifdef PC	
	/* Handle Windows */
	if (isDouble == 1)
	{
		/* Double precision */
		dspevd(ch1, ch2, &N, A, W, Z, &N, work, &lwork, iwork, &liwork, &info);
	}
	else
	{		
		/* Use single precision routine */
#ifdef SINGLE_SUPPORT		
		sspevd(ch1, ch2, &N, A, W, Z, &N, work, &lwork, iwork, &liwork, &info);
#else
		mexErrMsgTxt("Error executing sspevd function on this MATLAB version.\nUse double precision to solve eigen value problem");
#endif
	}
	/* End for handling Windows */
#else	
	/* Handle other OS */
	if (isDouble == 1)
	{	
		/* Double precision */
		dspevd_(ch1, ch2, &N, A, W, Z, &N, work, &lwork, iwork, &liwork, &info);
	}
	else
	{		
		/* Use single precision routine */
#ifdef SINGLE_SUPPORT
		sspevd_(ch1, ch2, &N, A, W, Z, &N, work, &lwork, iwork, &liwork, &info);
#else
		mexErrMsgTxt("Error executing sspevd_ function on this MATLAB version.\nUse double precision to solve eigen value problem");
#endif
	}
	/* End for handling other OS */
#endif
	/* End for calling subroutine to find eigen values and vectors */
	
	if (info == 1)
	{
		mexErrMsgTxt("Error executing dspevd/sspevd function.");
	}
	
	/* Free memory */
	if (optional == 1) 
	{
		mxFree(A);
	}	
	mxFree(work);
	mxFree(iwork);
}


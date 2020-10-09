/****************************************************************************************************************************
		Function to compute selected eigen values and vectors of real symmetric matrix using packed storage method in LAPACK
*****************************************************************************************************************************/

#include "mex.h"
#include "string.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Declare Variables */
	char *ch1 = "V", *ch2 = "I", *ch3 = "L";	
	void *A, *W, *Z, *work;
	double ABSTOLD = 0.0, VLD = 0.0, VUD = 0.0;
	float ABSTOLF = 0.0, VLF = 0.0, VUF = 0.0;
	
#ifdef MWSIZE
	mwSize N, *iwork, *ifail, info = 0, numElem, optional, IL, IU, M, num_eigs, isDouble = 1;
#else
	int N, *iwork, *ifail, info = 0, numElem, optional, IL, IU, M, num_eigs, isDouble = 1;
#endif	

	/****** Do Error check *************/
	if (nrhs != 4) 
	{
		mexErrMsgTxt("icatb_eig_symm_sel function accepts four input arguments.\nThe usage of the function is [eigen_vectors, eigen_values] = icatb_eig_symm_sel(A, N, num_eigs, opt);\nwhere A is the lower triangular matrix without zeros entered as a vector, N is the original dimension, num_eigs is the number of eigen values and opt is parameter to operate on copy of the matrix or the original array.");
	}
	
#ifdef MWSIZE	
	/* Inputs Matrix (A), order (N), number of eigen values (num_eigs), parameter (opt) to operate on copy or original matrix*/   
	numElem = (mwSize) mxGetNumberOfElements(prhs[0]);  
	N = (mwSize) mxGetScalar(prhs[1]); 
	num_eigs = (mwSize) mxGetScalar(prhs[2]);	
	/* Parameter for using the copy or original array to do operations*/	
	optional = (mwSize) mxGetScalar(prhs[3]);
#else
	numElem = (int) mxGetNumberOfElements(prhs[0]);
	N = (int) mxGetScalar(prhs[1]);
	num_eigs = (int) mxGetScalar(prhs[2]);	
	optional = (int) mxGetScalar(prhs[3]);	
#endif	
		
	
	if (numElem != (N*(N + 1)/2))
	{
		mexErrMsgTxt("Number of elements of input matrix (A) must be equal to N*(N+1)/2"); 
	}
		
	
	if (num_eigs > N)
	{
		mexErrMsgTxt("Number of eigen values desired is greater than the number of rows of matrix");
	}
	
	if (nlhs > 2) 
	{
		mexErrMsgTxt("Maximum allowed output arguments is 2."); 
	}    
	/********* End for doing error check *********/
	
	/* Eigen values desired */
	IL = N - num_eigs + 1;
	IU = N; 	
	
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
		work = (double *) mxCalloc(8*N, sizeof(double));		
		
		/* Pointer to eigen vectors */
		plhs[0] = mxCreateDoubleMatrix(N, num_eigs, mxREAL);
		Z = mxGetPr(plhs[0]);			
		
		/* Pointer to eigen values */
		plhs[1] = mxCreateDoubleMatrix(1, num_eigs, mxREAL);
		W = mxGetPr(plhs[1]);		

		/* Tolerance */		
		ABSTOLD = mxGetEps();		
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
		work = (float *) mxCalloc(8*N, sizeof(float));
		
		/* Pointer to eigen vectors */
		plhs[0] = mxCreateNumericMatrix(N, num_eigs, mxSINGLE_CLASS, mxREAL);	
		Z = (float *) mxGetData(plhs[0]);			
		
		/* Pointer to eigen values */
		plhs[1] = mxCreateNumericMatrix(1, num_eigs, mxSINGLE_CLASS, mxREAL);	
		W = (float *) mxGetData(plhs[1]); 
		
		ABSTOLF = (float) mxGetEps();		
	}		

	
#ifdef MWSIZE
	/* Pointer to iwork */
	iwork = (mwSize *) mxCalloc(5*N, sizeof(mwSize)); 	
	/* Pointer to ifail */
	ifail = (mwSize *) mxCalloc(N, sizeof(mwSize));		
#else	
	/* Pointer to iwork */
	iwork = (int *) mxCalloc(5*N, sizeof(int)); 	
	/* Pointer to ifail */
	ifail = (int *) mxCalloc(N, sizeof(int));		
#endif	
	
	/* Call subroutine to find eigen values and vectors */
#ifdef PC	
	/* Handle Windows */
	if (isDouble == 1)
	{	
		/* Double precision */
		dspevx(ch1, ch2, ch3, &N, A, &VLD, &VUD, &IL, &IU, &ABSTOLD, &M, W, Z, &N, work, iwork, ifail, &info);	
	}
	else
	{
		/* Use single precision routine */
#ifdef SINGLE_SUPPORT		
		sspevx(ch1, ch2, ch3, &N, A, &VLF, &VUF, &IL, &IU, &ABSTOLF, &M, W, Z, &N, work, iwork, ifail, &info);			
#else
		mexErrMsgTxt("Error executing sspevx function on this MATLAB version.\nUse double precision to solve eigen value problem");
#endif
	}
	/* End for handling Windows */
#else
	/* Handle other OS */
	if (isDouble == 1)
	{
		dspevx_(ch1, ch2, ch3, &N, A, &VLD, &VUD, &IL, &IU, &ABSTOLD, &M, W, Z, &N, work, iwork, ifail, &info);	
	}
	else
	{
		/* Use single precision routine */
#ifdef SINGLE_SUPPORT
		sspevx_(ch1, ch2, ch3, &N, A, &VLF, &VUF, &IL, &IU, &ABSTOLF, &M, W, Z, &N, work, iwork, ifail, &info);	
#else
		mexErrMsgTxt("Error executing sspevx_ function on this MATLAB version.\nUse double precision to solve eigen value problem");
#endif
	}
	/* End for handling other OS */
#endif
	/* End for calling subroutine to find eigen values and vectors */
	
	if (M != num_eigs)
	{
		mexPrintf("%s%d\t%s%d\n", "No. of eigen values estimated: ", M, "No. of eigen values desired: ", num_eigs);
		mexErrMsgTxt("Error executing dspevx function\n");
	}
	
	if (info == 1) 
	{
		mexErrMsgTxt("Error executing dspevx/sspevx function.");
	}
	
	/* Free memory */
	if (optional == 1) 
	{
		mxFree(A);
	}	
	mxFree(work);
	mxFree(iwork);
	mxFree(ifail);

}






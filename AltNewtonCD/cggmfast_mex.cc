
// Compile on OS X (10.8.2): 
//   brew install eigen
//   mex gcrf_newton.cc gcrf_newton_mex.cc -I/usr/local/include/eigen3/
// 
// Compile on Linux (Ubuntu 12.10): 
//   sudo apt-get install libeigen3-dev
//   mex gcrf_newton.cc gcrf_newton_mex.cc -I/usr/include/eigen3  

#include <vector>
#include <mex.h>
#include <Eigen/Dense>
#include "cggmfast.h"

using Eigen::Map;
using Eigen::VectorXi;
using Eigen::MatrixXd;
using Eigen::VectorXd;

enum InputArguments {
	INPUT_Y,
	INPUT_X,
	INPUT_LAMBDA_Y,
	INPUT_LAMBDA_X,
	INPUT_OPTIONS,
};

enum OutputArguments {
	OUTPUT_LAMBDA,
	OUTPUT_THETA,
	OUTPUT_HISTORY
};

// Constants needed to create the history structure
const char *kHistoryFieldNames[] = 
	{"objval", "time", "active_set_size", "subgrad", "l1norm"};
const int kNumHistoryFields = 
	sizeof(kHistoryFieldNames)/sizeof(*kHistoryFieldNames);
const mwSize kHistoryDims[5] = {1, 1, 1, 1, 1};
enum HistoryFields {
	HISTORY_OBJVAL,
	HISTORY_TIME,
	HISTORY_ACTIVE_SET_SIZE,
	HISTORY_SUBGRAD,
	HISTORY_L1NORM
};

// Functions for marshalling data to Matlab 
mxArray* getMexArray(const std::vector<double>& v) {
	mxArray* a = mxCreateDoubleMatrix(1, v.size(), mxREAL);
	std::copy(v.begin(), v.end(), mxGetPr(a));
	return a;
}

mxArray* getMexArray(const MatrixXd& A) {
	mxArray* B = mxCreateDoubleMatrix(A.rows(), A.cols(), mxREAL);
	memcpy(mxGetPr(B), A.data(), A.rows()*A.cols()*mxGetElementSize(B));
	return B;
}


// Interface between Matlab and C++
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

	// Validates input data
	if (nrhs < 4) {
		mexErrMsgIdAndTxt("CGGMfast:arguments", "Not enough input arguments.");
	}
	if (nrhs > 5) {
		mexErrMsgIdAndTxt("CGGMfast:arguments", "Too many input arguments.");
	}
	long n = mxGetM(prhs[INPUT_Y]);
	long p = mxGetN(prhs[INPUT_X]);
	long q = mxGetN(prhs[INPUT_Y]);
	if (mxGetM(prhs[INPUT_X]) != n) {
		mexErrMsgIdAndTxt("CGGMfast:arguments", "X and Y must have same n.");
	}

	// Inputs options
	double lambda_y = mxGetScalar(prhs[INPUT_LAMBDA_Y]);
	double lambda_x = mxGetScalar(prhs[INPUT_LAMBDA_X]);
	CGGMOptions options;
	if (nrhs > INPUT_OPTIONS) {
		for (int i = 0; i < mxGetNumberOfFields(prhs[INPUT_OPTIONS]); i++) {
			const char* name = mxGetFieldNameByNumber(prhs[INPUT_OPTIONS], i);
			mxArray* value = mxGetFieldByNumber(prhs[INPUT_OPTIONS], 0, i);
			if (!strcmp(name, "quiet")) {
				options.quiet = mxGetScalar(value);
			} else if (!strcmp(name, "max_outer_iters")) {
				options.max_outer_iters = (int)mxGetScalar(value);
			} else if (!strcmp(name, "max_ls_iters")) {
				options.max_ls_iters = (int)mxGetScalar(value);
			} else if (!strcmp(name, "sigma")) {
				options.sigma = mxGetScalar(value);
			} else if (!strcmp(name, "alpha")) {
				options.alpha = mxGetScalar(value);
			} else if (!strcmp(name, "beta")) {
				options.beta = mxGetScalar(value);
			} else if (!strcmp(name, "tol")) {
				options.tol = mxGetScalar(value);
			} else if (!strcmp(name, "cd_tol")) {
				options.cd_tol = mxGetScalar(value);
			}
			// Ignores unknown fields
		}
	}
	
	// Inputs data matrices
	double scaling = 1.0/sqrt(n);
	MatrixXd Y = Map<MatrixXd>(mxGetPr(prhs[INPUT_Y]), n, q);
	MatrixXd X = Map<MatrixXd>(mxGetPr(prhs[INPUT_X]), n, p);
	VectorXd Y_mean = Y.colwise().mean();
	VectorXd X_mean = X.colwise().mean();
	Y.rowwise() -= Y_mean.transpose();
	X.rowwise() -= X_mean.transpose();
	Y *= scaling;
	X *= scaling;

	// Initializes parameters
	SpMatrixC Lambda(q, q);
	Lambda.reserve(VectorXi::Constant(q,1));
	for (long i = 0; i < q; ++i) {
		Lambda.insert(i,i) = 1.0/(0.01 + Y.col(i).dot(Y.col(i)));
	}
	SpMatrixC Theta(p, q);

	CGGMStats stats;
	CGGMfast(Y, X, lambda_y, lambda_x, options, Lambda, Theta, &stats);
	
	// Outputs parameters
	plhs[OUTPUT_LAMBDA] = mxCreateSparse(q,q,Lambda.nonZeros(),mxREAL);
	double* lambda_vals = mxGetPr(plhs[OUTPUT_LAMBDA]);
	mwIndex* lambda_irs = mxGetIr(plhs[OUTPUT_LAMBDA]);
	mwIndex* lambda_jcs = mxGetJc(plhs[OUTPUT_LAMBDA]);
	mwIndex k = 0;
	for (int j = 0; j < Lambda.outerSize(); ++j) { // loops through cols
		lambda_jcs[j] = k;
		for (SpMatrixC::InnerIterator it(Lambda,j); it; ++it) {
			lambda_vals[k] = it.value();
			lambda_irs[k] = (mwSize)it.row();
			k++;
		}
	}
	lambda_jcs[q] = k;

	plhs[OUTPUT_THETA] = mxCreateSparse(p,q,Theta.nonZeros(),mxREAL);
	double* theta_vals = mxGetPr(plhs[OUTPUT_THETA]);
	mwIndex* theta_irs = mxGetIr(plhs[OUTPUT_THETA]);
	mwIndex* theta_jcs = mxGetJc(plhs[OUTPUT_THETA]);
	k = 0;
	for (int j = 0; j < Theta.outerSize(); ++j) { // loops thru cols
		theta_jcs[j] = k;
		for (SpMatrixC::InnerIterator it(Theta,j); it; ++it) {
			theta_vals[k] = it.value();
			theta_irs[k] = (mwSize)it.row();
			k++;
		}
	}
	theta_jcs[q] = k;
		

	// Outputs history
	mxArray* history = mxCreateStructArray(
		2, kHistoryDims, kNumHistoryFields, kHistoryFieldNames);
	mxSetFieldByNumber(history, 0, HISTORY_OBJVAL, getMexArray(stats.objval));
	mxSetFieldByNumber(history, 0, HISTORY_TIME, getMexArray(stats.time));
	mxSetFieldByNumber(history, 0, HISTORY_ACTIVE_SET_SIZE, 
		getMexArray(stats.active_set_size));
	mxSetFieldByNumber(history, 0, HISTORY_SUBGRAD, 
		getMexArray(stats.subgrad));
	mxSetFieldByNumber(history, 0, HISTORY_L1NORM, getMexArray(stats.l1norm));

	plhs[OUTPUT_HISTORY] = history;
}

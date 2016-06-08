#if defined(LANG_M) || defined(MATLAB_MEX_FILE)
#include <mex.h>
#define PRINTF mexPrintf
#else
#include <stdio.h>
#define PRINTF printf
#endif

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <stdint.h>
#include <sys/time.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include "pseudo.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SimplicialLLT;
using Eigen::Success;
using Eigen::Upper;

using namespace std;

// Returns difference in seconds
double toddiff(struct timeval *start, struct timeval *end) {
	long long tstart = start->tv_sec * 1000000 + start->tv_usec;
	long long tend = end->tv_sec * 1000000 + end->tv_usec;
	return ((double)(tend - tstart))/1000000.0;
}

double SoftThreshold(double a, double kappa) {
	return fmax(0, a-kappa) - fmax(0, -a-kappa);
}

void print(const SpMatrixC& A) {
	PRINTF("p: %ld, nnz: %ld\n", A.outerSize(), A.nonZeros());
	for (long k = 0; k < A.outerSize(); ++k) {
		for (InIter it(A,k); it; ++it) {
			PRINTF("%ld %ld %.10lf\n", it.row()+1, it.col()+1, it.value());
		}
	}
}

double L1Norm(const SpMatrixC& A) {
	double result = 0;
	for (long k = 0; k < A.outerSize(); ++k) {
		for (InIter it(A,k); it; ++it) {
			result += fabs(it.value());
		}
	}
	return result;
}

double L1SubGrad(double x, double g, double lambda) {
	if (x > 0) {
		return g + lambda;
	} else if (x < 0) {
		return g - lambda;
	} else {
		return fmax(fabs(g) - lambda, 0);
	}
}

void lasso_cd(
		const VectorXd& y, 
		const MatrixXd& X, 
		vector<double>& lambdaReg,
		VectorXd& w, // initialize before calling
		double& obj,
		PseudoOptions& options) {
	long n = X.rows();
	long p = X.cols();
	
	vector<double> Xi2(p);
	for (long i = 0; i < p; i++) {
		Xi2[i] = X.col(i).dot(X.col(i));
	}

	vector<long> ix_list(p);
	for (long i = 0; i < p; i++) {
		ix_list[i] = i;
	}
	
	w.resize(p);
	for (long i = 0; i < p; i++) {
		w(i) = 0;
	}

	VectorXd resid = y - X * w; 
	for (int iter = 0; iter < options.max_iters; iter++) {
		random_shuffle(ix_list.begin(), ix_list.end());
		for (long ix = 0; ix < p; ix++) {
			long i = ix_list[ix];
			double lamb = lambdaReg[i];
			double old_wi = w(i);

			VectorXd vi = resid + old_wi*X.col(i);
			double deltai = X.col(i).dot(vi);
			if (deltai > lamb) {
				w(i) = (deltai - lamb) / Xi2[i];
			} else if (deltai < -lamb) {
				w(i) = (deltai + lamb) / Xi2[i];
			} else {
				w(i) = 0;
			}
			resid = resid + (old_wi - w(i))*X.col(i);
		}
	}
	resid = y - X * w;
	obj = 0.5*resid.dot(resid);
	for (long ix = 0; ix < p; ix++) {
		obj += fabs(lambdaReg[ix] * w[ix]);
	}
}

void Pseudo(
		const MatrixXd& Y, 
		const MatrixXd& X, 
		double lambda_y,
		double lambda_x,
		PseudoOptions& options,
		SpMatrixC& Lambda,
		SpMatrixC& Theta,
		PseudoStats* stats) {
	if (!options.quiet) {
		PRINTF("Eigen instruction sets: %s\n", Eigen::SimdInstructionSetsInUse());
	}

	srand(1);
	long p = Theta.rows();
	long q = Theta.cols();
	long n = Y.rows();
	struct timeval start_time, current_time;
	gettimeofday(&start_time, NULL);

	MatrixXd Yscaled = 1.0/sqrt(n) * Y;
	MatrixXd Xscaled = 1.0/sqrt(n) * X;

	vector<double> lambdaReg(p+q-1, lambda_x);
	for (long i = 0; i < q-1; i++) {
		lambdaReg[i] = lambda_y;
	}
	vector<Triplet> Betas_triplets;
	vector<double> Vars(q);
	MatrixXd preds(n, p+q-1);
	preds.rightCols(p) = Xscaled;
	for (long i = 0; i < q; i++) {
		VectorXd Beta(p+q-1); // TODO- SparseVector
		double obj = 0;

		// TODO- pull lasso_cd in to avoid moving data around
		// TODO- compute Xi2 only once
		// TODO- use SAFE rule
		preds.leftCols(i) = Yscaled.leftCols(i);
		for (long j = i + 1; j < q; j++) {
			preds.col(j-1) = Yscaled.col(j);
		}
		VectorXd outs = Yscaled.col(i);
		lasso_cd(outs, preds, lambdaReg, Beta, obj, options);
		
		VectorXd hatYi = preds * Beta;
		VectorXd errs = Y.col(i) - hatYi;
		double errs_mean = errs.mean();
		VectorXd centered_errs = (errs.array() - errs_mean).matrix();
		double Var = (1.0/n) * centered_errs.dot(centered_errs);
		for (long ix = 0; ix < p+q-1; ix++) {
			if (Beta(ix) != 0) {
				Betas_triplets.push_back(Triplet(ix,i, Beta(ix)));
			}
		}
		Vars[i] = Var;
	}
	SpMatrixC Betas_sp(p+q-1, q);
	Betas_sp.setFromTriplets(Betas_triplets.begin(), Betas_triplets.end());

	vector<Triplet> Omega_triplets;
	for (long i = 0; i < q; i++) {
		Omega_triplets.push_back(Triplet(i, i, 1.0 / Vars[i]));
		for (InIter it(Betas_sp, i); it.index() < q-1; ++it) {
			long j = it.index();
			long real_j = j;
			if (j >= i) {
				real_j = j + 1;
			}
			double Betas_ji = it.value();
			double Omega_ij = -1 * Vars[i] * Betas_ji;
			Omega_triplets.push_back(Triplet(i,real_j,Omega_ij));
		}
	}
	SpMatrixC Omega_sp(q, q);
	Omega_sp.setFromTriplets(Omega_triplets.begin(), Omega_triplets.end());

	SpMatrixC Omega_sp_T = Omega_sp.transpose();
	Omega_sp = 0.5*(Omega_sp + Omega_sp_T);
	vector<double> Omega_iis(q);
	if (options.diag_dominant) {
		for (long i = 0; i < q; i++) {
			double dom = 1e-2;
			double Omega_ii = 0;
			for (InIter it(Omega_sp,i); it; ++it) {
				if (it.index() == i) {
					Omega_ii = it.value();
				} else {
					dom += it.value();
				}
			}
			double new_Omega_ii = fmax(Omega_ii, dom);
			Omega_sp.coeffRef(i,i) = new_Omega_ii;
			Omega_iis[i] = new_Omega_ii;
		}
	} else {
		for (long i = 0; i < q; i++) {
			Omega_iis[i] = Omega_sp.coeffRef(i,i);
		}
	}
	Lambda = Omega_sp;

	vector<Triplet> Theta_triplets;
	for (long i = 0; i < q; i++) {
		for (InIter it(Betas_sp, i); it; ++it) {
			long j = it.index();
			if (j < q-1) {
				continue;
			}
			long real_j = j - (q - 1);
			double Theta_ji = -1 * Omega_iis[i] * it.value();
			Theta_triplets.push_back(Triplet(real_j, i, Theta_ji));
		}
	}
	Theta.setFromTriplets(Theta_triplets.begin(), Theta_triplets.end());
	
	gettimeofday(&current_time, NULL);
	stats->time.push_back(toddiff(&start_time, &current_time));
}


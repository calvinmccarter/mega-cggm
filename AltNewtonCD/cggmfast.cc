#if defined(LANG_M) || defined(MATLAB_MEX_FILE)
#include <mex.h>
#define PRINTF mexPrintf
#else
#include <stdio.h>
#define PRINTF printf
#endif

#include <cstdlib>
#include <iostream>
#include <vector>
#include <stdint.h>
#include <sys/time.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include "cggmfast.h"

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

struct TimeThetaCD {
        TimeThetaCD() :
    cd(0), qr(0) {}
    double cd;
    double qr;
};

struct TimeLambdaCD {
        TimeLambdaCD() :
    cd(0), ls(0) {}
    double cd;
    double ls; // includes update to R
};

struct LambdaState {
	LambdaState(
		const double logdetLambda_,
		const double trRtQ_)
		: logdetLambda(logdetLambda_),
	      trRtQ(trRtQ_) {}
	double logdetLambda;
	double trRtQ; // tr(Sigma*Theta'*Sxx*Theta)
};

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

// Works with upper diag
double L1NormUpperDiag(const SpMatrixC& A) {
	double result = 0;
	for (long k = 0; k < A.outerSize(); ++k) {
		for (InIter it(A,k); it; ++it) {
			if (it.row() < it.col()) {
				result += 2*fabs(it.value());
			} else if (it.row() == it.col()) {
				result += fabs(it.value());
			} else {
				continue;
			}
		}
	}
	return result;
}

// Works with either upper diag or symmetric A
double L1NormOffDiag(const SpMatrixC& A) {
	double result = 0;
	for (long k = 0; k < A.outerSize(); ++k) {
		for (InIter it(A,k); it; ++it) {
			if (it.row() < it.col()) {
				result += 2*fabs(it.value());
			} else {
				continue;
			}
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

// Interacts with Lambda as symmetric, not upper diag
double logdet(const SpMatrixC& Lambda) {
	SimplicialLLT<SpMatrixC> cholesky(Lambda);
	const SpMatrixC matrixL = cholesky.matrixL();
	return 2*matrixL.diagonal().array().log().sum();
}

double logdet(const SpMatrixC& matrixL, bool dummy) {
	return 2*matrixL.diagonal().array().log().sum();
}

// Computes trace(SA) = trace(X'*Y*A), where A is sparse matrix
double traceProduct(
		const SpMatrixC& A, const MatrixXd& X, const MatrixXd& Y) {
	long p = A.rows();
	long q = A.cols();
	long n = X.rows();
	// Complexity: nnz(A)*n
	
	double result = 0;
	for (long k = 0; k < A.outerSize(); ++k) {
		for (InIter it(A,k); it; ++it) {
			long i = it.row();
			long j = it.col();
			result += it.value() * X.col(i).dot(Y.col(j));
		}
	}
	return result;
}

// Computes trace(A' * B) = sum_ij(A_ij, B_ij)
double traceProduct(const MatrixXd& A, const MatrixXd& B) {
	return (A.cwiseProduct(B)).sum();
	/*
	long m = A.rows();
	long n = A.cols();
	// Complexity: mn
	
	double result = 0;
	for (long j = 0; j < n; ++j) { // assuming ColMajor order
		for (long i = 0; i < m; ++i) {
			result += A(i,j) * B(i,j);
		}
	}
	return result;
	*/
}

// Computes trace(A' * B)
double traceProduct(const SpMatrixC& A, const MatrixXd& B) {
	long m = A.rows();
	long n = A.cols();
	
	double result = 0;
	for (long k = 0; k < n; ++k) {
		for (InIter it(A,k); it; ++it) {
			long i = it.row();
			long j = it.col();
			result += it.value() * B(i, j);
		}
	}
	return result;
}

// Objective function to minimize
// Interacts with Lambda as symmetric, not upper diag
// Fast version
double Objective(
		const MatrixXd& Syy,
		const MatrixXd& Sxy,
		const MatrixXd& Sxx,
		double lambda_y,
		double lambda_x,
		const SpMatrixC& Lambda,
		const SpMatrixC& Theta,
		const MatrixXd& Q,
		const MatrixXd& R) {
	// f = -logdet(Lambda) + tr(Syy*Lambda + 2Sxy'*Theta + R'*Q)
	return -logdet(Lambda) + traceProduct(Lambda, Syy) +
		2*traceProduct(Theta, Sxy) + traceProduct(R, Q) +
		lambda_x*L1Norm(Theta) + lambda_y*L1NormOffDiag(Lambda);
}	

// Returns subgradient
// Adds 0s to Theta at locations of free set
double ThetaActiveSet(
		const MatrixXd& Y,
		const MatrixXd& X,
		SpMatrixC& Theta,
		const MatrixXd& R,
		double lambda,
		const CGGMOptions& options) {
	long n_x = X.rows();
	long n_y = Y.rows();
	long n_o = min(n_x, n_y);
	long p = Theta.rows();
	long q = Theta.cols();
	double subgrad = 0.0;
	MatrixXd G;
	G.noalias() = (2.0 / n_o)*X.topRows(n_o).transpose()*Y.topRows(n_o) + 
		(2.0 / sqrt(n_x)) * X.transpose() * R;

	vector<Triplet> triplets;
	triplets.reserve(Theta.nonZeros());
	if (options.refit) {
		for (long j = 0; j < q; j++) {
			InIter it(Theta,j);
			for (long i = 0; i < p; i++) {
				double Theta_ij = 0;
				if (it && it.row() == i) {
					Theta_ij = it.value();
					++it;
				}
				if (Theta_ij != 0) {
					triplets.push_back(Triplet(i, j, Theta_ij));
					subgrad += fabs(L1SubGrad(Theta_ij, G(i,j), lambda));
				}
			}
		}
		Theta.setFromTriplets(triplets.begin(), triplets.end());
		return subgrad;
	} else {
		for (long j = 0; j < q; j++) {
			InIter it(Theta,j);
			for (long i = 0; i < p; i++) {
				double Theta_ij = 0;
				if (it && it.row() == i) {
					Theta_ij = it.value();
					++it;
				}
				if (Theta_ij != 0 || fabs(G(i,j)) > lambda) {
					triplets.push_back(Triplet(i, j, Theta_ij));
					subgrad += fabs(L1SubGrad(Theta_ij, G(i,j), lambda));
				}
			}
		}
		Theta.setFromTriplets(triplets.begin(), triplets.end());
		return subgrad;
	}
}

// Returns subgradient
// Interacts with Lambda as upper diag (i <= j)
// Adds 0s to Lambda at locations of free set
double LambdaActiveSet(
		const MatrixXd& Syy,
		const MatrixXd& Sxy,
		const MatrixXd& Sxx,
		SpMatrixC& Lambda,
		const MatrixXd& Sigma,
		const MatrixXd& Psi,
		double lambda,
		const CGGMOptions& options) {
	long q = Lambda.cols();
	double subgrad = 0.0;
	
	vector<Triplet> triplets;
	triplets.reserve(Lambda.nonZeros());
	MatrixXd G = Syy - Sigma - Psi;
	if (options.refit) {	
		for (long j = 0; j < q; j++) {
			InIter it(Lambda, j);

			for (long i = 0; i < j; i++) {
				double Lambda_ij = 0;
				if (it && it.row() == i) {
					Lambda_ij = it.value();
					++it;
				}
				if (Lambda_ij != 0) {
					triplets.push_back(Triplet(i, j, Lambda_ij));
					subgrad += 2*fabs(L1SubGrad(Lambda_ij, G(i,j), lambda));
				}
			}
			triplets.push_back(Triplet(j, j, it.value()));
			subgrad += fabs(G(j,j));
		}
		Lambda.setFromTriplets(triplets.begin(), triplets.end());
		return subgrad;
	} else {
		for (long j = 0; j < q; j++) {
			InIter it(Lambda, j);

			for (long i = 0; i < j; i++) {
				double Lambda_ij = 0;
				if (it && it.row() == i) {
					Lambda_ij = it.value();
					++it;
				}
				if (Lambda_ij != 0 || fabs(G(i,j)) > lambda) {
					triplets.push_back(Triplet(i, j, Lambda_ij));
					subgrad += 2*fabs(L1SubGrad(Lambda_ij, G(i,j), lambda));
				}
			}
			triplets.push_back(Triplet(j, j, it.value()));
			subgrad += fabs(G(j,j));
		}
		Lambda.setFromTriplets(triplets.begin(), triplets.end());
		return subgrad;
	}
}

// Coordinate descent for Theta
// Updates Theta, Q, and R
// Fast version
void ThetaCoordinateDescent(
		const MatrixXd& Syy, // debugging
		const MatrixXd& Sxy, 
		const MatrixXd& Sxx,
		const MatrixXd& X,
		SpMatrixC& Theta,
		const MatrixXd& Sigma,
		MatrixXd& Q,
		MatrixXd& R,
		double lambda,
		const CGGMOptions& options,
		TimeThetaCD& time_report_theta) {
	struct timeval start_time, end_time;
	gettimeofday(&start_time, NULL);

	long n_x = X.rows();
	long q = Theta.cols();
	MatrixXd V;
	V.noalias() = Theta * Sigma;
	vector<long> columns;
	columns.reserve(q);
	for (long j = 0; j < q; ++j) {
		columns.push_back(j);
	}

	for (int iter = 0; iter < options.max_Theta_iters; ++iter) {
		// Shuffles columns
		for (long ix1 = 0; ix1 < q; ++ix1) {
			long ix2 = ix1 + rand() % (q - ix1);
			long tmp = columns[ix1];
			columns[ix1] = columns[ix2];
			columns[ix2] = tmp;
		}

		for (long j_ix = 0; j_ix < q; ++j_ix) {
			long j = columns[j_ix];
			long nnz_j = Theta.col(j).nonZeros();

			// Shuffles within columns
			vector<InIter*> itThetas;
			itThetas.reserve(nnz_j);
			InIter itTheta(Theta, j);
			for (; itTheta; ++itTheta) {
				itThetas.push_back(new InIter(itTheta));
			}
			for (long ix1 = 0; ix1 < nnz_j; ++ix1) {
				long ix2 = ix1 + rand() % (nnz_j - ix1);
				InIter* tmp = itThetas[ix1];
				itThetas[ix1] = itThetas[ix2];
				itThetas[ix2] = tmp;
			}

			// Performs coordinate descent
			for (long ix = 0; ix < nnz_j; ++ix) {
				InIter* itPtr = itThetas[ix];
				long i = itPtr->row();
				double Theta_ij = itPtr->value();

				double a = 2*Sigma(j,j)*Sxx(i,i);
				double b = 2*Sxy(i,j) + 2*Sxx.col(i).dot(V.col(j));
				double c = Theta_ij;
				double mu = -c + SoftThreshold(c - b/a, lambda/a);
				itPtr->valueRef() += mu;
				V.row(i) += mu*Sigma.row(j);
			}

			for (long ix = 0; ix < nnz_j; ++ix) {
				delete itThetas.back();
				itThetas.pop_back();
			}
		}
	}
	gettimeofday(&end_time, NULL);
	time_report_theta.cd += toddiff(&start_time, &end_time);

	gettimeofday(&start_time, NULL);
	Q.noalias() = (1.0 / sqrt(n_x)) * X * Theta;
	R.noalias() = Q * Sigma;
	gettimeofday(&end_time, NULL);
	time_report_theta.qr += toddiff(&start_time, &end_time);
}

// Coordinate descent for Lambda
// Interacts with Lambda as upper diag (i <= j)
// Updates Lambda, Sigma, R, and lState
// Fast version
void LambdaCoordinateDescent(
		const MatrixXd& Syy,
		const MatrixXd& Sxy, // debugging
		const MatrixXd& Sxx, // debugging
		SpMatrixC& Lambda,
		SpMatrixC& Theta,  // debugging
		MatrixXd& Sigma,
		const MatrixXd& Q,
		MatrixXd& R,
		const MatrixXd& Psi,
		double lambda,
		LambdaState& lState,
		const CGGMOptions& options,
		TimeLambdaCD& time_report_lambda) {
	
	struct timeval start_time, end_time;
	gettimeofday(&start_time, NULL);
	
	// Sets up constants
	long q = Lambda.rows();
	MatrixXd Iq = MatrixXd::Identity(q, q);
	
	// Sets up updated matrices
	SpMatrixC Delta = Lambda; // all 0s
	for (long k = 0; k < Delta.outerSize(); ++k) {
		for (InIter it(Delta,k); it; ++it) {
			it.valueRef() = 0.0;
		}
	}
	MatrixXd U;
	U.noalias() = Delta * Sigma;

	// Columns to shuffle
	vector<long> columns;
	columns.reserve(q);
	for (long j = 0; j < q; ++j) {
		columns.push_back(j);
	}

	for (int iter = 0; iter < options.max_Lambda_iters; ++iter) {
		// Shuffles columns
		for (long ix1 = 0; ix1 < q; ++ix1) {
			long ix2 = ix1 + rand() % (q - ix1);
			long tmp = columns[ix1];
			columns[ix1] = columns[ix2];
			columns[ix2] = tmp;
		}
		for (long j_ix = 0; j_ix < q; ++j_ix) {
			long j = columns[j_ix];
			long nnz_j = Lambda.col(j).nonZeros();

			// Shuffles within columns
			vector<InIter*> itLambdas;
			vector<InIter*> itDeltas;
			itLambdas.reserve(nnz_j);
			itDeltas.reserve(nnz_j);
			InIter itLambda(Lambda, j);
			InIter itDelta(Delta, j);
			for (; itLambda && itDelta; ++itLambda, ++itDelta) {
				itLambdas.push_back(new InIter(itLambda));
				itDeltas.push_back(new InIter(itDelta));
			}
			for (long ix1 = 0; ix1 < nnz_j; ++ix1) {
				long ix2 = ix1 + rand() % (nnz_j - ix1);
				InIter* tmp = itLambdas[ix1];
				itLambdas[ix1] = itLambdas[ix2];
				itLambdas[ix2] = tmp;
				tmp = itDeltas[ix1];
				itDeltas[ix1] = itDeltas[ix2];
				itDeltas[ix2] = tmp;
			}

			// Performs coordinate descent
			for (long ix = 0; ix < nnz_j; ++ix) {
				InIter* itLambdaPtr = itLambdas[ix];
				InIter* itDeltaPtr = itDeltas[ix];
				long i = itLambdaPtr->row();
				double Lambda_ij = itLambdaPtr->value();
				double Delta_ij = itDeltaPtr->value();
				
				if (i == j) { // diagonal elements
					double a = Sigma(i,i)*Sigma(i,i) + 2*Sigma(i,i)*Psi(i,i);
					double b = Syy(i,i) - Sigma(i,i) - Psi(i,i) 
						+ Sigma.row(i)*U.col(i) + 2*Psi.row(i)*U.col(i);
					double mu = -b/a;
					itDeltaPtr->valueRef() += mu;
					U.row(i) += mu*Sigma.row(i);
				} else { // off-diagonal elements
					double a = Sigma(i,j)*Sigma(i,j) + Sigma(i,i)*Sigma(j,j) 
						+ Sigma(i,i)*Psi(j,j) + 2*Sigma(i,j)*Psi(i,j) 
						+ Sigma(j,j)*Psi(i,i);
					double b = Syy(i,j) - Sigma(i,j) + Sigma.row(i)*U.col(j)
						- Psi(i,j) + Psi.row(i)*U.col(j) + Psi.row(j)*U.col(i);
					double c = Lambda_ij + Delta_ij;  
					double mu = -c + SoftThreshold(c - b/a, lambda/a);
					itDeltaPtr->valueRef() += mu;
					U.row(i) += mu*Sigma.row(j);
					U.row(j) += mu*Sigma.row(i);
				}
			}

			for (long ix = 0; ix < nnz_j; ++ix) {
				delete itLambdas.back();
				delete itDeltas.back();
				itLambdas.pop_back();
				itDeltas.pop_back();
			}
		}
	}
	gettimeofday(&end_time, NULL);
	time_report_lambda.cd += toddiff(&start_time, &end_time);
	
	// Backtracking line search
	gettimeofday(&start_time, NULL);
	double alpha = options.alpha;
	bool success = false;
	SpMatrixC Lambda_sym;
	Lambda_sym = Lambda.selfadjointView<Upper>();
	SpMatrixC Delta_sym;
	Delta_sym = Delta.selfadjointView<Upper>();
	SpMatrixC Lambda_alpha; // symmetric, not upper diag
	double logdetLambda_alpha, trRtQ_alpha;

	double trGradDelta = traceProduct(Delta, Syy) 
		- traceProduct(Delta, Sigma) - traceProduct(Delta, Psi);
	SpMatrixC LambdaPlusDelta = Lambda_sym + Delta_sym;
	double RHS = options.sigma*(trGradDelta 
		+ lambda*(L1NormOffDiag(LambdaPlusDelta) - L1NormOffDiag(Lambda)));

	for (int lsiter = 0; lsiter < options.max_ls_iters; ++lsiter) {
		Lambda_alpha = Lambda_sym + alpha*Delta_sym;
		SimplicialLLT<SpMatrixC> cholesky(Lambda_alpha);
		if (cholesky.info() != Success) {
			if (!options.quiet) {
				PRINTF("   line search %d, alpha=%f, not PD\n", lsiter, alpha);
			}
			alpha *= options.beta;
			continue;
		}
		logdetLambda_alpha = logdet(cholesky.matrixL(), 0);
		MatrixXd R_alpha = (cholesky.solve(Q.transpose())).transpose();
		trRtQ_alpha = traceProduct(R_alpha, Q);
		double LHS = -logdetLambda_alpha + lState.logdetLambda 
			+ alpha*traceProduct(Delta_sym, Syy) + trRtQ_alpha - lState.trRtQ
			+ lambda*(L1NormOffDiag(Lambda_alpha) - L1NormOffDiag(Lambda_sym));

		if (LHS <= alpha*RHS) {
			success = true;
			Lambda = Lambda + alpha*Delta;
			Sigma = cholesky.solve(Iq);
			R = R_alpha;
			lState = LambdaState(logdetLambda_alpha, trRtQ_alpha);
            /*
			if (!options.quiet) {
				PRINTF(
                    "   line search %d, alpha=%f, sufficient decrease=%f\n",
				    lsiter, alpha, LHS);
			}
            */
			break;
		} else if (!options.quiet) {
			PRINTF("   line search %d, alpha=%f, insufficient decrease=%f\n",
				lsiter, alpha, LHS);
		}
		alpha *= options.beta;
	}
	if (success) {
	} else if (!options.quiet) {
		PRINTF("  Lambda line search failed\n");
	}
	gettimeofday(&end_time, NULL);
	time_report_lambda.ls += toddiff(&start_time, &end_time);
}

void CGGMfast(
		const MatrixXd& Y, 
		const MatrixXd& X, 
		double lambda_y,
		double lambda_x,
		CGGMOptions& options,
		SpMatrixC& Lambda,
		SpMatrixC& Theta,
		CGGMStats* stats) {
	if (!options.quiet) {
		PRINTF("Eigen instruction sets: %s\n", Eigen::SimdInstructionSetsInUse());
	}
#if defined(LANG_M) || defined(MATLAB_MEX_FILE)
  	mexEvalString("drawnow;");  
#endif

	srand(1);
	long p = Theta.rows();
	long q = Theta.cols();
	long n_y = Y.rows();
	long n_x = X.rows();
	long n_o = min(n_y, n_x); 
	struct timeval start_time, current_time;
	gettimeofday(&start_time, NULL);

	MatrixXd Syy, Sxy, Sxx;
	Syy.noalias() = (1.0 / n_y) * (Y.transpose()*Y);
	Sxy.noalias() = (1.0 / n_o) * (X.topRows(n_o).transpose()*Y.topRows(n_o));
	Sxx.noalias() = (1.0 / n_x) * (X.transpose()*X);

	SimplicialLLT<SpMatrixC> cholesky(Lambda);
	if (cholesky.info() != Success) { 
		if (!options.quiet) PRINTF("Lambda0 not positive definite\n");
		return;
	}
	MatrixXd Sigma, Q, R;
	Sigma = cholesky.solve(MatrixXd::Identity(q,q));
	Q.noalias() = (1.0 / sqrt(n_x)) * X * Theta;
	R.noalias() = Q * Sigma;
	LambdaState lState(logdet(Lambda), traceProduct(Q, R));
	struct timeval mini_start_time, mini_end_time;
	
	for (int t_outer = 0; t_outer < options.max_outer_iters; ++t_outer) {
		MatrixXd Psi;
		// Calculate the subgradient and free sets
		gettimeofday(&mini_start_time, NULL);
    	double subgrad_theta = ThetaActiveSet(
			Y, X, Theta, R, lambda_x, options);
		gettimeofday(&mini_end_time, NULL);
		double theta_active_time = toddiff(&mini_start_time, &mini_end_time);

		gettimeofday(&mini_start_time, NULL);
		Psi.noalias() = R.transpose() * R;
		double subgrad_lambda = LambdaActiveSet(
			Syy, Sxy, Sxx, Lambda, Sigma, Psi, lambda_y, options);
		gettimeofday(&mini_end_time, NULL);
		double lambda_active_time = toddiff(&mini_start_time, &mini_end_time);

    	double l1_norm_theta = L1Norm(Theta);
		double l1_norm_lambda = L1NormUpperDiag(Lambda);
		
		double subgrad = subgrad_theta + subgrad_lambda;
		double l1_norm = l1_norm_theta + l1_norm_lambda;
		if (t_outer > 0) {
			stats->subgrad.push_back(subgrad);
			stats->l1norm.push_back(l1_norm);
		}

    	if (subgrad < options.tol*l1_norm) {
			if (!options.quiet) {
				PRINTF("Converged, subgrad=%f, norm=%f\n", subgrad, l1_norm);
			}
			break;
		}

		// CD for Lambda, also updates Sigma, R, and lState
		struct TimeLambdaCD time_report_lambda;
		gettimeofday(&mini_start_time, NULL);
		LambdaCoordinateDescent(
			Syy, Sxy, Sxx, Lambda, Theta, 
			Sigma, Q, R, Psi, lambda_y, lState, 
			options, time_report_lambda);
		gettimeofday(&mini_end_time, NULL);
		double lambda_cd_time =	toddiff(&mini_start_time, &mini_end_time);

		// CD for Theta, also updates Q and R
		struct TimeThetaCD time_report_theta;
		gettimeofday(&mini_start_time, NULL);
		ThetaCoordinateDescent(
			Syy, Sxy, Sxx, X, Theta, Sigma, Q, R, lambda_x, 
			options, time_report_theta);
		lState.trRtQ = traceProduct(R, Q);
		gettimeofday(&mini_end_time, NULL);
		double theta_cd_time = toddiff(&mini_start_time, &mini_end_time);


		if (!options.quiet) {
			PRINTF("Iteration %d, Lambda(subgrad=%f,active=%ld,norm=%f\n",
	 			t_outer, subgrad_lambda, Lambda.nonZeros(), l1_norm_lambda);
			PRINTF("             Theta(subgrad=%f,active=%ld,norm=%f\n",
				subgrad_theta, Theta.nonZeros(), l1_norm_theta);
			fflush(stdout);
		}

		SpMatrixC Lambda_sym;
		Lambda_sym = Lambda.selfadjointView<Upper>();
		double f = Objective(
			Syy, Sxy, Sxx, lambda_y, lambda_x, 
			Lambda_sym, Theta, Q, R);	

		stats->objval.push_back(f);
		gettimeofday(&current_time, NULL);
		stats->time.push_back(toddiff(&start_time, &current_time));
		long active_set_size = Theta.nonZeros() + Lambda.nonZeros();
		stats->active_set_size.push_back((double)active_set_size);
		stats->active_theta.push_back((double)Theta.nonZeros());
		stats->active_lambda.push_back((double)Lambda.nonZeros());

		// Output iteration time breakdown
		stats->time_lambda_active.push_back(lambda_active_time);
		stats->time_theta_active.push_back(theta_active_time);
		stats->time_lambda_cd.push_back(lambda_cd_time);
		stats->time_lambda_cd_cd.push_back(time_report_lambda.cd);
		stats->time_lambda_cd_ls.push_back(time_report_lambda.ls);
		stats->time_theta_cd.push_back(theta_cd_time);
		stats->time_theta_cd_cd.push_back(time_report_theta.cd);
		stats->time_theta_cd_qr.push_back(time_report_theta.qr);
	} // outer Newton iteration 

	// Return symmetric matrix
	SpMatrixC Lambda_sym;
	Lambda_sym = Lambda.selfadjointView<Upper>();
	Lambda = Lambda_sym;
}


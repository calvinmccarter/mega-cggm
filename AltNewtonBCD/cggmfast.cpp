#include "cggmfast.h"
#include "util.h"

#if defined(LANG_M) || defined(MATLAB_MEX_FILE)
#include <mex.h>
#define PRINTF mexPrintf
#else
#include <stdio.h>
#define PRINTF printf
#endif

#include <algorithm>
#include <assert.h>
#include <exception>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <vector>

using namespace std;

struct TimeThetaCD {
        TimeThetaCD() :
    prep(0), inv(0), s(0), vp(0), cd(0), qr(0) {}

    double prep;
	double inv;
	double s;
	double vp;
    double cd;
    double qr;
};

struct TimeLambdaCD {
        TimeLambdaCD() :
    prep(0), prep_z(0), prep_q(0), prep_d(0), cd(0), cd_prep(0), ls(0) {}
    double prep;
	double prep_z;
	double prep_q;
	double prep_d;
    double cd;
	double cd_prep;
    double ls; // includes update to R
};

struct LambdaState {
	LambdaState(
		const double logdetLambda_,
		const double trRtQ_,
		const int numNoBacktracking_)
		: logdetLambda(logdetLambda_),
	      trRtQ(trRtQ_),
	      numNoBacktracking(numNoBacktracking_) {}
	double logdetLambda;
	double trRtQ; // tr(Sigma*Theta'*Sxx*Theta)
	int numNoBacktracking;
};

double SoftThreshold(double a, double kappa) {
	return max(0.0, a-kappa) - max(0.0, -a-kappa);
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

// returns trace(Delta * inv(Lambda))
double traceDeltaSigma(
		const smat_t &Delta, const smat_t &Lambda, double tol) {
	double result = 0;
	if (Delta.is_symmetric != 1) {
		#pragma omp parallel for reduction(+:result) schedule(dynamic,64)
		for (long i = 0; i < Delta.p; i++) {
			vector<double> Sigma_col(Delta.p, 0);
			Lambda.ComputeInv(Sigma_col, i, tol);
			for (long idx = Delta.row_ptr[i];
					idx < Delta.row_ptr[i+1]; idx++) {
				long j = Delta.col_idx[idx];
				if (i == j) {
					result += Delta.values[idx] * Sigma_col[j];
				} else {
					result += 2 * Delta.values[idx] * Sigma_col[j];
				}
			}
		}
	} else {
		#pragma omp parallel for reduction(+:result) schedule(dynamic,64)
		for (long i = 0; i < Delta.p; i++) {
			vector<double> Sigma_col(Delta.p, 0);
			Lambda.ComputeInv(Sigma_col, i, tol);
			for (long idx = Delta.row_ptr[i];
					idx < Delta.row_ptr[i+1]; idx++) {
				long j = Delta.col_idx[idx];
				result += Delta.values[idx] * Sigma_col[j];
			}
		}
	}
	return result;
}

// returns trace(inv(Lambda_alpha) * (Q^T Q)), Q is column-major
double traceRtQ_alpha(
		const smat_t &Lambda_alpha,
		const vector< vector<double> > &Q,
		double tol) {
	long q = Lambda_alpha.p;
	long n = Q[0].size();

	// Compute R_alpha
	vector< vector<double> > R_alpha(q, vector<double>(n, 0));
	#pragma omp parallel for schedule(static)
	for (long i = 0; i < n; i++) {
		vector<double> Q_i(q, 0);
		for (long j = 0; j < q; j++) {
			Q_i[j] = Q[j][i];
		}
		vector<double> R_i(q, 0);
		Lambda_alpha.ComputeAinvb(Q_i, R_i, tol);
		for (long j = 0; j < q; j++) {
			R_alpha[j][i] = R_i[j];
		}
	}

	double result = 0;
	#pragma omp parallel for reduction(+:result) schedule(static)
	for (long i = 0; i < q; i++) {
		for (long nn = 0; nn < n; nn++) {
			result += R_alpha[i][nn] * Q[i][nn];
		}
	}
	return result;
}

double Objective(
		const vector< vector<double> > &Y,
		const vector< vector<double> > &X,
		double lambda_y,
		double lambda_x,
		const smat_t &Lambda,
		const sparse_t &Theta,
		const vector< vector<double> > &Q,
		const vector< vector<double> > &R,
		double tol) {
	long n_x = X[0].size();
	long n_y = Y[0].size();
	long n_o = min(n_x, n_y);
	double logdet;
	Lambda.ComputeLogdet(logdet, tol);
	//PRINTF("logdet=%.10f, trLambdaYtY=%.10f, 2trThetaXtY=%.10f, trRtQ=%.10f \n",
	//	logdet, Lambda.traceProduct(Y, Y), 2*Theta.traceProduct(X,Y),
	//	traceProduct(R, Q));	
	return -1*logdet + (1.0/n_y)*Lambda.traceProduct(Y, Y) +
		2*(1.0/n_o)*Theta.traceProduct(X, Y) + traceProduct(R, Q) +
		lambda_x*Theta.L1Norm() + lambda_y*Lambda.L1NormOffDiag();
}

double LogLik(
		const vector< vector<double> > &Y,
		const vector< vector<double> > &X,
		const smat_t &Lambda,
		const sparse_t &Theta,
		const vector< vector<double> > &Q,
		const vector< vector<double> > &R,
		double tol) {
	long n_x = X[0].size();
	long n_y = Y[0].size();
	long n_o = min(n_x, n_y);
	double logdet;
	Lambda.ComputeLogdet(logdet, tol);
	//PRINTF("logdet=%.10f, trLambdaYtY=%.10f, 2trThetaXtY=%.10f, trRtQ=%.10f \n",
	//	logdet, Lambda.traceProduct(Y, Y), 2*Theta.traceProduct(X,Y),
	//	traceProduct(R, Q));	
	double mean_ll = logdet - (1.0/n_y)*Lambda.traceProduct(Y, Y) -
		2*(1.0/n_o)*Theta.traceProduct(X, Y) - traceProduct(R, Q);
	return n_o*mean_ll;
}

// Assumes non-symmetrized Lambda
long CountNonzeros(const smat_t &Lambda, const sparse_t &Theta) {
	long Lambda_nnz = 0;
	for (long idx = 0; idx < Lambda.nnz; idx++) {
		if (Lambda.values[idx] != 0) {
			Lambda_nnz++;
		}
	}
	long Theta_nnz = 0;
	for (long idx = 0; idx < Theta.nnz; idx++) {
		if (Theta.values[idx] != 0) {
			Theta_nnz++;
		}
	}
	return Lambda_nnz + Theta_nnz;
}

double ThetaActiveSet(
		const vector< vector<double> > &Y,
		const vector< vector<double> > &X,
		double lambda_x,
		sparse_t &Theta,
		const vector< vector<double> > &R,
		const CGGMOptions &options) {
	double subgrad = 0;
	vector<Triplet> triplets;
	triplets.reserve(Theta.nnz);
	long n_x = X[0].size();
	long n_y = Y[0].size();
	long n_o = min(n_x, n_y);
	long p = X.size();
	long q = Theta.q;
	
	#pragma omp parallel for \
		shared(triplets) reduction(+:subgrad) schedule(static)
	for (long i = 0; i < p; i++) {
		long idx = Theta.row_ptr[i];
		for (long j = 0; j < q; j++) {
			double XtY_ij = 0;
			for (long nn = 0; nn < n_o; nn++) {
				XtY_ij += X[i][nn] * Y[j][nn];
			}
			double XtR_ij = 0;
			for (long nn = 0; nn < n_x; nn++) {
				XtR_ij += X[i][nn] * R[j][nn];
			}
			double G_ij = (2.0 / n_o) * XtY_ij + (2.0 / sqrt(n_x))*XtR_ij;
			double Theta_ij = 0;
			if (idx < Theta.row_ptr[i+1] && Theta.col_idx[idx] == j) {
				Theta_ij = Theta.values[idx];
				idx++;
			}
			if (Theta_ij != 0 || (!options.refit && fabs(G_ij) > lambda_x) ) {
				subgrad += fabs(L1SubGrad(Theta_ij, G_ij, lambda_x));
				Triplet triplet(i, j, Theta_ij);
				#pragma omp critical
				triplets.push_back(triplet);
			}
		}
	}

	Theta.setTriplets(p, q, triplets);
	return subgrad;
}

double LambdaActiveSet(
		const vector< vector<double> > &Y,
		double lambda_y,
		smat_t &Lambda,
		const vector< vector<double> > &R,
		CGGMOptions &options) {
	long n_y = Y[0].size();
	long n_x = R[0].size();
	long q = Lambda.p;

	double subgrad = 0;
	vector<Triplet> triplets;
	triplets.reserve(Lambda.nnz);

	#pragma omp parallel for \
		shared(triplets) reduction(+:subgrad) schedule(static)
	for (long i = 0; i < q; i++) {
		vector<double> Sigma_i(q, 0);
		Lambda.ComputeInv(Sigma_i, i, options.grad_tol);
		
		// Off-diagonal elts:
		vector<double> G_i(i, 0);
		for (long j = 0; j < i; j++) {
			double Syy_ij = 0;
			for (long nn = 0; nn < n_y; nn++) {
				Syy_ij += Y[i][nn] * Y[j][nn];
			}
			Syy_ij = (1.0 / n_y) * Syy_ij;
			double Psi_ij = 0;
			for (long nn = 0; nn < n_x; nn++) {
				Psi_ij += R[i][nn] * R[j][nn];
			}
			G_i[j] = Syy_ij - Sigma_i[j] - Psi_ij;
		}

		long idx = Lambda.row_ptr[i];		
		for (long j = 0; j < i; j++) {
			double Lambda_ij = 0;
			if (idx < Lambda.row_ptr[i+1] && Lambda.col_idx[idx] == j) {
				Lambda_ij = Lambda.values[idx];
				idx++;
			}
			double G_ij = G_i[j];
			if (Lambda_ij != 0 || (!options.refit && fabs(G_ij) > lambda_y) ) {
				subgrad += 2*fabs(L1SubGrad(Lambda_ij,G_ij,lambda_y));
				Triplet triplet(i, j, Lambda_ij);
				#pragma omp critical
				triplets.push_back(triplet);
			}
		}

		// Diagonal elts:
		double Syy_ii = 0;
		for (long nn = 0; nn < n_y; nn++) {
			Syy_ii += Y[i][nn] * Y[i][nn];
		}
		Syy_ii = (1.0 / n_y) * Syy_ii;
		double Psi_ii = 0;
		for (long nn = 0; nn < n_x; nn++) {
			Psi_ii += R[i][nn] * R[i][nn];
		}
		double G_ii = Syy_ii - Sigma_i[i] - Psi_ii;
		subgrad += fabs(G_ii);
		Triplet triplet(i, i, Lambda.values[idx]);
		#pragma omp critical
		triplets.push_back(triplet);
	}
	Lambda.setTriplets(q, triplets);
	return subgrad;
}

void ThetaCoordinateDescent(
		const vector< vector<double> > &Y,
		const vector< vector<double> > &X,
		const double lambda_y, // debugging
		const double lambda_x,
		const smat_t &Lambda,
		sparse_t &Theta,
		vector< vector<double> > &Q,
		vector< vector<double> > &R,
		const vector<long> & block_ix,
		const CGGMOptions &options,
		TimeThetaCD &time_report_theta) {
	long n_x = X[0].size();
	long n_y = Y[0].size();
	long n_o = min(n_x, n_y);
	long p = Theta.p;
	long q = Theta.q;
	struct timeval start_time, end_time;
	struct timeval start_time_prep, end_time_prep;

	vector< vector<long> > block_list;
	get_block_list(block_ix, block_list);
	long num_blocks = block_list.size();
	
	vector<long> blocks_permuted(num_blocks, 0);
	for (long i = 0; i < num_blocks; i++) {
		blocks_permuted[i] = i;
	}
	random_shuffle(blocks_permuted.begin(), blocks_permuted.end());

	vector<long> rows_permuted(p, 0);
	for (long i = 0; i < p; i++) {
		rows_permuted[i] = i;
	}
	
	for (vector<long>::iterator block_it = blocks_permuted.begin();
			block_it != blocks_permuted.end(); ++block_it) {
		long curr_block = *block_it;
		long curr_block_size = block_list[curr_block].size();

		// Compute nonzero rows of V for current block
		// All nonzero rows *must be* marked as nonzero
		// Truly zero rows can be marked as nonzero, but reduces speedup
		long nnz_V_block = 0;
		vector<long> nz_V_block_ix(p, -1);
		for (long i = 0; i < p; i++) {
			// Iff row of Theta is empty, row of V is empty
			if (Theta.row_ptr[i] != Theta.row_ptr[i+1]) {
				// Maps actual row index to row in condensed V:
				nz_V_block_ix[i] = nnz_V_block;
				nnz_V_block += 1;
			}
		}
		vector<long> nz_V_block(nnz_V_block);
		for (long i = 0; i < p; i++) {
			// Maps row index for condensed V to actual row index
			long i_ix = nz_V_block_ix[i];
			if (i_ix != -1) {
				nz_V_block[i_ix] = i;
			}
		}

		gettimeofday(&start_time, NULL);
		gettimeofday(&start_time_prep, NULL);
		
		// Compute Sigma for current block
		vector< vector<double> > Sigma_block( // column-major
				curr_block_size, vector<double>(q, 0));
		#pragma omp parallel for schedule(dynamic,64)
		for (long block_col_ix = 0; block_col_ix < curr_block_size;
				block_col_ix++) {
			Lambda.ComputeInv(Sigma_block[block_col_ix],
					block_list[curr_block][block_col_ix], options.grad_tol);
		}

		gettimeofday(&end_time_prep, NULL);
		time_report_theta.inv += toddiff(&start_time_prep, &end_time_prep);
		gettimeofday(&start_time_prep, NULL);

		// Compute V = Theta*Sigma for current block	
		gettimeofday(&start_time_prep, NULL);
		vector< vector<double> > V_block( // column-major
				curr_block_size, vector<double>(nnz_V_block, 0));
	    #pragma omp parallel for schedule(static)
		for (long j = 0; j < curr_block_size; j++) {
			for (long i_ix = 0; i_ix < nnz_V_block; i_ix++) {
				long i = nz_V_block[i_ix];
				double tmp = 0;
				for (long idx = Theta.row_ptr[i]; 
						idx < Theta.row_ptr[i+1]; idx++) {
					long k = Theta.col_idx[idx];
					tmp += Theta.values[idx] * Sigma_block[j][k];
				}
				V_block[j][i_ix] = tmp;
			}
		}

		gettimeofday(&end_time, NULL);
		gettimeofday(&end_time_prep, NULL);
		time_report_theta.prep += toddiff(&start_time, &end_time);
		time_report_theta.vp += toddiff(&start_time_prep, &end_time_prep);

		random_shuffle(rows_permuted.begin(), rows_permuted.end());
		for (vector<long>::iterator row_it = rows_permuted.begin();
				row_it != rows_permuted.end(); ++row_it) {
			long curr_row = *row_it;
			long curr_row_ix = nz_V_block_ix[curr_row];
			
			gettimeofday(&start_time, NULL);

			// Compute current free set size, and skip if free set empty
			long curr_free_size = 0;
			long idx_cnt = Theta.row_ptr[curr_row];
			long block_col_ix_cnt = 0;
			while (idx_cnt < Theta.row_ptr[curr_row+1] 
					&& block_col_ix_cnt < curr_block_size) {
				long free_col = Theta.col_idx[idx_cnt];
				long block_col = block_list[curr_block][block_col_ix_cnt];
				if (free_col < block_col) {
					idx_cnt++;
				} else if (free_col > block_col) {
					block_col_ix_cnt++;
				} else {
					curr_free_size++;
					idx_cnt++;
					block_col_ix_cnt++;
				}
			}
			if (curr_free_size == 0) { 
				continue;
			}
			
			// Compute free set on (current row, current block)
			vector<long> idx_intersect(curr_free_size, 0);
			vector<long> block_col_ix_intersect(curr_free_size, 0);
			long curr_free_ix = 0;
			long idx = Theta.row_ptr[curr_row];
			long block_col_ix = 0;
			while (idx < Theta.row_ptr[curr_row+1] 
					&& block_col_ix < curr_block_size) {
				long free_col = Theta.col_idx[idx];
				long block_col = block_list[curr_block][block_col_ix];
				if (free_col < block_col) {
					idx++;
				} else if (free_col > block_col) {
					block_col_ix++;
				} else {
					idx_intersect[curr_free_ix] = idx;
					block_col_ix_intersect[curr_free_ix] = block_col_ix;
					curr_free_ix++;
					idx++;
					block_col_ix++;
				}
			}
			
			gettimeofday(&start_time_prep, NULL);

			// Compute Sxx_i
			vector<double> Sxx_i(nnz_V_block, 0);
			#pragma omp parallel for schedule(static)
			for (long j_ix = 0; j_ix < nnz_V_block; j_ix++) {
				long j = nz_V_block[j_ix];
				double Sxx_ij = 0;
				for (long nn = 0; nn < n_x; nn++) {
					Sxx_ij += X[curr_row][nn] * X[j][nn];
				}
				Sxx_i[j_ix] = (1.0 / n_x) * Sxx_ij;
			}
			double Sxx_ii = Sxx_i[curr_row_ix];

			// Compute Sxy_ij for free set
		 	vector<double> Sxy(curr_free_size, 0);
			#pragma omp parallel for schedule(static)
			for (long free_ix = 0; free_ix < curr_free_size; free_ix++) {
				long j = Theta.col_idx[idx_intersect[free_ix]];
				double Sxy_ij = 0;
				for (long nn = 0; nn < n_o; nn++) {
					Sxy_ij += X[curr_row][nn] * Y[j][nn];
				}
				Sxy[free_ix] = (1.0 / n_o) * Sxy_ij;
			}

			gettimeofday(&end_time_prep, NULL);
			time_report_theta.s += toddiff(&start_time_prep, &end_time_prep);
			gettimeofday(&start_time_prep, NULL);

			// Compute P_ij for free set
			vector<double> P(curr_free_size, 0);
			#pragma omp parallel for schedule(static)
			for (long free_ix = 0; free_ix < curr_free_size; free_ix++) {
				long block_col_ix = block_col_ix_intersect[free_ix];
				// j == block_list[curr_block][block_col_ix]
				//   == Theta.col_idx[idx_intersect[free_ix]]
				double P_ij = 0;
				for (long j_ix = 0; j_ix < nnz_V_block; j_ix++) {
					P_ij += Sxx_i[j_ix] * V_block[block_col_ix][j_ix];
				}
				P_ij -= Sxx_i[curr_row_ix] * V_block[block_col_ix][curr_row_ix];
				P[free_ix] = P_ij;
			}

			gettimeofday(&end_time, NULL);
			time_report_theta.prep += toddiff(&start_time, &end_time);
			gettimeofday(&start_time, NULL);
			gettimeofday(&end_time_prep, NULL);
			time_report_theta.vp += toddiff(&start_time_prep, &end_time_prep);

			vector<long> row_permute(curr_free_size, 0);
			for (long ix = 0; ix < curr_free_size; ix++) {
				row_permute[ix] = ix;
			}
			random_shuffle(row_permute.begin(), row_permute.end());
			// Perform Theta coordinate descent
			for (long inner = 0; inner < options.Theta_inner_iters; inner++) {
				for (long free_ix = 0; free_ix < curr_free_size; free_ix++) {
					long idx = idx_intersect[free_ix];
					long block_col_ix = block_col_ix_intersect[free_ix];
					long j = Theta.col_idx[idx];
					
					double P_ij = P[free_ix];
					double V_ij = V_block[block_col_ix][curr_row_ix];
					double SV_ij = P_ij + Sxx_ii*V_ij;
					double Sxy_ij = Sxy[free_ix];

					double a = 2*Sigma_block[block_col_ix][j]*Sxx_ii;
					double b = 2*Sxy_ij + 2*SV_ij;
					double c = Theta.values[idx];
					double mu = -c + SoftThreshold(c - b/a, lambda_x/a);

					// Update Theta and V
					Theta.values[idx] += mu;

					#pragma omp parallel for schedule(static)
					for (long t = 0; t < curr_block_size; t++) {
						V_block[t][curr_row_ix] += mu * Sigma_block[t][j];
					}
				}  // loops over coordinates in free set
			}  // loops over Theta inner iter
			gettimeofday(&end_time, NULL);
			time_report_theta.cd += toddiff(&start_time, &end_time);
		}  // loops over rows
	}  // loops over blocks
	
	gettimeofday(&start_time, NULL);
	// Update Q := X*Theta
	Theta.computeLeftProduct(X, Q);
	double scaling = 1.0 / sqrt(n_x);
	#pragma omp parallel for schedule(static)
	for (long j = 0; j < q; j++) {
		for (long i = 0; i < n_x; i++) {
			Q[j][i] = scaling * Q[j][i];
		}
	}
	
	// Update R := Q * Sigma
	#pragma omp parallel for schedule(static)
	for (long i = 0; i < n_x; i++) {
		vector<double> Q_i(q, 0);
		for (long j = 0; j < q; j++) {
			Q_i[j] = Q[j][i];
		}
		vector<double> R_i(q, 0);
		Lambda.ComputeAinvb(Q_i, R_i, options.obj_tol);
		for (long j = 0; j < q; j++) {
			R[j][i] = R_i[j];
		}
	}
	gettimeofday(&end_time, NULL);
	time_report_theta.qr += toddiff(&start_time, &end_time);
}

void LambdaCoordinateDescent(
		vector< vector<double> > &Y,
		vector< vector<double> > &X,
		const double lambda_y,
		const double lambda_x, // debugging
		smat_t &Lambda,
		sparse_t &Theta, // debugging
		vector< vector<double> > &Q,
		vector< vector<double> > &R,
		const vector<long> &block_ix,
		LambdaState &lState,
		const CGGMOptions &options,
		TimeLambdaCD &time_report_lambda) {
	long n_x = X[0].size();
	long n_y = Y[0].size();
	long q = Lambda.p;
	struct timeval start_time, end_time;
	struct timeval start_time_prep, end_time_prep;
	struct timeval start_time_cd, end_time_cd;
	vector< vector<long> > block_list;
	get_block_list(block_ix, block_list);
	long num_blocks = block_list.size();

	vector<long> blocks_z_permuted;
	blocks_z_permuted.reserve(num_blocks);
	for (long i = 0; i < num_blocks; i++) {
		blocks_z_permuted.push_back(i);
	}
	vector<long> blocks_q_permuted = blocks_z_permuted;

	random_shuffle(blocks_z_permuted.begin(), blocks_z_permuted.end());

	// Initialize Delta to all-0 matrix
	smat_t Delta;
	Delta.patternfrom(Lambda);

	// Declare variables to avoid memory re-allocation
	vector<double> Syy, F, G, H;
	vector<long> free_permuted;
	vector<long> row_freeset, idx_freeset, z_ix_freeset, q_ix_freeset;

	for (vector<long>::iterator block_z_it = blocks_z_permuted.begin();
			block_z_it != blocks_z_permuted.end(); ++block_z_it) {
		long block_z = *block_z_it;
		long block_z_size = block_list[block_z].size();
	
		gettimeofday(&start_time, NULL);
		gettimeofday(&start_time_prep, NULL);
		// Compute Sigma and U := Delta*Sigma for block_z
		vector< vector<double> > Sigma_block_z(q, 
			vector<double>(block_z_size, 0)); // row-major
		vector< vector<double> > U_block_z(block_z_size); // column-major
		#pragma omp parallel for schedule(static)
		for (long block_col_ix = 0; block_col_ix < block_z_size;
				block_col_ix++) {
			long j = block_list[block_z][block_col_ix];
			// Compute Sigma column
			vector<double> Sigma_col(q, 0);
			Lambda.ComputeInv(Sigma_col, j,	options.hess_tol);
			for (long i = 0; i < q; i++) {
				Sigma_block_z[i][block_col_ix] = Sigma_col[i];
			}
			// Compute U column
			U_block_z[block_col_ix].resize(q);
			Delta.ComputeAx(Sigma_col, U_block_z[block_col_ix]);
		}
		//PRINTF("finished Sigma_block_z, U_block_z \n");
		//fflush(stdout);

		// Compute Psi := R^T * R for block_z
		vector< vector<double> > Psi_block_z(q); // row-major
		#pragma omp parallel for schedule(static)
		for (long i = 0; i < q; i++) {
			Psi_block_z[i].resize(block_z_size);
			for (long z_ix = 0; z_ix < block_z_size; z_ix++) {
				long j = block_list[block_z][z_ix];
				double tmp = 0;
				for (long nn = 0; nn < n_x; nn++) {
					tmp += R[i][nn] * R[j][nn];
				}
				Psi_block_z[i][z_ix] = tmp;
			}
		}
		//PRINTF("finished Psi_block_z \n");
		//fflush(stdout);
		gettimeofday(&end_time, NULL);
		gettimeofday(&end_time_prep, NULL);
		time_report_lambda.prep += toddiff(&start_time, &end_time);
		time_report_lambda.prep_z += toddiff(&start_time_prep, &end_time_prep);
		
		random_shuffle(blocks_q_permuted.begin(), blocks_q_permuted.end());
		for (vector<long>::iterator block_q_it = blocks_q_permuted.begin();
				block_q_it != blocks_q_permuted.end(); ++block_q_it) {
			long block_q = *block_q_it;
			long block_q_size = block_list[block_q].size();

			// Compute boundary nodes for (block_z, block_q)
			// Also find mapping from column j to q_ix (location in block_q)
			vector<bool> active_cols_in_q(q, false);
			long freeset_size = 0;
			for (long z_ix = 0; z_ix < block_z_size; z_ix++) {
				long i = block_list[block_z][z_ix];
				for (long idx = Lambda.row_ptr[i];
						idx < Lambda.row_ptr[i+1]; idx++) {
					long j = Lambda.col_idx[idx];
					if (block_ix[j] == block_q) {
						active_cols_in_q[j] = true;
						freeset_size++;
					}
				}
			}
			vector<bool> is_boundary(block_q_size, false);  // boundary nodes
			vector<long> col_to_q_ix(q, -1);  // mapping
			long curr_num_boundary = 0;
			for (long q_ix = 0; q_ix < block_q_size; q_ix++) {
				long j = block_list[block_q][q_ix];
				if (active_cols_in_q[j]) {
					is_boundary[q_ix] = true;
					curr_num_boundary++;
				}
				col_to_q_ix[j] = q_ix;
			}

			//if (!options.quiet && block_z != block_q && curr_num_boundary) {
			//	PRINTF("block_z:%ld block_q:%ld curr_num_boundary:%ld \n",
			//		block_z, block_q, curr_num_boundary);
			//}
			
			gettimeofday(&start_time, NULL);
			gettimeofday(&start_time_prep, NULL);
			// Compute Sigma and U for boundary nodes in block_q
			vector< vector<double> > U_block_q(block_q_size); // column-major
			vector< vector<double> > Sigma_block_q(q, 
				vector<double>(block_q_size, 0)); // row-major
			if (block_z == block_q) {
				Sigma_block_q = Sigma_block_z;
				U_block_q = U_block_z;
			} else {
				#pragma omp parallel for schedule(dynamic,64)
				for (long q_ix = 0; q_ix < block_q_size; q_ix++) {
					U_block_q[q_ix].resize(q);
					if (is_boundary[q_ix]) {
						long j = block_list[block_q][q_ix];
						// Compute Sigma column
						vector<double> Sigma_col(q, 0);
						Lambda.ComputeInv(Sigma_col, j,	options.hess_tol);
						for (long i = 0; i < q; i++) {
							Sigma_block_q[i][q_ix] = Sigma_col[i];
						}
						// Compute U column
						Delta.ComputeAx(Sigma_col, U_block_q[q_ix]);
					}
				}
			}
			//PRINTF("finished Sigma_block_q, U_block_q \n");
			//fflush(stdout);

			// Compute Psi for boundary nodes in block_q
			vector< vector<double> > Psi_block_q(q); // row-major
			if (block_z == block_q) {
				Psi_block_q = Psi_block_z;
			} else {
				#pragma omp parallel for schedule(static)
				for (long i = 0; i < q; i++) {
					Psi_block_q[i].resize(block_q_size);
					for (long q_ix = 0; q_ix < block_q_size; q_ix++) {
						if (is_boundary[q_ix]) {
							long j = block_list[block_q][q_ix];
							double tmp = 0;
							for (long nn = 0; nn < n_x; nn++) {
								tmp += R[i][nn] * R[j][nn];
							}
							Psi_block_q[i][q_ix] = tmp;
						}
					}
				}
			}
			gettimeofday(&end_time_prep, NULL);
			time_report_lambda.prep_q += toddiff(&start_time_prep, &end_time_prep);
			//PRINTF("finished Psi_block_q \n");
			//fflush(stdout);

			gettimeofday(&start_time_prep, NULL);
			// Compute free set and mapping
			row_freeset.resize(freeset_size); // row, access into Lambda.row_ptr
			idx_freeset.resize(freeset_size); // idx, access into Lambda.col_idx
			z_ix_freeset.resize(freeset_size); // access into Sigma_block_z etc
			q_ix_freeset.resize(freeset_size); // access into Sigma_block_q etc
			long freeset_counter = 0;
			for (long z_ix = 0; z_ix < block_z_size; z_ix++) {
				long i = block_list[block_z][z_ix];
				for (long idx = Lambda.row_ptr[i]; 
						idx < Lambda.row_ptr[i+1]; idx++) {
					long j = Lambda.col_idx[idx];
					long q_ix = col_to_q_ix[j];
					if (block_ix[j] == block_q) {
						row_freeset[freeset_counter] = i;
						idx_freeset[freeset_counter] = idx;
						z_ix_freeset[freeset_counter] = z_ix;
						q_ix_freeset[freeset_counter] = q_ix;
						freeset_counter++;
					}
				}
			}
			//PRINTF("finished freeset \n");
			//fflush(stdout);

			// Compute F_ij, G_ij, and H_ij for current free set
			F.resize(freeset_size);
			G.resize(freeset_size);
			H.resize(freeset_size);
			#pragma omp parallel for schedule(static)
			for (long free_ix = 0; free_ix < freeset_size; free_ix++) {
				long z_ix = z_ix_freeset[free_ix];
				long q_ix = q_ix_freeset[free_ix];

				double F_ij = 0;
				double G_ij = 0;
				double H_ij = 0;				
				for (long k = 0; k < Lambda.p; k++) {
					if (block_ix[k] == block_z || block_ix[k] == block_q) {
						continue;
					}
					F_ij += Sigma_block_z[k][z_ix] * U_block_q[q_ix][k];
					G_ij += Psi_block_z[k][z_ix] * U_block_q[q_ix][k];
					H_ij += Psi_block_q[k][q_ix] * U_block_z[z_ix][k];
				}
				F[free_ix] = F_ij;
				G[free_ix] = G_ij;
				H[free_ix] = H_ij;
			}
			//PRINTF("finished F, G, H \n");
			//fflush(stdout);

			// Compute Syy_ij for current free set
			Syy.resize(freeset_size);
			#pragma omp parallel for schedule(static)
			for (long free_ix = 0; free_ix < freeset_size; free_ix++) {
				long i = row_freeset[free_ix];
				long j = Lambda.col_idx[idx_freeset[free_ix]];
				double Syy_ij = 0;
				for (long nn = 0; nn < n_y; nn++) {
					Syy_ij += Y[i][nn] * Y[j][nn];
				}
				Syy[free_ix] = (1.0 / n_y) * Syy_ij;
			}
			gettimeofday(&end_time, NULL);
			gettimeofday(&end_time_prep, NULL);
			time_report_lambda.prep += toddiff(&start_time, &end_time);
			time_report_lambda.prep_d += toddiff(&start_time_prep, &end_time_prep);
			//PRINTF("finished Syy \n");
			//fflush(stdout);
	
			// Random shuffle within block
			free_permuted.resize(freeset_size);
			for (long i = 0; i < freeset_size; i++) {
				free_permuted[i] = i;
			}
			random_shuffle(free_permuted.begin(), free_permuted.end());
			
			gettimeofday(&start_time, NULL);
			// Perform Lambda coordinate descent
			//omp_set_dynamic(true);
			for (long inner=0; inner < options.Lambda_inner_iters; inner++) {
				for (long perm_ix = 0; perm_ix < freeset_size; perm_ix++) {
					gettimeofday(&start_time_cd, NULL);
					long free_ix = free_permuted[perm_ix];
					long i = row_freeset[free_ix];
					long idx = idx_freeset[free_ix];
					long j = Lambda.col_idx[idx];
					long z_ix = z_ix_freeset[free_ix];
					long q_ix = q_ix_freeset[free_ix];
					
					double Syy_ij = Syy[free_ix];
					double Lambda_ij = Lambda.values[idx];
					double Delta_ij = Delta.values[idx];
					double Sigma_ii = Sigma_block_z[i][z_ix];
					double Sigma_ij = Sigma_block_q[i][q_ix];
					double Sigma_jj = Sigma_block_q[j][q_ix];
					double Psi_ii = Psi_block_z[i][z_ix];
					double Psi_ij = Psi_block_q[i][q_ix];
					double Psi_jj = Psi_block_q[j][q_ix];
					
					double Sigma_i_U_j = F[free_ix];
					double Psi_i_U_j = G[free_ix];
					double Psi_j_U_i = H[free_ix];
					# pragma omp parallel for \
						reduction(+:Sigma_i_U_j, Psi_i_U_j, Psi_j_U_i) \
						schedule(static)
					for (long cz_ix = 0; cz_ix < block_z_size; cz_ix++) {
						long cz = block_list[block_z][cz_ix];
						Sigma_i_U_j += 
							Sigma_block_z[cz][z_ix] * U_block_q[q_ix][cz];
						Psi_i_U_j +=
							Psi_block_z[cz][z_ix] * U_block_q[q_ix][cz];
						Psi_j_U_i +=
							Psi_block_q[cz][q_ix] * U_block_z[z_ix][cz];
					}
					# pragma omp parallel for \
						reduction(+:Sigma_i_U_j, Psi_i_U_j, Psi_j_U_i) \
						schedule(static)
					for (long cq_ix = 0; cq_ix < block_q_size; cq_ix++) {
						long cq = block_list[block_q][cq_ix];
						Sigma_i_U_j += 
							Sigma_block_z[cq][z_ix] * U_block_q[q_ix][cq];
						Psi_i_U_j +=
							Psi_block_z[cq][z_ix] * U_block_q[q_ix][cq];
						Psi_j_U_i +=
							Psi_block_q[cq][q_ix] * U_block_z[z_ix][cq];
					}
					gettimeofday(&end_time_cd, NULL);
					time_report_lambda.cd_prep += toddiff(
						&start_time_cd, &end_time_cd);
					
					if (i == j) { // diagonal elements
						double a = Sigma_ii*Sigma_ii + 2*Sigma_ii*Psi_ii;
						double b = Syy_ij - Sigma_ij - Psi_ij
							+ Sigma_i_U_j + 2*Psi_i_U_j;
						double mu = -b/a;
						
						Delta.values[idx] += mu;
						#pragma omp parallel for schedule(static)
						for (long t = 0; t < block_z_size; t++) {
							U_block_z[t][i] += mu * Sigma_block_z[i][t];
							U_block_q[t][i] = U_block_z[t][i];
						}
					} else { // off-diagonal elements
						double a = Sigma_ij*Sigma_ij + Sigma_ii*Sigma_jj
							+ Sigma_ii*Psi_jj + 2*Sigma_ij*Psi_ij
							+ Sigma_jj*Psi_ii;
						double b = Syy_ij - Sigma_ij - Psi_ij
							+ Sigma_i_U_j + Psi_i_U_j + Psi_j_U_i;
						double c = Lambda_ij + Delta_ij;
						double mu = -c + SoftThreshold(c - b/a, lambda_y/a);
						
						Delta.values[idx] += mu;
						#pragma omp parallel for schedule(static)
						for (long t = 0; t < block_z_size; t++) {
							U_block_z[t][i] += mu * Sigma_block_z[j][t];
							U_block_z[t][j] += mu * Sigma_block_z[i][t];
						}
						#pragma omp parallel for schedule(static)
						for (long t = 0; t < block_q_size; t++) {
							U_block_q[t][i] += mu * Sigma_block_q[j][t];
							U_block_q[t][j] += mu * Sigma_block_q[i][t];
						}
					}
				}  // loops over coordinates in free set
				//omp_set_dynamic(false);
				//PRINTF("finished updating coordinates in free set \n");
				//fflush(stdout);	
			}  // loops over Lambda inner iter
			gettimeofday(&end_time, NULL);
			time_report_lambda.cd += toddiff(&start_time, &end_time);
		}  // loops over block_q
	}  // loops over over block_z

	
	// Backtracking line search
	gettimeofday(&start_time, NULL); // time_report.ls

	// Only checks for PD-ness if previous lsiters was 0
	if (lState.numNoBacktracking > 2) {
		smat_t Lambda_alpha(Lambda, Delta, 1.0);
		double logdetLambda_alpha;
		int flag_pd = Lambda_alpha.ComputeLogdet(
			logdetLambda_alpha, options.obj_tol);
		if (flag_pd) {
			Lambda.copyfrom(Lambda_alpha);
			double trRtQ_alpha = traceRtQ_alpha(
					Lambda_alpha, Q, options.obj_tol);
			// Update R := Q * Sigma
			#pragma omp parallel for schedule(static)
			for (long i = 0; i < n_x; i++) {
				vector<double> Q_i(q, 0);
				for (long j = 0; j < q; j++) {
					Q_i[j] = Q[j][i];
				}
				vector<double> R_i(q, 0);
				Lambda.ComputeAinvb(Q_i, R_i, options.obj_tol);
				for (long j = 0; j < q; j++) {
					R[j][i] = R_i[j];
				}
			}
			lState = LambdaState(
				logdetLambda_alpha, trRtQ_alpha, lState.numNoBacktracking + 1);
			if (!options.quiet) {
				PRINTF("   line search PD\n");
			}
			gettimeofday(&end_time, NULL);
			time_report_lambda.ls += toddiff(&start_time, &end_time);
			return;
		}
	}

	double alpha = 1;
	bool success = false;
	double trDeltaSigma = traceDeltaSigma(Delta, Lambda, options.obj_tol);
	double trDeltaSyy = (1.0 / n_y) * Delta.traceProduct(Y, Y);
	double trDeltaRtR = Delta.traceProduct(R, R);
	double trGradDelta = trDeltaSyy - trDeltaSigma - trDeltaRtR;
	smat_t LambdaPlusDelta(Lambda, Delta, 1);
	double RHS = options.sigma*(trGradDelta + 
		lambda_y*(LambdaPlusDelta.L1NormOffDiag() - Lambda.L1NormOffDiag()));

	for (int lsiter = 0; lsiter < options.max_ls_iters; lsiter++) {
		smat_t Lambda_alpha(Lambda, Delta, alpha);
		double logdetLambda_alpha;
		int flag_pd = Lambda_alpha.ComputeLogdet(
			logdetLambda_alpha, options.obj_tol);
		if (!flag_pd) {
			if (!options.quiet) {
				PRINTF("   line search %d, alpha=%f, not PD\n", lsiter, alpha);
			}
			alpha *= options.beta;
			continue;
		}

		double trRtQ_alpha = traceRtQ_alpha(Lambda_alpha, Q, options.obj_tol);
		double LHS = -logdetLambda_alpha + lState.logdetLambda
			+ alpha*trDeltaSyy + trRtQ_alpha - lState.trRtQ
			+ lambda_y*(Lambda_alpha.L1NormOffDiag() - Lambda.L1NormOffDiag());

		if (LHS <= alpha*RHS) {
			/*
			if (!options.quiet) {
				PRINTF("   lsearch %d, alpha=%f, sufficient decrease=%f\n",
					lsiter, alpha, LHS);
			}
			*/
			Lambda.copyfrom(Lambda_alpha);
			// Update R := Q * Sigma
			#pragma omp parallel for schedule(static)
			for (long i = 0; i < n_x; i++) {
				vector<double> Q_i(q, 0);
				for (long j = 0; j < q; j++) {
					Q_i[j] = Q[j][i];
				}
				vector<double> R_i(q, 0);
				Lambda.ComputeAinvb(Q_i, R_i, options.obj_tol);
				for (long j = 0; j < q; j++) {
					R[j][i] = R_i[j];
				}
			}
			int numNoBacktracking = 0;
			if (lsiter == 0) {
				numNoBacktracking = lState.numNoBacktracking + 1;
			}
			lState = LambdaState(
				logdetLambda_alpha, trRtQ_alpha, numNoBacktracking);
			success = true;
			break;
		} 
		if (!options.quiet) {
			PRINTF("   lsearch %d, alpha=%f, insufficient decrease=%f\n",
				lsiter, alpha, LHS);
		}
		alpha *= options.beta;
	}
	if (!success && !options.quiet) {
		PRINTF("  Lambda line search failed\n");
	}
	gettimeofday(&end_time, NULL);
	time_report_lambda.ls += toddiff(&start_time, &end_time);
}
	

void CGGMfast(
		vector< vector<double> > &Y,
		vector< vector<double> > &X,
		double lambda_y,
		double lambda_x,
		CGGMOptions &options,
		smat_t &Lambda,
		sparse_t &Theta,
		CGGMStats &stats) {
	long n_x = X[0].size();
	long n_y = Y[0].size();
	long p = X.size();
	long q = Y.size();

	// Set number of blocks with limited memory
	// 16 = 8 bytes per double
	// 300 = 10 nonzeros per q * 3 (i,j,val) * 10x bigger initial nonzeros
	double memory_left = 1.0e6*options.memory_usage/8 
		- n_x*p - 3*n_y*q - 300*q;
	if (memory_left < 0) {
		fprintf(stderr, "memory_left < 0\n");
		exit(1);
	}
	bool auto_blocks_Theta = (options.num_blocks_Theta == 0);
	bool auto_blocks_Lambda = (options.num_blocks_Lambda == 0);
	if (options.num_blocks_Theta < 1) {
		options.num_blocks_Theta = 1 + long((p+q)*q / memory_left);
	}
	if (options.num_blocks_Lambda < 1) {
		options.num_blocks_Lambda = 1 + long(7*q*q / memory_left);
	}
	if (!options.quiet) {
		fprintf(stdout, "min_blocks_Theta:%d min_blocks_Lambda:%d \n", 
			options.num_blocks_Theta, options.num_blocks_Lambda);
		fprintf(stdout, "auto_blocks_Theta:%d auto_blocks_Lambda:%d \n",
			auto_blocks_Theta, auto_blocks_Lambda);
	}
	if (options.num_blocks_Theta < 0) {
		fprintf(stderr, "min_blocks_Theta < 0\n");
		exit(1);
	}
	if (options.num_blocks_Lambda < 0) {
		fprintf(stderr, "min_blocks_Lambda < 0\n");
		exit(1);
	}


	srand(1);
	int num_threads = min(options.max_threads, omp_get_max_threads());
	omp_set_num_threads(num_threads);
	omp_set_dynamic(true);
	if (!options.quiet) {
		PRINTF("num threads:%ld \n", num_threads);
	}
	fflush(stdout);

	struct timeval start_time, current_time;
	gettimeofday(&start_time, NULL);

	// Initialize Q^T := X*Theta
	vector< vector<double> > Q(q, vector<double>(n_x, 0));
	Theta.computeLeftProduct(X, Q);
	double scaling = 1.0 / sqrt(n_x);
	#pragma omp parallel for schedule(static)
	for (long j = 0; j < q; j++) {
		for (long i = 0; i < n_x; i++) {
			Q[j][i] = scaling * Q[j][i];
		}
	}
	
	// Initialize R^T := Q * Sigma
	vector< vector<double> > R(q, vector<double>(n_x, 0));
	#pragma omp parallel for schedule(static)
	for (long i = 0; i < n_x; i++) {
		vector<double> Q_i(q, 0);
		for (long j = 0; j < q; j++) {
			Q_i[j] = Q[j][i];
		}
		vector<double> R_i(q, 0);
		Lambda.ComputeAinvb(Q_i, R_i, options.obj_tol);
		for (long j = 0; j < q; j++) {
			R[j][i] = R_i[j];
		}
	}

	double initLogdet;
	Lambda.ComputeLogdet(initLogdet, options.obj_tol);
	LambdaState lState(initLogdet, traceProduct(Q, R), 0);
	struct timeval mini_start_time, mini_end_time;

	for (int tOuter = 0; tOuter < options.max_outer_iters; tOuter++) {
		gettimeofday(&mini_start_time, NULL);
		double subgradTheta = ThetaActiveSet(Y, X, lambda_x, Theta, R, options);
		gettimeofday(&mini_end_time, NULL);
		double theta_active_time = toddiff(&mini_start_time, &mini_end_time);
		if (!options.quiet) {
			PRINTF("Iteration %d: finished Theta active %ld \n",
				tOuter, Theta.nnz);
			fflush(stdout);
		}

		gettimeofday(&mini_start_time, NULL);
		double subgradLambda = LambdaActiveSet(
			Y, lambda_y, Lambda, R, options);
		gettimeofday(&mini_end_time, NULL);
		double lambda_active_time = toddiff(&mini_start_time, &mini_end_time);
		if (!options.quiet) {
			PRINTF("Iteration %d: finished Lambda active %ld \n",
				tOuter, Lambda.nnz);
			fflush(stdout);
		}

		double normTheta = Theta.L1Norm();
		double normLambdaOff = Lambda.L1NormOffDiag();
		double normLambda = Lambda.L1Norm();
		double subgrad = subgradTheta + subgradLambda;
		double l1norm = normTheta + normLambda;
		if (tOuter > 0) {
			stats.subgrad.push_back(subgrad);
			stats.l1norm.push_back(l1norm);
		}
	
		if (subgrad < options.tol*l1norm) {
			if (!options.quiet) {
				PRINTF("Converged, subgrad=%f, norm=%f\n", subgrad, l1norm);
			}
			break;
		}

		// Graph clustering
		gettimeofday(&mini_start_time, NULL);
		vector<long> block_ix_Lambda(q, 0);
		vector<long> block_ix_Theta(q, 0);
		long boundary_Lambda, boundary_Theta;
		if (auto_blocks_Lambda) {
			boundary_Lambda = Lambda.auto_clustering(
				block_ix_Lambda, options.num_blocks_Lambda, n_y);
		} else {
			boundary_Lambda = Lambda.clustering(
				block_ix_Lambda, options.num_blocks_Lambda);
		}
		if (auto_blocks_Theta) {
			boundary_Theta = Theta.autoClusteringColumns(
				block_ix_Theta, options.num_blocks_Theta, n_x);
		} else {
			boundary_Theta = Theta.clusteringColumns(
				block_ix_Theta, options.num_blocks_Theta);
		}
		long blocks_Lambda = -1;
		for (long i = 0; i < block_ix_Lambda.size(); i++) {
			blocks_Lambda = max(blocks_Lambda, block_ix_Lambda[i]);
		}
		blocks_Lambda += 1;
		long blocks_Theta = -1;
		for (long i = 0; i < block_ix_Theta.size(); i++) {
			blocks_Theta = max(blocks_Theta, block_ix_Theta[i]);
		}
		blocks_Theta += 1;
		gettimeofday(&mini_end_time, NULL);
		double clustering_time = toddiff(&mini_start_time, &mini_end_time);

		// Lambda CD
		if (!options.quiet) {
			PRINTF("Iteration %d: starting Lambda CD \n", tOuter);
			fflush(stdout);
		}
		struct TimeLambdaCD time_report_lambda;
		gettimeofday(&mini_start_time, NULL);
		srand(1);
		LambdaCoordinateDescent(
			Y, X, lambda_y, lambda_x, Lambda, Theta, Q, R,
			block_ix_Lambda, lState, options, time_report_lambda);
		gettimeofday(&mini_end_time, NULL);
		double lambda_cd_time =	toddiff(&mini_start_time, &mini_end_time);

		// Theta CD
		if (!options.quiet) {
			PRINTF("Iteration %d: starting Theta CD \n", tOuter);
			fflush(stdout);
		}
		struct TimeThetaCD time_report_theta;
		gettimeofday(&mini_start_time, NULL);
		srand(1);
		ThetaCoordinateDescent(
			Y, X, lambda_y, lambda_x, Lambda, Theta, Q, R,
			block_ix_Theta, options, time_report_theta);
		lState.trRtQ = traceProduct(R, Q);
		gettimeofday(&mini_end_time, NULL);
		double theta_cd_time = toddiff(&mini_start_time, &mini_end_time);

		if (!options.quiet) {
			PRINTF("Iteration %d, Lambda(subgrad=%f,active=%ld,norm=%f)\n",
				tOuter, subgradLambda, Lambda.nnz, normLambda);
			PRINTF("             Theta(subgrad=%f,active=%ld,norm=%f)\n",
				subgradTheta, Theta.nnz, normTheta);
			PRINTF("             Hessian tol=%e\n", options.hess_tol);
			fflush(stdout);
		}
	
		gettimeofday(&mini_start_time, NULL);
		double f = Objective(
			Y, X, lambda_y, lambda_x, Lambda, Theta, Q, R, options.obj_tol);
		gettimeofday(&mini_end_time, NULL);
		double objective_time = toddiff(&mini_start_time, &mini_end_time);
		
		gettimeofday(&current_time, NULL);
		stats.objval.push_back(f);
		stats.time.push_back(toddiff(&start_time, &current_time));
		stats.active_set_size.push_back((double)(Theta.nnz + Lambda.nnz));
		stats.active_theta.push_back((double)Theta.nnz);
		stats.active_lambda.push_back((double)Lambda.nnz);
		
		// Output iteration time breakdown
		stats.time_lambda_active.push_back(lambda_active_time);
		stats.time_theta_active.push_back(theta_active_time);
		stats.time_clustering.push_back(clustering_time);
		stats.time_objective.push_back(objective_time);

		stats.time_theta_cd.push_back(theta_cd_time);
		stats.time_theta_cd_prep.push_back(time_report_theta.prep);
		stats.time_theta_cd_inv.push_back(time_report_theta.inv);
		stats.time_theta_cd_s.push_back(time_report_theta.s);
		stats.time_theta_cd_vp.push_back(time_report_theta.vp);
		stats.time_theta_cd_cd.push_back(time_report_theta.cd);
		stats.time_theta_cd_qr.push_back(time_report_theta.qr);
		
		stats.time_lambda_cd.push_back(lambda_cd_time);
		stats.time_lambda_cd_prep.push_back(time_report_lambda.prep);
		stats.time_lambda_cd_prep_z.push_back(time_report_lambda.prep_z);
		stats.time_lambda_cd_prep_q.push_back(time_report_lambda.prep_q);
		stats.time_lambda_cd_prep_d.push_back(time_report_lambda.prep_d);
		stats.time_lambda_cd_cd.push_back(time_report_lambda.cd);
		stats.time_lambda_cd_cd_prep.push_back(time_report_lambda.cd_prep);
		double cd_apply = time_report_lambda.cd - time_report_lambda.cd_prep;
		stats.time_lambda_cd_cd_apply.push_back(cd_apply);
		stats.time_lambda_cd_ls.push_back(time_report_lambda.ls);

		stats.boundary_lambda.push_back(boundary_Lambda);
		stats.boundary_theta.push_back(boundary_Theta);
		stats.blocks_lambda.push_back(blocks_Lambda);
		stats.blocks_theta.push_back(blocks_Theta);

	}
	// Model selection statistics
	stats.loglik = LogLik(Y, X, Lambda, Theta, Q, R, options.obj_tol);
	long n_o = min(n_x, n_y);
	long nonzeros = CountNonzeros(Lambda, Theta);
	double dims = p*q + q*(q+1)/2;
	stats.AIC = 2*nonzeros - 2*stats.loglik;
	stats.BIC = log(n_o)*nonzeros - 2*stats.loglik;
	stats.eBIC = log(n_o)*nonzeros + 4*log(dims)*nonzeros - 2*stats.loglik;

	smat_t Lambda_sym;
	Lambda_sym.symmetricfrom(Lambda);
	Lambda.copyfrom(Lambda_sym);
}


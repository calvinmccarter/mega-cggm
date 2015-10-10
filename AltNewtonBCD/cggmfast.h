#ifndef _CGGMFAST_H_
#define _CGGMFAST_H_

#include "smat.h"
#include "sparse.h"

#include <vector>
#include <time.h>

using std::vector;

struct CGGMOptions {
		CGGMOptions() : 
	quiet(false),
	max_outer_iters(50),
	Theta_inner_iters(1),
	Lambda_inner_iters(1),
    max_ls_iters(15),
    beta(0.5),
    sigma(1.0e-4),
    tol(1.0e-2),
	obj_tol(1.0e-13),
	grad_tol(1.0e-10),
	hess_tol(1.0e-8),
 	num_blocks_Theta(-1), // <0: min with memory, 0: auto, >0: user input
 	num_blocks_Lambda(-1), // <0: min with memory, 0: auto, >0: user input
	memory_usage(32000), // default 32 Gb 
	max_threads(1) {}
	
	// Whether to print info messages
	bool quiet;

	// Maximum iterations 
	int max_outer_iters;
	int Theta_inner_iters;
	int Lambda_inner_iters;
	int max_ls_iters;

	// Line search parameters
	double beta;
	double sigma;

	// Tolerance parameter for outer iterations 
	double tol;

	// Conjugate gradient parameters
	double obj_tol; // also for updating Q and R
	double grad_tol;
	double hess_tol;

	// Limited-memory parameters
	long num_blocks_Theta;
	long num_blocks_Lambda;
	long memory_usage;

	// Parallelization parameters
	int max_threads;
};

struct CGGMStats {
	vector<double> objval;
	vector<double> time;
	vector<double> active_set_size;
	vector<double> active_theta;
	vector<double> active_lambda;
	vector<double> subgrad;
	vector<double> l1norm;

	vector<long> boundary_lambda;
	vector<long> boundary_theta;
	vector<long> blocks_lambda;
	vector<long> blocks_theta;
    
	vector<double> time_objective;
	vector<double> time_clustering;
    vector<double> time_theta_active;
    vector<double> time_lambda_active;
    
    vector<double> time_theta_cd;
	// Splitting time_theta_cd:
    vector<double> time_theta_cd_prep; // V, Sigma, Sxx, P, Sxy
    vector<double> time_theta_cd_cd; // updates to Theta, V
    vector<double> time_theta_cd_qr; // updates to QR
	// Splitting time_theta_cd_prep:
	vector<double> time_theta_cd_inv; // Sigma
	vector<double> time_theta_cd_s; // Sxx, Sxy
	vector<double> time_theta_cd_vp; // V, P

    vector<double> time_lambda_cd;
    vector<double> time_lambda_cd_prep;
    vector<double> time_lambda_cd_prep_z;
    vector<double> time_lambda_cd_prep_q;
    vector<double> time_lambda_cd_prep_d;
    vector<double> time_lambda_cd_cd;
	vector<double> time_lambda_cd_cd_prep;
	vector<double> time_lambda_cd_cd_apply;
    vector<double> time_lambda_cd_ls;
};

extern "C" {
void CGGMfast(
	vector< vector<double> > &Y,
	vector< vector<double> > &X,
	double lambda_y,
	double lambda_x,
	CGGMOptions &options,
	smat_t &Lambda,
	sparse_t &Theta,
	CGGMStats &stats);
}

#endif

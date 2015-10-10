// An interface to BigQUIC with raw data as input.
//
// See the README file for more information.

//#include <metis.h>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <string>
#include <vector>
#include "cggmfast.h"
#include "util.h"

using namespace std;

void exit_with_help() {
	printf(
		"Usage: ./hugecggm_run [options] "
		"n_x p n_y q Y_filename X_filename Lambda_filename Theta_filename stats_filename\n"
		"sparse Lambda and Theta format (1-based array indices):\n"
		"    num_rows num_columns num_nonzeros\n"
		"    row column value\n"
		"options:\n"
		"    -y lambda_y : set the regularization parameter lambda_y (default 0.5)\n"    
		"    -x lambda_x : set the regularization parameter lambda_x (default 0.5)\n"    
		"    -L Lambda0_filename: filename with initial Lambda\n"
		"    -T Theta0_filename: filename with initial Theta\n"
		"    -v verbose: show information or not (0 or 1)\n"
		"    -i max_outer_iters: max number of outer iterations\n"
		"    -s sigma: backtracking termination criterion\n"
		"    -q tol: tolerance for terminating outer loop\n"
		"    -o obj_tol: CG tolerance for calculating objective function\n"
		"    -g grad_tol: CG tolerance for calculating gradient\n"
		"    -h hess_tol: CG tolerance for calculating hessian\n"
		"    -l num_blocks_Lambda: number of blocks for Lambda CD\n"
		"    -t num_blocks_Theta: number of blocks for Theta CD\n"
		"    -m memory_usage: memory capacity in MB\n"  
		"    -n threads : set the max number of threads\n"    
	);
	exit(1);
}

int main(int argc, char **argv) {
	double lambda_y = 0.5;
	double lambda_x = 0.5;
	int num_reqd_args = 9;

	CGGMOptions options;
	if (argc < 1 + num_reqd_args) {
		fprintf(stderr,"not enough arguments\n");
		exit_with_help();
	}
	vector<string> cmdargs(argv + 1, argv + argc);
	int num_args = cmdargs.size();
	int num_opts_and_vals = num_args - num_reqd_args;
	int num_opts = (int) num_opts_and_vals / 2;
	if (num_opts_and_vals % 2 != 0) {
		fprintf(stderr,"option is missing a value\n");
		exit_with_help();
	}

	string Lambda0_filename = "";
	string Theta0_filename = "";

	for (int i = 0; i < num_opts; i++) {
		if (cmdargs[2*i][0] != '-') {
			fprintf(stderr,"incorrect option format\n");
			exit_with_help();
		}
		switch (cmdargs[2*i][1]) {
			case 'y':
				lambda_y = atof(cmdargs[2*i+1].c_str());
				break;
			case 'x':
				lambda_x = atof(cmdargs[2*i+1].c_str());
				break;
			case 'L':
				Lambda0_filename = cmdargs[2*i+1];
				break;
			case 'T':
				Theta0_filename = cmdargs[2*i+1];
				break;
			case 'v':
				options.quiet = (atoi(cmdargs[2*i+1].c_str()) == 0);
				break;
			case 'i':
				options.max_outer_iters = atoi(cmdargs[2*i+1].c_str());
				break;
			case 's':
				options.sigma = atof(cmdargs[2*i+1].c_str());
				break;
			case 'q':
				options.tol = atof(cmdargs[2*i+1].c_str());
				break;
			case 'o':
				options.obj_tol = atof(cmdargs[2*i+1].c_str());
				break;
			case 'g':
				options.grad_tol = atof(cmdargs[2*i+1].c_str());
				break;
			case 'h':
				options.hess_tol = atof(cmdargs[2*i+1].c_str());
				break;
			case 'l':
				options.num_blocks_Lambda = atol(cmdargs[2*i+1].c_str());
				break;
			case 't':
				options.num_blocks_Theta = atol(cmdargs[2*i+1].c_str());
				break;
			case 'm':
				options.memory_usage = atol(cmdargs[2*i+1].c_str());
				break;
			case 'n':
				options.max_threads = atoi(cmdargs[2*i+1].c_str());
				break;
			default:
				fprintf(stderr,"unknown option: -%c\n", cmdargs[2*i][1]);
				exit_with_help();
				break;
		}
	}

	long n_x = atol(cmdargs[num_args-9].c_str());
	long p = atol(cmdargs[num_args-8].c_str());
	long n_y = atol(cmdargs[num_args-7].c_str());
	long q = atol(cmdargs[num_args-6].c_str());

	string Y_filename = cmdargs[num_args-5];
	string X_filename = cmdargs[num_args-4];
	string Lambda_filename = cmdargs[num_args-3];
	string Theta_filename = cmdargs[num_args-2];
	string stats_filename = cmdargs[num_args-1];

	if (!options.quiet) {
		fprintf(stdout,"n_x=%li p=%li n_y=%li q=%li Yf=%s Xf=%s Lf=%s Tf=%s sf=%s \n",
			n_x, p, n_y, q, Y_filename.c_str(), X_filename.c_str(), 
			Lambda_filename.c_str(), Theta_filename.c_str(), stats_filename.c_str());
	}

	// Read input data from file
	vector< vector<double> > Y(q, vector<double>(n_y, 0));
	vector< vector<double> > X(p, vector<double>(n_x, 0));
	double val;
	ifstream ifY(Y_filename.c_str(), ifstream::in);
	for (long i = 0; i < n_y; i++) {
		for (long j = 0; j < q; j++) {
			if (!ifY.good()) {
				fprintf(stderr, "line %ld column %ld \n", i, j);
				fprintf(stderr, "error reading Y_file\n");
				exit_with_help();
			}
			ifY >> val;
			Y[j][i] = val;
		}
	}
	ifY.close();
	ifstream ifX(X_filename.c_str(), ifstream::in);
	for (long i = 0; i < n_x; i++) {
		for (long j = 0; j < p; j++) {
			if(!ifX.good()) {
				fprintf(stderr, "error reading X_file\n");
				exit_with_help();
			}
			ifX >> val;
			X[j][i] = val;
		}
	}
	ifX.close();
	if (!options.quiet) {
		fprintf(stdout, "finished reading data\n");
	}

	// Center and scale by 1/sqrt(n)
	//double scaling = 1.0/sqrt(n);
	#pragma omp parallel for schedule(static)
	for (long i = 0; i < q; i++) {
		double mean = 0;
		for (long nn = 0; nn < n_y; nn++) {
			mean += Y[i][nn];
		}
		mean = mean/n_y;
		for (long nn = 0; nn < n_y; nn++) {
			Y[i][nn] = (Y[i][nn] - mean);
		}
	}
	#pragma omp parallel for schedule(static)
	for (long i = 0; i < p; i++) {
		double mean = 0;
		for (long nn = 0; nn < n_x; nn++) {
			mean += X[i][nn];
		}
		mean = mean/n_x;
		for (long nn = 0; nn < n_x; nn++) {
			X[i][nn] = (X[i][nn] - mean);
		}
	}
	if (!options.quiet) {
		fprintf(stdout, "finished normalizing data\n");
	}

	// Initialize parameters
	smat_t Lambda;
	vector<Triplet> triplets;
	for (long i = 0; i < q; i++) {
		double tmp = 0;
		for (long nn = 0; nn < n_y; nn++) {
			tmp += Y[i][nn] * Y[i][nn];
		}
		if (tmp == 0) {
			fprintf(stderr, "an output variable has 0 variance\n");
			exit(1);
		}
		triplets.push_back(Triplet(i, i, 1.0/(0.01 + (1.0/n_y)*tmp)));
	}
	Lambda.setTriplets(q, triplets);
	sparse_t Theta(p, q);

	// Initialize Lambda0 if specified by user
	if (!Lambda0_filename.empty()) {
		ifstream ifL(Lambda0_filename.c_str(), ifstream::in);
		long L0_p, L0_q, L0_nnz;
		ifL >> L0_p >> L0_q >> L0_nnz;
		if (L0_p != q || L0_q != q) {
			fprintf(stderr, "error reading Lambda0_file\n");
			exit(1);
		}
		vector<Triplet> triplets;
		long i, j;
		double val;
		for (long n = 0; n < L0_nnz; n++) {
			ifL >> i >> j >> val;
			if (!ifL.good()) {
				fprintf(stderr, "error reading Lambda0_file\n");
				exit(1);
			}
			if (i >= j) {
				triplets.push_back(Triplet(i-1, j-1, val));
			}
		}
		Lambda.setTriplets(q, triplets);
		ifL.close();
	}

	// Initialize Theta0 if specified by user
	if (!Theta0_filename.empty()) {
		ifstream ifT(Theta0_filename.c_str(), ifstream::in);
		long T0_p, T0_q, T0_nnz;
		ifT >> T0_p >> T0_q >> T0_nnz;
		if (T0_p != p || T0_q != q) {
			fprintf(stderr, "error reading Theta0_file\n");
			exit(1);
		}
		vector<Triplet> triplets;
		long i, j;
		double val;
		for (long n = 0; n < T0_nnz; n++) {
			ifT >> i >> j >> val;
			if (!ifT.good()) {
				fprintf(stderr, "error reading Theta0_file\n");
				exit(1);
			}
			triplets.push_back(Triplet(i-1, j-1, val));
		}
		Theta.setTriplets(p, q, triplets);
		ifT.close();
	}

	// Run optimization
	fflush(stdout);
	CGGMStats stats;
	CGGMfast(Y, X, lambda_y, lambda_x, options, Lambda, Theta, stats);

	// Output sparse Lambda
	long true_Lambda_nnz = 0;
	for (long idx = 0; idx < Lambda.nnz; idx++) {
		if (Lambda.values[idx] != 0) {
			true_Lambda_nnz++;
		}
	}
	ofstream fL(Lambda_filename.c_str(), ofstream::out);
	fL.precision(12);
	fL << q << " " << q << " " << true_Lambda_nnz << endl;
	for (long i = 0; i < q; i++) {
		for (long idx = Lambda.row_ptr[i]; idx < Lambda.row_ptr[i+1]; idx++) {
			if (Lambda.values[idx] != 0) {
				fL << i+1 << " " << Lambda.col_idx[idx]+1 << " " 
					<< Lambda.values[idx] << endl;
			}
		}
	}
	fL.close();

	// Output sparse Theta
	long true_Theta_nnz = 0;
	for (long idx = 0; idx < Theta.nnz; idx++) {
		if (Theta.values[idx] != 0) {
			true_Theta_nnz++;
		}
	}
	ofstream fT(Theta_filename.c_str(), ofstream::out);
	fT.precision(12);
	fT << p << " " << q << " " << true_Theta_nnz << endl;
	for (long i = 0; i < p; i++) {
		for (long idx = Theta.row_ptr[i]; idx < Theta.row_ptr[i+1]; idx++) {
			if (Theta.values[idx] != 0) {
				fT << i+1 << " " << Theta.col_idx[idx]+1 << " " 
					<< Theta.values[idx] << endl;
			}
		}
	}
	fT.close();

	// Output stats
	ofstream fS(stats_filename.c_str(), ofstream::out);
	fS.precision(13);
	
	fS << "objval ";
	for (int i = 0; i < stats.objval.size(); i++) {
		fS << stats.objval[i] << " ";
	}
	fS << endl;
	fS << "time ";
	for (int i = 0; i < stats.time.size(); i++) {
		fS << stats.time[i] << " ";
	}
	fS << endl;
	fS << "active_set_size ";
	for (int i = 0; i < stats.active_set_size.size(); i++) {
		fS << stats.active_set_size[i] << " ";
	}
	fS << endl;
	fS << "active_theta ";
	for (int i = 0; i < stats.active_theta.size(); i++) {
		fS << stats.active_theta[i] << " ";
	}
	fS << endl;
	fS << "active_lambda ";
	for (int i = 0; i < stats.active_lambda.size(); i++) {
		fS << stats.active_lambda[i] << " ";
	}
	fS << endl;
	fS << "subgrad ";
	for (int i = 0; i < stats.subgrad.size(); i++) {
		fS << stats.subgrad[i] << " ";
	}
	fS << endl;
	fS << "l1norm ";
	for (int i = 0; i < stats.l1norm.size(); i++) {
		fS << stats.l1norm[i] << " ";
	}
	fS << endl;
	fS << "time_objective ";
	for (int i = 0; i < stats.time_objective.size(); i++) {
		fS << stats.time_objective[i] << " ";
	}
	fS << endl;
	fS << "time_clustering ";
	for (int i = 0; i < stats.time_clustering.size(); i++) {
		fS << stats.time_clustering[i] << " ";
	}
	fS << endl;
	fS << "time_theta_active ";
	for (int i = 0; i < stats.time_theta_active.size(); i++) {
		fS << stats.time_theta_active[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_active ";
	for (int i = 0; i < stats.time_lambda_active.size(); i++) {
		fS << stats.time_lambda_active[i] << " ";
	}
	fS << endl;
	fS << "time_theta_cd ";
	for (int i = 0; i < stats.time_theta_cd.size(); i++) {
		fS << stats.time_theta_cd[i] << " ";
	}
	fS << endl;
	fS << "time_theta_cd_prep ";
	for (int i = 0; i < stats.time_theta_cd_prep.size(); i++) {
		fS << stats.time_theta_cd_prep[i] << " ";
	}
	fS << endl;
	fS << "time_theta_cd_cd ";
	for (int i = 0; i < stats.time_theta_cd_cd.size(); i++) {
		fS << stats.time_theta_cd_cd[i] << " ";
	}
	fS << endl;
	fS << "time_theta_cd_qr ";
	for (int i = 0; i < stats.time_theta_cd_qr.size(); i++) {
		fS << stats.time_theta_cd_qr[i] << " ";
	}
	fS << endl;
	fS << "time_theta_cd_inv ";
	for (int i = 0; i < stats.time_theta_cd_inv.size(); i++) {
		fS << stats.time_theta_cd_inv[i] << " ";
	}
	fS << endl;
	fS << "time_theta_cd_s ";
	for (int i = 0; i < stats.time_theta_cd_s.size(); i++) {
		fS << stats.time_theta_cd_s[i] << " ";
	}
	fS << endl;
	fS << "time_theta_cd_vp ";
	for (int i = 0; i < stats.time_theta_cd_vp.size(); i++) {
		fS << stats.time_theta_cd_vp[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd ";
	for (int i = 0; i < stats.time_lambda_cd.size(); i++) {
		fS << stats.time_lambda_cd[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_prep ";
	for (int i = 0; i < stats.time_lambda_cd_prep.size(); i++) {
		fS << stats.time_lambda_cd_prep[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_prep_z ";
	for (int i = 0; i < stats.time_lambda_cd_prep_z.size(); i++) {
		fS << stats.time_lambda_cd_prep_z[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_prep_q ";
	for (int i = 0; i < stats.time_lambda_cd_prep_q.size(); i++) {
		fS << stats.time_lambda_cd_prep_q[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_prep_d ";
	for (int i = 0; i < stats.time_lambda_cd_prep_d.size(); i++) {
		fS << stats.time_lambda_cd_prep_d[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_cd ";
	for (int i = 0; i < stats.time_lambda_cd_cd.size(); i++) {
		fS << stats.time_lambda_cd_cd[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_cd_prep ";
	for (int i = 0; i < stats.time_lambda_cd_cd_prep.size(); i++) {
		fS << stats.time_lambda_cd_cd_prep[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_cd_apply ";
	for (int i = 0; i < stats.time_lambda_cd_cd_apply.size(); i++) {
		fS << stats.time_lambda_cd_cd_apply[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_ls ";
	for (int i = 0; i < stats.time_lambda_cd_ls.size(); i++) {
		fS << stats.time_lambda_cd_ls[i] << " ";
	}
	fS << endl;
	fS << "boundary_lambda ";
	for (int i = 0; i < stats.boundary_lambda.size(); i++) {
		fS << stats.boundary_lambda[i] << " ";
	}
	fS << endl;
	fS << "boundary_theta ";
	for (int i = 0; i < stats.boundary_theta.size(); i++) {
		fS << stats.boundary_theta[i] << " ";
	}
	fS << endl;
	fS << "blocks_lambda ";
	for (int i = 0; i < stats.blocks_lambda.size(); i++) {
		fS << stats.blocks_lambda[i] << " ";
	}
	fS << endl;
	fS << "blocks_theta ";
	for (int i = 0; i < stats.blocks_theta.size(); i++) {
		fS << stats.blocks_theta[i] << " ";
	}
	fS << endl;

	fS.close();

	// TODO- use hdf5 format
	return 0;
}


#include "sparse.h"

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
#include <limits>
#include <math.h>
#include <metis.h>
#include <omp.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <time.h>

using namespace std;

sparse_t::sparse_t() {
	p = 0;
	q = 0;
	nnz = 0;
}

sparse_t::sparse_t(long p_, long q_) {
	sparse_t();
	p = p_;
	q = q_;
	nnz = 0;
	row_ptr.resize(p+1);
}

sparse_t::sparse_t(sparse_t &X) {
	sparse_t();
	p = X.p;
	q = X.q;
	nnz = X.nnz;
	values = X.values;
	row_ptr = X.row_ptr;
	col_idx = X.col_idx;
}

void sparse_t::setZeros(long p_, long q_) {
	p = p_;
	q = q_;
	nnz = 0;
	values.clear();
	row_ptr.clear();
	row_ptr.resize(p+1);
	col_idx.clear();
}

void sparse_t::setTriplets(long p_,
						   long q_,
						   const vector<Triplet> &triplets_) {
	// assumes max(rows) < p_, max(cols) < q_
	p = p_;
	q = q_;
	vector<Triplet> triplets = triplets_;
	nnz = triplets.size();
	
	sort(triplets.begin(), triplets.end(), by_row_then_col);

	// Inserts from (row-major) sorted triplets
	row_ptr.clear();
	row_ptr.resize(p+1);
	col_idx.resize(nnz);
	values.resize(nnz);

	long prev_row = -1;
	for (long n = 0; n < nnz; n++) {
		long curr_row = triplets[n].row;
		for (long idx = prev_row + 1; idx <= curr_row; ++idx) {
			row_ptr[idx] = n;
		}
		prev_row = curr_row;

		col_idx[n] = triplets[n].col;
		values[n] = triplets[n].val;
	}
	for (long idx = prev_row + 1; idx <= p; ++idx) {
		row_ptr[idx] = nnz;
	}
}

void sparse_t::setSortedTriplets(long p_,
						   long q_,
						   const vector<Triplet> &triplets) {
	// assumes max(rows) < p_, max(cols) < q_
	p = p_;
	q = q_;
	nnz = triplets.size();

	// Inserts from (row-major) sorted triplets
	row_ptr.clear();
	row_ptr.resize(p+1);
	col_idx.resize(nnz);
	values.resize(nnz);

	long prev_row = -1;
	for (long n = 0; n < nnz; n++) {
		long curr_row = triplets[n].row;
		for (long idx = prev_row + 1; idx <= curr_row; ++idx) {
			row_ptr[idx] = n;
		}
		prev_row = curr_row;

		col_idx[n] = triplets[n].col;
		values[n] = triplets[n].val;
	}
	for (long idx = prev_row + 1; idx <= p; ++idx) {
		row_ptr[idx] = nnz;
	}
}

void sparse_t::copyfrom(const sparse_t &X) {
	p = X.p;
	q = X.q;
	nnz = X.nnz;
	values = X.values;
	row_ptr = X.row_ptr;
	col_idx = X.col_idx;
}

void sparse_t::transposefrom(const sparse_t &X) {
	p = X.q;
	q = X.p;
	nnz = X.nnz;
	
	// Transposes X's triplets, then sorts by-row-then-column
	vector<Triplet> triplets;
	X.getTriplets(triplets);
    for (long n = 0; n < nnz; n++) {
        long tmp = triplets[n].row;
        triplets[n].row = triplets[n].col;
        triplets[n].col = tmp;
    }
	sort(triplets.begin(), triplets.end(), by_row_then_col);
	
	// Inserts from (row-major) sorted triplets
	row_ptr.clear();
	row_ptr.resize(p+1);
	col_idx.clear();
	col_idx.resize(nnz);
	values.clear();
	values.resize(nnz);

	long prev_row = -1;
	for (long n = 0; n < nnz; n++) {
		long curr_row = triplets[n].row;
		for (long idx = prev_row + 1; idx <= curr_row; ++idx) {
			row_ptr[idx] = n;
		}
		prev_row = curr_row;

		col_idx[n] = triplets[n].col;
		values[n] = triplets[n].val;
	}
	for (long idx = prev_row + 1; idx <= p; ++idx) {
		row_ptr[idx] = nnz;
	}
}

void sparse_t::productfrom(const sparse_t &X, const sparse_t &Y) {
	p = X.p;
	q = Y.q;
	vector<Triplet> triplets;

	// elt (i,j) uses ith row of X and jth row of Ytrans
	sparse_t Ytrans;
	Ytrans.transposefrom(Y);

	// Iterates by rows, then columns
	for (long i = 0; i < p; i++) {
		for (long j = 0; j < q; j++) {
			// Performs sparse dot product, similar to merge from mergesort
			double result = 0;
			long idx_X = X.row_ptr[i];
			long idx_end_X = X.row_ptr[i+1];
			long idx_Y = Y.row_ptr[i];
			long idx_end_Y = Y.row_ptr[i+1];
			while (idx_X < idx_end_X && idx_Y < idx_end_Y) {
				if (X.col_idx[idx_X] == Y.col_idx[idx_Y]) {
					result += X.values[idx_X] * Y.values[idx_Y];
					idx_X++;
					idx_Y++;
				} else if (X.col_idx[idx_X] < Y.col_idx[idx_Y]) {
					idx_X++;
				} else {
					idx_Y++;
				}
			}
			if (result != 0) {
				triplets.push_back(Triplet(i, j, result));
			}
		}
	}

	// Inserts from (row-major) sorted triplets
	row_ptr.clear();
	row_ptr.resize(p+1);
	col_idx.clear();
	col_idx.resize(nnz);
	values.clear();
	values.resize(nnz);

	long prev_row = -1;
	for (long n = 0; n < nnz; n++) {
		long curr_row = triplets[n].row;
		for (long idx = prev_row + 1; idx <= curr_row; ++idx) {
			row_ptr[idx] = n;
		}
		prev_row = curr_row;

		col_idx[n] = triplets[n].col;
		values[n] = triplets[n].val;
	}
	for (long idx = prev_row + 1; idx <= p; ++idx) {
		row_ptr[idx] = nnz;
	}
}

bool sparse_t::isZero() const {
	return (nnz == 0);
}

double sparse_t::L1Norm() const {
	double result = 0;
	for (long i = 0; i < p; i++) {
		for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
			result += fabs(values[idx]);
		}
	}
	return result;
}

void sparse_t::getTriplets(vector<Triplet> &triplets) const {
	triplets.clear();
	triplets.reserve(nnz);
	for (long i = 0; i < p; i++) {
		for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
			triplets.push_back(Triplet(i, col_idx[idx], values[idx]));
		}
	}
}
	
void sparse_t::print(FILE* fp) const {
	//fprintf(fp, "p: %ld, q: %ld, nnz: %ld\n", p, q, nnz);
	PRINTF("p: %ld, q: %ld, nnz: %ld\n", p, q, nnz);
	for (long i = 0; i < p; i++) {
		for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
			//fprintf(fp, "%ld %ld %.10lf\n",
			//	i+1, col_idx[idx]+1, values[idx]);
			PRINTF("%ld %ld %.10lf\n",
				i+1, col_idx[idx]+1, values[idx]);

		}
	}
}

// B(n,q) = A(n,p) * This(p,q), where A and B are column-major 
// A == X', B == Q', This == Theta
void sparse_t::computeLeftProduct(
		const vector< vector<double> > &A, vector< vector<double> > &B) const {
	long n = A[0].size();
	B.resize(q);
	for (long i = 0; i < q; i++) {
		B[i].resize(n);
	}
	
	sparse_t tran;
	tran.transposefrom(*this);

    #pragma omp parallel for schedule(static)
	for (long i = 0; i < q; i++) {
		for (long j = 0; j < n; j++) {
			double tmp = 0;
			for (long idx = tran.row_ptr[i]; idx < tran.row_ptr[i+1]; idx++) {
				long k = tran.col_idx[idx];
				tmp += A[k][j] * tran.values[idx];
			}
			B[i][j] = tmp;
		}
	}
}

// B(p,n) = This(p,q) * A(q,n), where A is column-major and B is row-major
void sparse_t::computeRightProduct(
		const vector< vector<double> > &A, vector< vector<double> > &B) const {
	long n = A.size();
	long q = A[0].size();
	if (q == 0) {
		B.resize(0);
		return;
	}

	vector<long> nz_V_block;
	nz_V_block.reserve(p);
	for (long i = 0; i < p; i++) {
		// Iff row of Theta is empty, row of V is empty
		if (row_ptr[i] != row_ptr[i+1]) {
			nz_V_block.push_back(i);
		}
	}
	long nnz_V_block = nz_V_block.size();

    #pragma omp parallel for schedule(dynamic,64)
	for (long i_ix = 0; i_ix < nnz_V_block; i_ix++) {
		long i = nz_V_block[i_ix];
		for (long j = 0; j < n; j++) {
			double tmp = 0;
			for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
				long k = col_idx[idx];
				tmp += values[idx] * A[j][k];
			}
			B[i][j] = tmp;
		}
	}
}

double sparse_t::traceProduct(const vector< vector<double> > &X) const { 
	double result = 0;
	for (long i = 0; i < p; i++) {
		for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
			result += values[idx] * X[i][col_idx[idx]];
		}
	}
	return result;
}

// Both X and Y are column-major
double sparse_t::traceProduct(const vector< vector<double> > &X,
							  const vector< vector<double> > &Y) const {
	long n_x = X[0].size();
	long n_y = Y[0].size();
	long n_o = min(n_x, n_y);
	double result = 0;
    #pragma omp parallel for reduction(+:result)
	for (long i = 0; i < p; i++) {
		for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
			long j = col_idx[idx];
			double xy_ij = 0;
			for (long nn = 0; nn < n_o; nn++) {
				xy_ij += X[i][nn] * Y[j][nn];
			}
			result += values[idx] * xy_ij;
		}
	}
	return result;
}

long sparse_t::clusteringColumns(vector<long> &block_ind, long nblock) const {
	block_ind.resize(q);
	for (long i = 0; i < q; i++) {
		block_ind[i] = 0;
	}
	if (nblock == 1) {
		long nnz_rows = 0;
		for (long i = 0; i < p; i++) {
			if (row_ptr[i] != row_ptr[i+1]) {
				nnz_rows++;
			}
		}
		return nnz_rows;
	}
	sparse_t tran;
	tran.transposefrom(*this);
	vector<Triplet> graph_triplets;
    #pragma omp parallel for shared(graph_triplets) schedule(dynamic,128)
	for (long i = 0; i < q; i++) {
		for (long j = 0; j < q; j++) {
			long idx_i = tran.row_ptr[i];
			long idx_j = tran.row_ptr[j];
			bool success = false;
			while (idx_i < tran.row_ptr[i+1] && idx_j < tran.row_ptr[j+1]) {
				if (tran.col_idx[idx_i] < tran.col_idx[idx_j]) {
					idx_i++;
				} else if (tran.col_idx[idx_i] > tran.col_idx[idx_j]) {
					idx_j++;
				} else {
					success = true;
					break;
				}
			}
			if (success) {
				# pragma omp critical
				graph_triplets.push_back(Triplet(i, j, 1));
			}
		}
	}
	sparse_t graph;
	graph.setTriplets(q, q, graph_triplets);
	
	// Metis inputs
	idx_t nvtxs = q;  //number of vertices
	idx_t ncon = 1;  // number of balancing constraints
	idx_t *xadj = (idx_t *)malloc(sizeof(idx_t)*(graph.p+1));
	for (long i = 0; i <= graph.p; i++) {
		xadj[i] = graph.row_ptr[i];
	}
	idx_t *adjncy = (idx_t *)malloc(sizeof(idx_t)*graph.nnz);
	for (long idx = 0; idx < graph.nnz; idx++) {
		adjncy[idx] = graph.col_idx[idx];
	}
	idx_t nparts = nblock;
	idx_t options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // edgecut minimization
	//options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM; // rand matching coarsening
	//options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_EDGE; // initial partitioning
	//options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM; // refinement algorithm
	//options[METIS_OPTION_NO2HOP] = 1;  // 2-hop is useful for power-law degree
	//options[METIS_OPTION_NCUTS] = 1; // default number of different partitions
	//options[METIS_OPTION_NITER] = 10; // default iters per level
	options[METIS_OPTION_UFACTOR] = 50; // default load tolerance
	//options[METIS_OPTION_MINCONN] = 0; // don't minimize deg subdomain graph
	//options[METIS_OPTION_CONTIG] = 0; // allow non-contiguous partitions
	options[METIS_OPTION_SEED] = 0; // random number generator seed
	options[METIS_OPTION_NUMBERING] = 0;  // partition indices start from 0
	options[METIS_OPTION_DBGLVL] = 0; // no debugging info
	
	// Metis outputs
	idx_t objval;  // stores edge-cut of solution
	idx_t *part = (idx_t *)malloc(sizeof(idx_t)*graph.p);
	
	// Run graph partitioning
	int flag = METIS_PartGraphKway(
		&nvtxs, &ncon, xadj, adjncy,
		NULL, NULL, NULL, &nparts, NULL,
		NULL, options, &objval, part);
	if (flag != METIS_OK) {
		PRINTF("Metis error %d \n", flag);
	}
	for (long i = 0; i < graph.p; i++) {
		block_ind[i] = part[i];
	}
	free(xadj);
	free(adjncy);
	free(part);

	// Compute number of boundary nodes
	vector< vector<long> > block_list;
	get_block_list(block_ind, block_list);
	vector< vector<bool> > is_boundary(nblock, vector<bool>(p, false));
	vector<long> num_boundary(nblock);
	for (long i = 0; i < p; i++) {
		for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
			long j = col_idx[idx];
			long block_j = block_ind[j];
			if (!is_boundary[block_j][i]) {
				num_boundary[block_j]++;
				is_boundary[block_j][i] = true;
			}
		}
	}
	long sum_boundary = 0;
	for (long i = 0; i < nblock; i++) {
		sum_boundary += num_boundary[i];
	}
	return sum_boundary;
}

long sparse_t::autoClusteringColumns(
		vector<long> &best_block_ind, long min_nblock, long n) const {
	long nnz_rows = 0;
	for (long i = 0; i < p; i++) {
		if (row_ptr[i] != row_ptr[i+1]) {
			nnz_rows++;
		}
	}

	long max_iters = 7;
	double best_cost = numeric_limits<double>::max();
	long best_boundary;

	long curr_k = min_nblock;
	for (long i = 0; i < max_iters; i++) {
		if (curr_k > q) {
			break;
		}
		vector<long> curr_block_ind;
		long curr_boundary = clusteringColumns(curr_block_ind, curr_k);
		double curr_cost = 
			curr_boundary*n*nnz_rows + 2*nnz*q/(double(curr_k));
		PRINTF("curr_k:%ld curr_boundary:%6ld curr_cost:%.4e \n",
			curr_k, curr_boundary, curr_cost);
		if (curr_cost > best_cost) {
			break;
		}
		best_boundary = curr_boundary;
		best_cost = curr_cost;
		best_block_ind = curr_block_ind;
		curr_k *= 2;
	}
	return best_boundary;
}

long sparse_t::clusteringRows(vector<long> &block_ind, long nblock) const {
	// XXX
	return 0;
}

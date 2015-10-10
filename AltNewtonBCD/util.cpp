#include "util.h"

#if defined(LANG_M) || defined(MATLAB_MEX_FILE)
#include <mex.h>
#define PRINTF mexPrintf
#else
#include <stdio.h>
#define PRINTF printf
#endif

#include <math.h>
#include <omp.h>
#include <sys/time.h>
#include <utility>
#include <vector>

using namespace std;

// Returns difference in seconds
double toddiff(struct timeval *start, struct timeval *end) {
	long long tstart = start->tv_sec * 1000000 + start->tv_usec;
	long long tend = end->tv_sec * 1000000 + end->tv_usec;
	return ((double)(tend - tstart))/1000000.0;
}

// Returns list of sorted indices for each block
void get_block_list(
		const vector<long> &block_ix, vector< vector<long> > &block_list) {
	long max_ix = -1;
	for (long i = 0; i < block_ix.size(); i++) {
		max_ix = max(max_ix, block_ix[i]);
	}
	long num_blocks = max_ix + 1;
	block_list.resize(num_blocks);
	for (long i = 0; i < num_blocks; i++) {
		block_list[i].reserve(block_ix.size()/num_blocks + 1);
		block_list[i].resize(0);
	}
	for (long i = 0; i < block_ix.size(); i++) {
		block_list[block_ix[i]].push_back(i);
	}
}

bool by_row_then_col(const Triplet &left, const Triplet &right) {
	if (left.row == right.row) {
		return left.col < right.col;
	}
	return left.row < right.row;
}

bool by_col_then_row(const Triplet &left, const Triplet &right) {
	if (left.col == right.col) {
		return left.row < right.row;
	}
	return left.col < right.col;
}

bool by_violation(const pair<Triplet, double> &left, 
		const pair<Triplet, double> &right) {
	return left.second > right.second;
}

void print(const vector<Triplet> &triplets) {
	for (long i = 0; i < triplets.size(); i++) {
		PRINTF("%ld %ld %.10lf\n", 
			triplets[i].row+1, triplets[i].col+1, triplets[i].val);
	}
}


// Row-major
void print(const vector< vector<double> > &A) {
	PRINTF("p: %ld", A.size());
	if (A.size()) {
		PRINTF(" q: %ld \n", A[0].size());
	}
	for (long i	= 0; i < A.size(); i++) {
		for (long j = 0; j < A[0].size(); j++) {
			PRINTF("%.2lf ", A[i][j]);
		}
		PRINTF("\n");
	}
}

double L1Norm(const vector< vector<double> > &A) {
	double result = 0;
	for (long i = 0; i < A.size(); i++) {
		for (long j = 0; j < A[i].size(); j++) {
			result += fabs(A[i][j]);
		}
	}
	return result;
}

void dense_plus(
		const vector< vector<double> > &A,
		const vector< vector<double> > &B,
		vector< vector<double> > &C) {
	long n = A.size();
	C.resize(n);
	for (long i = 0; i < n; i++) {
		C[i].resize(A[i].size());
	}
	for (long i = 0; i < n; i++) {
		for (long j = 0; j < A[i].size(); j++) {
			C[i][j] = A[i][j] + B[i][j];
		}
	}
}


// C(n,q) = A(n,p) * B(p,q), all row-major
void dense_times(
		const vector< vector<double> > &A,
		const vector< vector<double> > &B,
		vector< vector<double> > &C) {
	long n = A.size();
	long p = B.size();
	if (n == 0 || p == 0) {
		C.clear();
		return;
	}

	long q = B[0].size();
	C.resize(n);
	for (long i = 0; i < n; i++) {
		C[i].resize(q);
	}
	for (long i = 0; i < n; i++) { // row
		for (long j = 0; j < q; j++) { // column
			double tmp = 0;
			for (long k = 0; k < p; k++) {
				// ith row of A, jth column of B
				tmp += A[i][k]*B[k][j];
			}
			C[i][j] = tmp;
		}
	}
}

double traceProduct(
		const vector< vector<double> > &A,
		const vector< vector<double> > &B) {
	double result = 0;
    #pragma omp parallel for reduction(+:result)
	for (long i = 0; i < A.size(); i++) {
        double tmp = 0;
		for (long j = 0; j < A[i].size(); j++) {
			tmp += A[i][j] * B[i][j];
		}
        result += tmp;
	}
	return result;	
}


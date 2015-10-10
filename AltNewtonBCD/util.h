#ifndef _UTIL_H_
#define _UTIL_H_

#include <sys/time.h>
#include <utility>
#include <vector>

double toddiff(struct timeval *start, struct timeval *end);

// Block operations

// Returns list of sorted indices for each block
void get_block_list(
		const std::vector<long> &block_ix, 
		std::vector< std::vector<long> > &block_list);

// Triplets
struct Triplet {
	Triplet(const long row_,
			const long col_,
			const double val_)
		: row(row_), col(col_), val(val_) {}
	long row;
	long col;
	double val;
};

bool by_row_then_col(const Triplet &left, const Triplet &right);

bool by_col_then_row(const Triplet &left, const Triplet &right);

bool by_violation(const std::pair<Triplet, double> &left, 
	const std::pair<Triplet, double> &right);

void print(const std::vector<Triplet> &triplets);


// Dense matrices

void print(const std::vector< std::vector<double> > &A);

double L1Norm(const std::vector< std::vector<double> > &A);

// C = A + B, all row-major
void dense_plus(
		const std::vector< std::vector<double> > &A,
		const std::vector< std::vector<double> > &B,
		std::vector< std::vector<double> > &C);

// C = A * B, all row-major
void dense_times(
		const std::vector< std::vector<double> > &A,
		const std::vector< std::vector<double> > &B,
		std::vector< std::vector<double> > &C);

// sum_ij( [A]_ij * [B]_ij ) == sum_ij( [A^T]_ij * [B^T]_ij )
double traceProduct(
		const std::vector< std::vector<double> > &A,
		const std::vector< std::vector<double> > &B);

#endif

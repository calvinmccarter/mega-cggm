#ifndef _SPARSE_H_
#define _SPARSE_H_

#include "util.h"

#include <vector>
#include <stdio.h>

using std::vector;

// Compressed sparse row format
class sparse_t {
	public:
		long p;
		long q;
		long nnz;
		vector<double> values;
		vector<long> row_ptr;
		vector<long> col_idx;

		// Constructors
		sparse_t();
		sparse_t(long p_, long q_);
		sparse_t(sparse_t &X);

		void setZeros(const long p_, const long q_);
		void setTriplets(
				const long p_,
				const long q_,
				const vector<Triplet> &triplets_);
        // Where triplets_ are already sorted by row, then column:
		void setSortedTriplets(
				const long p_,
				const long q_,
				const vector<Triplet> &triplets_);
		void copyfrom(const sparse_t &X);

		// Transpose and sparse product
		// aliasing-unsafe:
		void transposefrom(const sparse_t &X);
		// aliasing-unsafe:
		void productfrom(const sparse_t &X, const sparse_t &Y);

		bool isZero() const;
		double L1Norm() const;
		void getTriplets(vector<Triplet> &triplets) const;
		void print(FILE* fp) const;

		// B = A * This, where A and B are column-major
		void computeLeftProduct(const vector< vector<double> > &A,
								vector< vector<double> > &B) const;
		// B = This * A, where A is column-major and B is row-major
		void computeRightProduct(const vector< vector<double> > &A,
								 vector< vector<double> > &B) const;

		// Trace of product: sum_ij(This_ij X_ij)
		double traceProduct(const vector< vector<double> > &X) const;
		// Trace of low-rank product: sum_ij(This_ij (X^T Y)_ij)
		double traceProduct(const vector< vector<double> > &X,
							const vector< vector<double> > &Y) const;

		long autoClusteringColumns(vector<long> &best_block_ind,
									long min_nblock, long n) const;
		long clusteringColumns(vector<long> &block_ind, long nblock) const;
		long clusteringRows(vector<long> &block_ind, long nblock) const;
};

#endif

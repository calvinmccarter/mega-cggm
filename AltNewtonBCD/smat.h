#ifndef _SMAT_H_
#define _SMAT_H_

#include "util.h"

#include <vector>
#include <stdio.h>

using std::vector;

class smat_t {
	public:
		long p;
		long nnz;
		vector<double> values;
		vector<long> row_ptr;
		vector<long> col_idx;
		int is_symmetric;

        vector<double> diagonal;

		// Constructors
		smat_t();
		smat_t(long p_); 
		smat_t(smat_t &X);
		// Constructs from X + alpha*D:
		smat_t(smat_t &X, smat_t &D, double alpha);

		void setIdentity(long p_);
		void setTriplets(long p_, const vector<Triplet> &triplets_);
		void reset();
		void copyfrom(const smat_t &X);
		// Makes a symmetric matrix from lower triangular matrix:
		void symmetricfrom(const smat_t &X);
		void patternfrom(const smat_t &X);

		unsigned long isDiag() const;
		double L1Norm() const;
		double L1NormOffDiag() const;
		void print(FILE *fp) const;

		void ComputeAx(const vector<double> &x, vector<double> &Ax) const;
		void ComputeAx(const vector<double> &x, vector<double> &Ax, long p_) const;
		void ComputeAx_omp(
            const vector<double> &x, vector<double> &Ax, long p_) const;

		int ComputeAinvb(vector<double> &b, vector<double> &x, 
			double tol) const;
		int ComputeAinvb(vector<double> &b, vector<double> &x,
			long p_, double tol) const;
		int ComputeAinvb_pcg(vector<double> &b, vector<double> &x,
			double tol) const;

		int ComputeAinvb_omp(vector<double> &b, vector<double> &x,
			long p_, double tol) const;
		int ComputeAinvb_pcg_omp(vector<double> &b, vector<double> &x,
			long p_, double tol) const;

		int ComputeInv(vector<double> &x, long col, double tol) const;

		int ComputeLogdet_serial(double &result, double tol) const;
        int ComputeLogdet(double &result, double tol) const;

		// Compute (A+V*V')x
		void ComputeAVVtx(const vector< vector<double> > &V, 
			const vector<double> &x, vector<double> &AVVtx) const;

		// Solve (A+V*V')x = b
		int ComputeInvAVVt(const vector< vector<double> > &V, 
			long col, vector<double> &x, double tol) const;

		// Trace of product: sum_ij(This_ij X_ij)
		double traceProduct(const vector< vector<double> > &X) const;
		// Trace of low-rank product: sum_ij(This_ij (X^T Y)_ij)
		double traceProduct(const vector< vector<double> > &X,
							const vector< vector<double> > &Y) const;
	
		long clustering(vector<long> &block_ind, long nblock) const;
		long auto_clustering(vector<long> &block_ind,
							 long min_nblock, long n) const;
};

#endif

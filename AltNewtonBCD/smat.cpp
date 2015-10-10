#include "smat.h"

#include "util.h"

#if defined(LANG_M) || defined(MATLAB_MEX_FILE)
#include <mex.h>
#define PRINTF mexPrintf
#else
#include <stdio.h>
#define PRINTF printf
#endif

#include <algorithm>
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

double innerproduct(vector<double> &x, vector<double> &y ) {
	long p = x.size();
	double tmp=0;
	for ( long i=0 ; i<p ; i++ )
		tmp += x[i]*y[i];
	return tmp;
}

void vector_plus(vector<double> &x,
				 vector<double> &y,
				 vector<double> &z,
				 double c) {
	long p = x.size();
	for ( long i=0 ; i<p ; i++ ) {
		x[i] = y[i] + z[i]*c;
	}
}

smat_t::smat_t() {
	is_symmetric = 0;
}

smat_t::smat_t(long p_) {
	smat_t();
	p = p_;
	is_symmetric = 0;
	diagonal.resize(p_);
}

smat_t::smat_t(smat_t &X) {
	smat_t();
	p = X.p;
	nnz = X.nnz;
	values = X.values;
	row_ptr = X.row_ptr;
	col_idx = X.col_idx;
	is_symmetric = 0;
	diagonal = X.diagonal;
}

smat_t::smat_t(smat_t &X, smat_t &D, double alpha) {
	smat_t();
	is_symmetric = 0;
	p = X.p;
	nnz  = 0;
	values.resize(X.nnz+D.nnz);
	row_ptr.resize(p+1);
	col_idx.resize(X.nnz+D.nnz);
	for ( long i=0 ; i<p ; i++)	{
		row_ptr[i] = nnz;
		long idx1 = X.row_ptr[i];
		long idx2 = D.row_ptr[i];
		while ( (idx1 <X.row_ptr[i+1]) && (idx2<D.row_ptr[i+1])) {
			if ( X.col_idx[idx1] < D.col_idx[idx2] ) {
				col_idx[nnz] = X.col_idx[idx1];
				values[nnz] = X.values[idx1];
				idx1++, nnz++;
			} else if ( X.col_idx[idx1] > D.col_idx[idx2]) {
				col_idx[nnz] = D.col_idx[idx2];
				values[nnz] = alpha*D.values[idx2];
				idx2++, nnz++;
			} else {
				col_idx[nnz] = X.col_idx[idx1];
				values[nnz] = X.values[idx1]+alpha*D.values[idx2];
				idx1++, idx2++, nnz++;
			}
		}
		if ( idx1 < X.row_ptr[i+1] ) {
			for ( ; idx1<X.row_ptr[i+1] ; idx1++, nnz++ ) {
				col_idx[nnz] = X.col_idx[idx1];
				values[nnz] = X.values[idx1];
			}
		} else if ( idx2 < D.row_ptr[i+1]) {
			for ( ; idx2<D.row_ptr[i+1] ; idx2++, nnz++ ) {
				col_idx[nnz] = D.col_idx[idx2];
				values[nnz] = alpha*D.values[idx2];
			}
		}
	}
	row_ptr[p] = nnz;

	diagonal.resize(p);
	for (long i = 0; i < p; i++) {
		diagonal[i] = X.diagonal[i] + alpha*D.diagonal[i];
	}
}

void smat_t::setIdentity(long p_) {
	p = p_;
	nnz = p;
	values.resize(p);
	row_ptr.resize(p+1);
	col_idx.resize(p);
	for ( long i=0 ; i<p ; i++) {
		values[i] =1;
		col_idx[i] = i;
		row_ptr[i] = i;
	}
	row_ptr[p] = p;

	diagonal.resize(p);
	for (long i = 0; i < p; i++) {
		diagonal[i] = 1;
	}
}

void smat_t::setTriplets(long p_, const vector<Triplet> &triplets_) {
	p = p_;
	nnz = triplets_.size();
	is_symmetric = 0;
	vector<Triplet> triplets = triplets_;
	
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
		if (triplets[n].col > curr_row) {
			is_symmetric = 1;
		}
	}
	for (long idx = prev_row + 1; idx <= p; ++idx) {
		row_ptr[idx] = nnz;
	}
	
	diagonal.resize(p);
	for (long i = 0; i < p; i++) {
		diagonal[i] = values[row_ptr[i+1]-1];
	}
}

void smat_t::reset() {
	values.clear();
	col_idx.clear();
	row_ptr.clear();
	nnz = 0;
	diagonal.clear();
	diagonal.resize(p);
}

void smat_t::copyfrom(const smat_t &X) {
	p = X.p;
	nnz = X.nnz;
	values = X.values;
	row_ptr = X.row_ptr;
	col_idx = X.col_idx;
	diagonal = X.diagonal;
}

// make a symmetric matrix from lower triangular matrix
void smat_t::symmetricfrom(const smat_t &X) {
	is_symmetric = 1;
	p=X.p;
	diagonal = X.diagonal;

	nnz = 0;
	row_ptr.resize(p+1,0);
	col_idx.resize(2*X.nnz);
	values.resize(2*X.nnz);
	vector<long> tmp(p);
	for ( long i=0 ; i<p ; i++ ) {
		tmp[i] = X.row_ptr[i];
	}
	for (long i=0 ; i<p ; i++) {
		row_ptr[i] = nnz;
		for ( long idx = X.row_ptr[i] ; idx<X.row_ptr[i+1] ; idx++ ) {
			col_idx[nnz] = X.col_idx[idx];
			values[nnz] = X.values[idx];
			nnz++;
		}
		for ( long j=i+1 ; j<p ; j++ ) {
			if ( tmp[j] < X.row_ptr[j+1] ) {
				if ( X.col_idx[tmp[j]] == i) {
					col_idx[nnz] = j;

					values[nnz] = X.values[tmp[j]];

					nnz++;
					tmp[j]++;
				}
			}
		}
	}
	row_ptr[p] = nnz;
}

void smat_t::patternfrom(const smat_t &X) {
	p = X.p;
	nnz = X.nnz;
	values.clear();
	values.resize(nnz);
	row_ptr = X.row_ptr;
	col_idx = X.col_idx;
	diagonal.clear();
	diagonal.resize(nnz);
}

unsigned long smat_t::isDiag() const {
	for (long i = 0; i < p; i++) {
		for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
			if (col_idx[idx] != i) {
				return 0;
			}
		}
	}
	return 1;
}

double smat_t::L1Norm() const {
	double result = 0;
	if (is_symmetric != 1) {
		for (long i = 0; i < p; i++) {
			for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
				if (col_idx[idx] == i) {
					result += fabs(values[idx]);
				} else {
					result += 2 * fabs(values[idx]);
				}
			}
		}
	} else {
		for (long i = 0; i < p; i++) {
			for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
				result += fabs(values[idx]);
			}
		}
	}
	return result;
}

double smat_t::L1NormOffDiag() const {
	double result = 0;
	if (is_symmetric != 1) {
		for (long i = 0; i < p; i++) {
			for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
				if (col_idx[idx] != i) {
					result += 2 * fabs(values[idx]);
				}
			}
		}
	} else {
		for (long i = 0; i < p; i++) {
			for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
				if (col_idx[idx] != i) {
					result += fabs(values[idx]);
				}
			}
		}
	}
	return result;
}

void smat_t::print(FILE *fp) const {
	//fprintf(fp, "p: %ld, nnz: %ld\n", p, nnz);
	PRINTF("p: %ld, nnz: %ld\n", p, nnz);
	for ( long i=0 ; i<p ; i++ ) {
		for ( long idx = row_ptr[i] ; idx<row_ptr[i+1] ; idx++ ) {
			//fprintf(fp, "%ld %ld %.10lf\n", 
			//	i+1, col_idx[idx]+1, values[idx]);
			PRINTF("%ld %ld %.10lf\n", 
				i+1, col_idx[idx]+1, values[idx]);

		}
	}
}

void smat_t::ComputeAx(const vector<double> &x, vector<double> &Ax) const {
	ComputeAx(x, Ax, p);
}

void smat_t::ComputeAx(const vector<double> &x, vector<double> &Ax, long p_) const {
	Ax.resize(p_);
	for ( long i=0 ; i<p_ ; i++ ) {
		Ax[i] = 0;
	}
	if ( is_symmetric != 1 ) {
		for (long i=0 ; i<p_ ; i++ ) {
			for ( long idx = row_ptr[i] ; idx < row_ptr[i+1] ; idx++) {
				long j = col_idx[idx];
				double v = values[idx];
				Ax[i] += v*x[j];
				if (i != j) {
					Ax[j] += v*x[i];
				}
			}
		}
	} else {
		for (long i=0 ; i<p_ ; i++ ) {
			double tmp=0;
			for ( long idx = row_ptr[i] ; idx < row_ptr[i+1] ; idx++) {
				long j = col_idx[idx];
				double v = values[idx];
				tmp += v*x[j];
			}
			Ax[i] = tmp;
		}
	}
}

void smat_t::ComputeAx_omp(
		const vector<double> &x, vector<double> &Ax, long p_) const {
	Ax.resize(p_);
	for ( long i=0 ; i<p_ ; i++ ) {
		Ax[i] = 0;
	}
	for (long i=0 ; i<p_ ; i++ ) {
		for ( long idx = row_ptr[i] ; idx < row_ptr[i+1] ; idx++) {
			long j = col_idx[idx];
			double v = values[idx];
			Ax[i] += v*x[j];
			if (i != j) {
				Ax[j] += v*x[i];
			}
		}
	}
}

int smat_t::ComputeInv(vector<double> &x, long col, double tol) const {
	vector<double> ei(p, 0);
	ei[col] = 1;
	return ComputeAinvb_pcg(ei, x, tol);
}	

int smat_t::ComputeAinvb(vector<double> &b,
						 vector<double> &x, 
						 double tol) const {
	return ComputeAinvb(b, x, p, tol);
}

int smat_t::ComputeAinvb(vector<double> &b,
						 vector<double> &x,
						 long p_,
						 double tol) const {
	vector<double>  r(p_), pp(p_), Ax_result(p_);
	int maxiter = 30;

	// Initial from x = 0
	for (long i = 0; i < p_; i++) {
		x[i] = 0;
	}
	r = b;
	pp = r;
	double r_norm = innerproduct(r,r);
	double initial = r_norm;
//			printf("initial rnorm: %lf\n", r_norm);
	if (r_norm < 1e-13) {
		return 0;
	}
	int iter;
	for (iter = 0; iter < maxiter; iter++) {
		ComputeAx(pp, Ax_result, p_);
		double alpha = innerproduct(r, r)/innerproduct(Ax_result, pp);
		vector_plus(x, x, pp, alpha);
		vector_plus(r, r, Ax_result, (-1)*alpha);
		double r_norm_now = innerproduct(r,r);
//				printf("residual: %lf\n", r_norm_now);
		if (r_norm_now < tol*initial) {
			break;
		}
		double beta = r_norm_now/r_norm;
		r_norm = r_norm_now;
		vector_plus(pp, r, pp, beta);
	}
	return (iter+1);
}

int smat_t::ComputeAinvb_pcg(
		vector<double> &b,
		vector<double> &x,
		double tol) const {
	vector<double> r(p), z(p), P(p), AP(p);
	int maxiter = 30;
	
	// Compute inv(M)
	vector<double> Minv(p);
	for (long i = 0; i < p; i++) {
		Minv[i] = 1.0 / diagonal[i];
		if (!isnormal(Minv[i])) {
			PRINTF("Minv failed\n");
			return 1;
		}
	}

	// Initialize from x = 0
	for (long i = 0; i < p; i++) {
		x[i] = 0;
	}
	r = b;
	for (long i = 0; i < p; i++) {
		z[i] = Minv[i] * r[i];
	}
	P = z;
	double r_norm_0 = innerproduct(r, r);
	if (r_norm_0 < 1e-13) {
		return 0;
	}
	int iter;
	for (iter = 0; iter < maxiter; iter++) {
		ComputeAx(P, AP, p);
		double rtz = innerproduct(r, z);
		double alpha = rtz/innerproduct(P,AP);

		vector_plus(x, x, P, alpha);
		vector_plus(r, r, AP, (-1)*alpha);
		double r_norm = innerproduct(r, r);
		if (r_norm < tol*r_norm_0) {
			break;
		}
		for (long i = 0; i < p; i++) {
			z[i] = Minv[i] * r[i];
		}
		double beta = innerproduct(r, z)/rtz;
		vector_plus(P, z, P, beta);
	}
	return (iter+1);
}	

int smat_t::ComputeAinvb_omp(
		vector<double> &b,
		vector<double> &x,
		long p_, double tol) const {
	vector<double>  r(p_), pp(p_), Ax_result(p_);
	int maxiter = 30;

	// Initial from x = 0
	for (int i=0; i < p_; i++) {
		x[i] = 0;
	}
	r = b;
	pp = r;
	double r_norm = innerproduct(r,r);
	double initial = r_norm;
	if (r_norm < 1e-13) {
		return 0;
	}
	int iter;
	for (iter = 0; iter < maxiter; iter++) {
		ComputeAx_omp(pp, Ax_result, p_);
		double alpha = innerproduct(r, r)/innerproduct(Ax_result, pp);
		vector_plus(x, x, pp, alpha);
		vector_plus(r, r, Ax_result, (-1)*alpha);
		double r_norm_now = innerproduct(r,r);
		if (r_norm_now < tol*initial) {
			break;
		}
		double beta = r_norm_now/r_norm;
		r_norm = r_norm_now;
		vector_plus(pp, r, pp, beta);
	}
	return (iter+1);
}

int smat_t::ComputeAinvb_pcg_omp(
		vector<double> &b,
		vector<double> &x,
		long p_, double tol) const {
	vector<double>  r(p_), z(p_), P(p_), AP(p_);
	int maxiter = 30;

	// Compute inv(M)
	vector<double> Minv(p);
	for (long i = 0; i < p; i++) {
		Minv[i] = 1.0 / diagonal[i];
		if (!isnormal(Minv[i])) {
			PRINTF("Minv failed\n");
			return 1;
		}
	}

	// Initialize from x = 0
	for (long i = 0; i < p_; i++) {
		x[i] = 0;
	}
	r = b;
	for (long i = 0; i < p_; i++) {
		z[i] = Minv[i] * r[i];
	}
	P = z;
	double r_norm_0 = innerproduct(r, r);
	if (r_norm_0 < 1e-13) {
		return 0;
	}
	int iter;
	for (iter = 0; iter < maxiter; iter++) {
		ComputeAx_omp(P, AP, p_);
		double rtz = innerproduct(r, z);
		double alpha = rtz/innerproduct(P, AP);
		vector_plus(x, x, P, alpha);
		vector_plus(r, r, AP, (-1)*alpha);
		double r_norm = innerproduct(r,r);
		if (r_norm < tol*r_norm_0) {
			break;
		}
		for (long i = 0; i < p_; i++) {
			z[i] = Minv[i] * r[i];
		}
		double beta = innerproduct(r, z)/rtz;
		vector_plus(P, z, P, beta);
	}
	return (iter+1);
}
int smat_t::ComputeLogdet_serial(double &result, double tol) const {
	double logdet = 0;
	// check if all diagonals are positive
	int pd = 1;
	long i;
	for ( i=0 ; i<p ; i++ ) {
		if ( row_ptr[i+1]-1 < 0) {
			break;
		}
		if ( col_idx[row_ptr[i+1]-1] != i) {
			break;
		}
		if ( values[row_ptr[i+1]-1] < 0) {
			break;
		}
	}
	if ( i<p ) {
		return 0;
	}

	logdet = log(values[0]);
	for ( long pend=1 ; pend<p ; pend++) {
		// compute
		double tmp = values[row_ptr[pend+1]-1];
		vector<double> b (pend,0);
		for ( long idx = row_ptr[pend] ; idx < row_ptr[pend+1]-1 ; idx++ ) {
			b[col_idx[idx]] = values[idx];
		}
		if ( pend == 1 ) {
			tmp = tmp - b[0]*b[0]/values[0];
		} else {
			vector<double> Ainvb (pend,0);
			ComputeAinvb(b, Ainvb, pend, tol);
			tmp = tmp - innerproduct(b, Ainvb);
		}

		if ( tmp <= 0) {
			return 0;
		}
		logdet += log(tmp);
//				if ( pend %100 == 0)
//					printf("logdet iter %d: %lf\n", pend, logdet);
	}
	result = logdet;
	return 1;
}

int smat_t::ComputeLogdet(double &result, double tol) const {
	double logdet = 0;
	// check if all diagonals are positive
	int pd = 1;
	long i;
	for (i = 0; i < p; i++) {
		if (row_ptr[i+1]-1 < 0) {
			break;
		}
		if (col_idx[row_ptr[i+1]-1] != i) {
			break;
		}
		if (values[row_ptr[i+1]-1] < 0) {
			break;
		}
	}
	if (i < p) {
		return 0;
	}

	logdet = 0;
	int errorflag = 0;
#pragma omp parallel for shared(errorflag) reduction(+:logdet) schedule(dynamic)
	for (long pend = 1; pend < p; pend++) {
		if (errorflag == 1) {
			continue;
		}
		// compute
		double tmp = values[row_ptr[pend+1]-1];
		vector<double> b(pend, 0);
		for (long idx = row_ptr[pend]; idx < row_ptr[pend+1]-1; idx++) {
			b[col_idx[idx]] = values[idx];
		}

		if (pend == 1) {
			tmp = tmp - b[0]*b[0]/values[0];
		} else {
			vector<double> Ainvb (pend,0);
			int numiter = ComputeAinvb_pcg_omp(b, Ainvb, pend, tol);
			tmp = tmp - innerproduct(b, Ainvb);
		}

		if (tmp <= tol) {
			errorflag = 1;
		}
		logdet += log(tmp);
	}
	logdet += log(values[0]);
	result = logdet;
	if (errorflag == 1) {
		return 0;
	} else {
		return 1;
	}
}

int smat_t::ComputeAVVtx(const vector< vector<double> > &V, 
			const vector<double> &x, vector<double> &AVVtx) const {
	long r = V[0].size();
	vector<double> Ax_result(p, 0.0);
	vector<double> Vtx_result(r, 0.0);
	vector<double> VVtx_result(p, 0.0);
	AVVtx.resize(p);

	ComputeAx(x, Ax_result);
	
	for (long i = 0; i < r; i++) {
		for (long j = 0; j < p; j++) {
			Vtx_result[i] += V[j][i] * x[j];
		}
	}
	for (long i = 0; i < p; i++) {
		for (long j = 0; j < r; j++) {
			VVtx_result[i] += V[i][j] * Vtx_result[j];
		}
	}
	for (long i = 0; i < p; i++) {
		AVVtx[i] = Ax_result[i] + VVtx_result[i];
	}
}


// Solve (A+V*V')x = b
// V is tall matrix in row-major format
int smat_t::ComputeInvAVVt(const vector< vector<double> > &V,
	long col, vector<double> &x, double tol) const {

	vector<double> b(p, 0);
	b[col] = 1;

	long z_dims = V[0].size();
	vector<double>  r(p), pp(p);
	vector<double> Ax_result(p), AVVtx_result(p);
	int maxiter = 50;

	// Initial from x = 0
	for (long i = 0; i < p; i++) {
		x[i] = 0;
	}
	r = b;
	pp = r;
	double r_norm = innerproduct(r,r);
	double initial = r_norm;
			//printf("initial rnorm: %lf\n", r_norm);
	if (r_norm < 1e-13) {
		return 0;
	}
	int iter;
	for (iter = 0; iter < maxiter; iter++) {
		ComputeAVVtx(V, pp, AVVtx_result);

		double alpha = innerproduct(r, r)/innerproduct(AVVtx_result, pp);
		vector_plus(x, x, pp, alpha);
		vector_plus(r, r, AVVtx_result, (-1)*alpha);
		double r_norm_now = innerproduct(r,r);
				//printf("residual: %lf\n", r_norm_now);
		if (r_norm_now < tol*initial) {
			break;
		}
		double beta = r_norm_now/r_norm;
		r_norm = r_norm_now;
		vector_plus(pp, r, pp, beta);
	}
	return (iter+1);
}

double smat_t::traceProduct(const vector< vector<double> > &X) const {
	double result = 0;
	if (is_symmetric != 1) {
		for (long i = 0; i < p; i++) {
			for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
				long j = col_idx[idx];
				if (i == j) {
					result += values[idx] * X[i][j];
				} else {
					result += 2 * values[idx] * X[i][j];
				}
			}
		}
	} else {
		for (long i = 0; i < p; i++) {
			for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
				long j = col_idx[idx];
				result += values[idx] * X[i][j];
			}
		}
	}	
	return result;
}

// Both X and Y are column-major
double smat_t::traceProduct(const vector< vector<double> > &X,
					const vector< vector<double> > &Y) const {
	long n = X[0].size();
	double result = 0;
	if (is_symmetric != 1) {
		#pragma omp parallel for reduction(+:result) schedule(static)
		for (long i = 0; i < p; i++) {
			for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
				long j = col_idx[idx];
				double xy_ij = 0;
				for (long nn = 0; nn < n; nn++) {
					xy_ij += X[i][nn] * Y[j][nn];
				}
				if (i == j) {
					result += values[idx] * xy_ij;
				} else {
					result += 2 * values[idx] * xy_ij;
				}
			}
		}
	} else {
		#pragma omp parallel for reduction(+:result) schedule(static)
		for (long i = 0; i < p; i++) {
			for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
				long j = col_idx[idx];
				double xy_ij = 0;
				for (long nn = 0; nn < n; nn++) {
					xy_ij += X[i][nn] * Y[j][nn];
				}
				result += values[idx] * xy_ij;
			}
		}
	}
	return result;
}

long smat_t::clustering(vector<long> &block_ind, long nblock) const {
	block_ind.resize(p);
	for (long i = 0; i < p; ++i) {
		block_ind[i] = 0;
	}
	if (nblock == 1) {
		return 0;
	}
	if (is_symmetric == 1) {
		PRINTF("Error: clustering symmetric format smat_t\n");
		return -1;
	}
	smat_t Lambda_sym;
	Lambda_sym.symmetricfrom(*this);	

	// Metis inputs
	idx_t nvtxs = p;  // number of vertices
	idx_t ncon = 1;   // number of balancing constraints
	idx_t *xadj = (idx_t *)malloc(sizeof(idx_t)*(Lambda_sym.p+1));
	for (long i = 0; i <= Lambda_sym.p ;i++) {
		xadj[i] = Lambda_sym.row_ptr[i];
	}
	idx_t *adjncy = (idx_t *)malloc(sizeof(idx_t)*Lambda_sym.nnz);
	for (long idx = 0; idx < Lambda_sym.nnz; idx++) {
		adjncy[idx] = Lambda_sym.col_idx[idx];
	}
 	idx_t nparts = nblock;
	idx_t options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // edgecut minimization
	//options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; // rand matching coarsening
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
	idx_t objval; // stores edge-cut of solution
	idx_t *part = (idx_t *)malloc(sizeof(idx_t)*p);

	// Run graph partitioning
	int flag = METIS_PartGraphKway(
		&nvtxs, &ncon, xadj, adjncy,
		NULL, NULL, NULL, &nparts, NULL,
		NULL, options, &objval, part);
	if (flag != METIS_OK) {
		PRINTF("Metis error %d \n", flag);
	}
	for (long i = 0; i < p; i++) {
		block_ind[i] = part[i];
	}
	free(xadj);
	free(adjncy);
	free(part);
	
	long num_boundary = 0;
	vector< vector<long> > block_list;
	get_block_list(block_ind, block_list);
	for (long block_z = 0; block_z < block_list.size(); block_z++) {
		vector<bool> is_boundary(p, false);
		long block_z_size = block_list[block_z].size();
		for (long i_ix = 0; i_ix < block_z_size; i_ix++) {
			long i = block_list[block_z][i_ix];
			for (long idx = row_ptr[i]; idx < row_ptr[i+1]; idx++) {
				long j = col_idx[idx];
				if (block_ind[j] != block_z) {
					if (!is_boundary[j]) {
						num_boundary++;
					}
					is_boundary[j] = true;
				}
			}
		}
	}
	return num_boundary;
}

// Returns total number of boundary variables
long smat_t::auto_clustering(
		vector<long> &best_block_ind, long min_nblock, long n) const {
	long max_iters = 7;
	long pcg = 20;
	double best_cost = numeric_limits<double>::max();
	long best_boundary;

	long curr_k = min_nblock;
	for (long i = 0; i < max_iters; i++) {
		if (curr_k > p) {
			break;
		}
		vector<long> curr_block_ind;
		long curr_boundary = clustering(curr_block_ind, curr_k);
		double curr_cost = 
			curr_k*curr_boundary*1.0e-2*(pcg*nnz + n*p) + 2*nnz*p/(double(curr_k));
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

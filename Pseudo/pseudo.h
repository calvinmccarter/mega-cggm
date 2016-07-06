#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double,Eigen::ColMajor,long int> SpMatrixC;
typedef Eigen::SparseMatrix<double,Eigen::RowMajor,long int> SpMatrixR;
typedef SpMatrixC::InnerIterator InIter;
typedef SpMatrixR::InnerIterator InIterR;
typedef Eigen::Triplet<double> Triplet;

struct PseudoOptions {
		PseudoOptions() : 
	quiet(false),
	max_iters(50),
    tol(1.0e-2),
    diag_dominant(true) {}
	
	// Whether to print info messages
	bool quiet;

	// Maximum Lasso iterations 
	int max_iters;

	// Tolerance parameter for Lasso
	double tol;

	// Make Lambda diagonal dominant
	bool diag_dominant; 

	// Use SAFE screening rule
	bool screen;
};

struct PseudoStats {
	std::vector<double> time;
};

void Pseudo(
	const Eigen::MatrixXd& Y,
	const Eigen::MatrixXd& X,
	double lambda_y,
	double lambda_x,
	PseudoOptions& options,
	SpMatrixC& Lambda,
	SpMatrixC& Theta,
	PseudoStats* stats);


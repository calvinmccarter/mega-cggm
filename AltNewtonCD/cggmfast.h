#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double,Eigen::ColMajor,long int> SpMatrixC;
typedef Eigen::SparseMatrix<double,Eigen::RowMajor,long int> SpMatrixR;
typedef SpMatrixC::InnerIterator InIter;
typedef SpMatrixR::InnerIterator InIterR;
typedef Eigen::Triplet<double> Triplet;

struct CGGMOptions {
		CGGMOptions() : 
	quiet(false),
	max_outer_iters(50),
	max_Theta_iters(1),
	max_Lambda_iters(1),
    max_ls_iters(10),
    alpha(1),
    beta(0.5),
    sigma(1.0e-4),
    tol(1.0e-2),
    cd_tol(0.05),
    refit(false) {}
	
	// Whether to print info messages
	bool quiet;

	// Maximum iterations 
	int max_outer_iters;
	int max_Theta_iters;
	int max_Lambda_iters;
	int max_ls_iters;

	// Line search parameters
	double alpha;
	double beta;
	double sigma;

	// Tolerance parameters 
	double tol;
	double cd_tol;

	// Algorithm parameters
	bool refit; // if true: only remove edges when fitting model
	// Only makes sense when providing Lambda0 and Theta0
};

struct CGGMStats {
	std::vector<double> objval;
	std::vector<double> time;
	std::vector<double> active_set_size;
	std::vector<double> active_theta;
	std::vector<double> active_lambda;
    std::vector<double> subgrad;
	std::vector<double> l1norm;

	std::vector<double> time_lambda_active;
	std::vector<double> time_theta_active;
	std::vector<double> time_theta_cd;
	std::vector<double> time_theta_cd_cd;
	std::vector<double> time_theta_cd_qr;
	std::vector<double> time_lambda_cd;
	std::vector<double> time_lambda_cd_cd;
	std::vector<double> time_lambda_cd_ls;
};

void CGGMfast(
	const Eigen::MatrixXd& Y,
	const Eigen::MatrixXd& X,
	double lambda_y,
	double lambda_x,
	CGGMOptions& options,
	SpMatrixC& Lambda,
	SpMatrixC& Theta,
	CGGMStats* stats);


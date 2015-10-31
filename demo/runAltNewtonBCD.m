function [Lambda, Theta, stats] = runAltNewtonCD(...
    Y, X, lambdaLambda, lambdaTheta, options)
% Inputs:
% Y (n samples x q dimensions)
% X (n samples x p dimensions)
% lambdaLambda (regularization for Lambda)
% lambdaTheta (regularization for Theta)
% [options] struct with the following options:
%   - Lambda0(none): q x q sparse matrix to initialize Lambda
%   - Theta0(none): p x q sparse matrix to initialize Theta
%   - verbose(0): show information or not (0 or 1)
%   - max_outer_iters(50): max number of outer iterations
%   - sigma(1e-4): backtracking termination criterion
%   - tol(1e-2): subgradient/L1norm tolerance for terminating outer loop
%   - obj_tol(1.0e-13): CG tolerance for calculating objective function
%   - grad_tol(1.0e-10): CG tolerance for calculating gradient
%   - hess_tol(1.0e-8): CG tolerance for calculating hessian
%   - num_blocks_Lambda(-1): number of blocks for Lambda BCD
%   - num_blocks_Theta(-1): number of blocks for Theta BCD
%   - memory_usage(32000): memory available to process in Mb
%   - max_threads(4): max number of threads to use 
%   - refit(0): refit selected model without adding any edges

    olddir = pwd;
    thisfunc = which(mfilename());
    thisdir = thisfunc(1:end-17);
    cd(thisdir);
  
    verbose = 0;
    max_outer_iters = 50;
    sigma = 1.0e-4;
    tol = 1.0e-2;
    obj_tol = 1.0e-13;
    grad_tol = 1.0e-10;
    hess_tol = 1.0e-8;
    num_blocks_Lambda = -1; % <0: min with memory, 0: auto, >0: user input
    num_blocks_Theta = -1; % <0: min with memory, 0: auto, >0: user input
    memory_usage = 32000; % default is 32 Gb
    max_threads = 4; % actual num threads is min(max_threads, omp_get_max_threads())
    refit = 0;
    L0_str = '';
    T0_str = '';
    dummy = randi(1e6);
    if exist('options', 'var')
        if isfield(options, 'verbose')
            verbose = options.verbose;
        end
        if isfield(options, 'max_outer_iters')
            max_outer_iters = options.max_outer_iters;
        end
        if isfield(options, 'sigma')
            refit = options.sigma;
        end
        if isfield(options, 'tol')
            tol = options.tol;
        end
        if isfield(options, 'obj_tol')
            obj_tol = options.obj_tol;
        end
        if isfield(options, 'grad_tol')
            grad_tol = options.grad_tol;
        end
        if isfield(options, 'hess_tol')
            hess_tol = options.hess_tol;
        end
        if isfield(options, 'num_blocks_Lambda')
            num_blocks_Lambda = options.num_blocks_Lambda;
        end
        if isfield(options, 'num_blocks_Theta')
            num_blocks_Theta = options.num_blocks_Theta;
        end
        if isfield(options, 'memory_usage')
            memory_usage = options.memory_usage;
        end
        if isfield(options, 'max_threads')
            max_threads = options.max_threads;
        end
        if isfield(options, 'Lambda0')
            Lambda0file = sprintf('Lambda0-dummy-%i.txt', dummy);
            sparse_to_txt(Lambda0file, options.Lambda0);
            L0_str = sprintf('-L \"%s\" ', Lambda0file);
        end
        if isfield(options, 'Theta0')
            Theta0file = sprintf('Theta0-dummy-%i.txt', dummy);
            sparse_to_txt(Theta0file, options.Theta0);
            T0_str = sprintf('-T \"%s\" ', Theta0file);
        end
        if isfield(options, 'refit')
            refit = options.refit;
        end
    end
    [n_y, q] = size(Y);
    [n_x, p] = size(X);

    Yfile = sprintf('Y-dummy-%i.txt', dummy);
    Xfile = sprintf('X-dummy-%i.txt', dummy);
    Lambdafile = sprintf('Lambda-dummy-%i.txt', dummy);
    Thetafile = sprintf('Theta-dummy-%i.txt', dummy);
    statsfile = sprintf('stats-dummy-%i.txt', dummy);
    dlmwrite(Yfile, Y, 'delimiter', ' ', 'precision', 10);
    dlmwrite(Xfile, X, 'delimiter', ' ', 'precision', 10);
    
    basic_str = sprintf('-y %g -x %g -v %i -i %i -s %g -q %g ', ...
        lambdaLambda, lambdaTheta, ...
        verbose, max_outer_iters, sigma, tol);
    extra_str = sprintf('-o %g -g %g -h %g -l %i -t %i -m %i -n %i -r %i ', ...
        obj_tol, grad_tol, hess_tol, num_blocks_Lambda, num_blocks_Theta, ...
        memory_usage, max_threads, refit);
    command_str = sprintf('%s %s %s %s %s %i %i %i %i %s %s %s %s %s', ...
        '../AltNewtonBCD/hugecggm_run', basic_str, extra_str, ...
        L0_str, T0_str, ...
        n_x, p, n_y, q, Yfile, Xfile, ...
        Lambdafile, Thetafile, statsfile);

    system(command_str);
    Lambda = txt_to_sparse(Lambdafile);
    Theta = txt_to_sparse(Thetafile);
    stats = txt_to_struct(statsfile);
    system(['rm ' Yfile ' ' Xfile ' ' Lambdafile ' ' Thetafile ' ' statsfile]);
    if L0_str
        system(['rm ' Lambda0file]);
    end
    if T0_str
        system(['rm ' Theta0file]);
    end
    cd(olddir);
end

function [Lambda, Theta, stats] = runAltNewtonCD(...
    Y, X, lambdaLambda, lambdaTheta, options)
% Inputs:
% Y (n samples x q dimensions)
% X (n samples x p dimensions)
% lambdaLambda (regularization for Lambda)
% lambdaTheta (regularization for Theta)
% [options] struct with the following options:
%   - verbose(0): show information or not (0 or 1)
%   - max_outer_iters(50): max number of outer iterations
%   - sigma(1e-4): backtracking termination criterion
%   - tol(1e-2): tolerance for terminating outer loop
%   - refit(0): refit selected model without penalty
    
    olddir = pwd;
    thisfunc = which(mfilename());
    thisdir = thisfunc(1:end-16);
    cd(thisdir);

    verbose = 0;
    max_outer_iters = 50;
    sigma = 1.0e-4;
    tol = 1.0e-2;
    refit = 0;

    if exist('options', 'var')
        if isfield(options, 'verbose')
            verbose = options.verbose;
        end
        if isfield(options, 'max_outer_iters')
            max_outer_iters = options.max_outer_iters;
        end
        if isfield(options, 'tol')
            tol = options.tol;
        end
        if isfield(options, 'refit')
            refit = options.refit;
        end

    end
    [n_y, q] = size(Y);
    [n_x, p] = size(X);

    dummy = randi(1e6);
    Yfile = sprintf('Y-dummy-%i.txt', dummy);
    Xfile = sprintf('X-dummy-%i.txt', dummy);
    Lambdafile = sprintf('Lambda-dummy-%i.txt', dummy);
    Thetafile = sprintf('Theta-dummy-%i.txt', dummy);
    statsfile = sprintf('stats-dummy-%i.txt', dummy);
    dlmwrite(Yfile, Y, 'delimiter', ' ', 'precision', 10);
    dlmwrite(Xfile, X, 'delimiter', ' ', 'precision', 10);
    option_str = sprintf('-y %f -x %f -v %i -i %i -s %f -q %f -r %i', ...
        lambdaLambda, lambdaTheta, ...
        verbose, max_outer_iters, sigma, tol, refit);
    command_str = sprintf('%s %s %i %i %i %i %s %s %s %s %s', ...
        '../AltNewtonCD/cggmfast_run', option_str, ...
        n_x, p, n_y, q, Yfile, Xfile, ...
        Lambdafile, Thetafile, statsfile);

    system(command_str);
    Lambda = txt_to_sparse(Lambdafile);
    Theta = txt_to_sparse(Thetafile);
    stats = txt_to_struct(statsfile);
    system(['rm ' Yfile ' ' Xfile ' ' Lambdafile ' ' Thetafile ' ' statsfile]);
    cd(olddir);
end

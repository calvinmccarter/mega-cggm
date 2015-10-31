function [ll] = loglik(Y, X, Lambda, Theta)
% returns loglikelihood(Y, X; Lambda, Theta), useful for cross-validation
    [n, p] = size(X);
    [n, q] = size(Y);
    Y_XB = Y - (Lambda \ (-X*Theta)')'; % (n, q)
    Y_XB_LambdaT = (Y_XB * Lambda)';
    ll = -sum(Y_XB_LambdaT(:) .* Y_XB(:)) + n*logdet(Lambda);
end

function [ld] = logdet(A)
    spA = sparse(A);
    [Lq, p] = chol(spA, 'lower');
    if p ~= 0
        warning('Taking determinant of non-PD matrix');
        ld = inf;
        return;
    end 
    ld = 2*sum(log(full(diag(Lq))));
end

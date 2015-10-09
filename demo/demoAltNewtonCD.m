system('(cd ../AltNewtonCD && make)');
close all;
X = dlmread('Xfile'); % (10 samples x 15 features)
Y = dlmread('Yfile'); % (10 samples x 12 features)
lambdaLambda = 0.1;
lambdaTheta = 0.2;

% vanilla example:
[Lambda_vanilla, Theta_vanilla, stats_vanilla] = runAltNewtonCD(...
    Y, X, lambdaLambda, lambdaTheta);
figure('name', 'vanilla');
subplot(1,2,1); imagesc(Lambda_vanilla);
subplot(1,2,2); imagesc(Theta_vanilla);
assert(abs(-10.74786 - stats_vanilla.objval(end)) < 1.0e-2);

% options example:
options.tol = 1.0e-5;
[Lambda_exact, Theta_exact, stats_exact] = runAltNewtonCD(...
    Y, X, lambdaLambda, lambdaTheta, options);
figure('name', 'high-accuracy');
subplot(1,2,1); imagesc(Lambda_exact);
subplot(1,2,2); imagesc(Theta_exact);

% semi-supervised example:
partialX = X(1:end-4,:); % Y has 10 samples, X has 6 samples
[Lambda_semi, Theta_semi, stats_semi] = runAltNewtonCD(...
    Y, partialX, lambdaLambda, lambdaTheta);
figure('name', 'missing X');
subplot(1,2,1); imagesc(Lambda_semi);
subplot(1,2,2); imagesc(Theta_semi);


[r,s] = system('(cd ../AltNewtonBCD && make)');
if r ~= 0
    fprintf('first go to AltNewtonBCD directory and run $ make\n');
    return;
end
close all;
X = dlmread('Xfile'); % (10 samples x 15 features)
Y = dlmread('Yfile'); % (10 samples x 12 features)
lambdaLambda = 0.1;
lambdaTheta = 0.2;

% vanilla example:
[Lambda_vanilla, Theta_vanilla, stats_vanilla] = runAltNewtonBCD(...
    Y, X, lambdaLambda, lambdaTheta);
figure('name', 'vanilla');
subplot(1,2,1); imagesc(Lambda_vanilla);
subplot(1,2,2); imagesc(Theta_vanilla);
assert(abs(-10.8671 - stats_vanilla.objval(end)) < 1.0e-2);

% options example:
options.tol = 1.0e-5;
[Lambda_exact, Theta_exact, stats_exact] = runAltNewtonBCD(...
    Y, X, lambdaLambda, lambdaTheta, options);
figure('name', 'high-accuracy');
subplot(1,2,1); imagesc(Lambda_exact);
subplot(1,2,2); imagesc(Theta_exact);

% initialized example
options.Lambda0 = Lambda_vanilla;
options.Theta0 = Theta_vanilla;
[Lambda_warm, Theta_warm, stats_warm] = runAltNewtonBCD(...
    Y, X, lambdaLambda, lambdaTheta, options);
figure('name', 'warm started');
subplot(1,2,1); imagesc(Lambda_warm);
subplot(1,2,2); imagesc(Theta_warm);


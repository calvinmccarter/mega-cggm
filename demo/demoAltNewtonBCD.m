[r,s] = system('(cd ../AltNewtonBCD && make)');
if r ~= 0
    fprintf('first go to AltNewtonBCD directory and run $ make\n');
    return;
end
close all; clear all;
X = dlmread('Xfile'); % (10 samples x 15 features)
Y = dlmread('Yfile'); % (10 samples x 12 features)
lambdaLambda = 0.1;
lambdaTheta = 0.2;

% vanilla example:
[Lambda_vanilla, Theta_vanilla, stats_vanilla] = runAltNewtonBCD(...
    Y, X, lambdaLambda, lambdaTheta);
figure('name', 'vanilla');
subplot(1,2,1); imagesc(Lambda_vanilla); colorbar;
title('Lambda'); set(gca,'dataAspectRatio',[1 1 1]);
subplot(1,2,2); imagesc(Theta_vanilla); colorbar;
title('Theta'); set(gca,'dataAspectRatio',[1 1 1]);

% options example:
options.tol = 1.0e-5;
[Lambda_exact, Theta_exact, stats_exact] = runAltNewtonBCD(...
    Y, X, lambdaLambda, lambdaTheta, options);
figure('name', 'high-accuracy');
subplot(1,2,1); imagesc(Lambda_exact);
title('Lambda'); set(gca,'dataAspectRatio',[1 1 1]);
subplot(1,2,2); imagesc(Theta_exact);
title('Theta'); set(gca,'dataAspectRatio',[1 1 1]);

% warm-starting example
options.Lambda0 = Lambda_vanilla;
options.Theta0 = Theta_vanilla;
[Lambda_warm, Theta_warm, stats_warm] = runAltNewtonBCD(...
    Y, X, 0.2*lambdaLambda, 0.5*lambdaTheta, options);
figure('name', 'warm started');
subplot(1,2,1); imagesc(Lambda_warm);
title('Lambda'); set(gca,'dataAspectRatio',[1 1 1]);
subplot(1,2,2); imagesc(Theta_warm);
title('Theta'); set(gca,'dataAspectRatio',[1 1 1]);

% refitting example
options.Lambda0 = Lambda_vanilla;
options.Theta0 = Theta_vanilla;
options.refit = 1; % don't add any edges
[Lambda_refit, Theta_refit, stats_refit] = runAltNewtonBCD(...
    Y, X, 0.2*lambdaLambda, lambdaTheta, options);
figure('name', 'refitting');
subplot(1,2,1); imagesc(Lambda_refit); colorbar;
title('Lambda'); set(gca,'dataAspectRatio',[1 1 1]);
subplot(1,2,2); imagesc(Theta_refit); colorbar;
title('Theta'); set(gca,'dataAspectRatio',[1 1 1]);

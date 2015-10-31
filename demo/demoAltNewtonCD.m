system('(cd ../AltNewtonCD && make)');
close all; clear all;
X = dlmread('Xfile'); % (10 samples x 15 features)
Y = dlmread('Yfile'); % (10 samples x 12 features)
lambdaLambda = 0.1;
lambdaTheta = 0.2;

% vanilla example:
[Lambda_vanilla, Theta_vanilla, stats_vanilla] = runAltNewtonCD(...
    Y, X, lambdaLambda, lambdaTheta);
figure('name', 'vanilla');
subplot(1,2,1); imagesc(Lambda_vanilla); colorbar;
title('Lambda'); set(gca,'dataAspectRatio',[1 1 1]);
subplot(1,2,2); imagesc(Theta_vanilla); colorbar;
title('Theta'); set(gca,'dataAspectRatio',[1 1 1]);

% options example:
options.tol = 1.0e-5;
[Lambda_vanilla, Theta_vanilla, stats_vanilla] = runAltNewtonCD(...
    Y, X, lambdaLambda, lambdaTheta, options);
figure('name', 'options');
subplot(1,2,1); imagesc(Lambda_vanilla);
title('Lambda'); set(gca,'dataAspectRatio',[1 1 1]);
subplot(1,2,2); imagesc(Theta_vanilla);
title('Theta'); set(gca,'dataAspectRatio',[1 1 1]);

% semi-supervised example:
partialX = X(1:end-4,:); % Y has 10 samples, X has 6 samples
[Lambda_semi, Theta_semi, stats_semi] = runAltNewtonCD(...
    Y, partialX, lambdaLambda, lambdaTheta);
figure('name', 'missing X');
subplot(1,2,1); imagesc(Lambda_semi);
title('Lambda'); set(gca,'dataAspectRatio',[1 1 1]);
subplot(1,2,2); imagesc(Theta_semi);
title('Theta'); set(gca,'dataAspectRatio',[1 1 1]);

% warm-starting example
options.Lambda0 = Lambda_vanilla;
options.Theta0 = Theta_vanilla;
options.refit = 0; % allowed to add edges
[Lambda_warm, Theta_warm, stats_warm] = runAltNewtonCD(...
    Y, X, 0.2*lambdaLambda, 0.5*lambdaTheta, options);
figure('name', 'warm-starting');
subplot(1,2,1); imagesc(Lambda_warm); 
title('Lambda'); set(gca,'dataAspectRatio',[1 1 1]);
subplot(1,2,2); imagesc(Theta_warm);
title('Theta'); set(gca,'dataAspectRatio',[1 1 1]);

% refitting example
options.Lambda0 = Lambda_vanilla;
options.Theta0 = Theta_vanilla;
options.refit = 1; % don't add any edges
[Lambda_refit, Theta_refit, stats_refit] = runAltNewtonCD(...
    Y, X, 0.2*lambdaLambda, lambdaTheta, options);
figure('name', 'refitting');
subplot(1,2,1); imagesc(Lambda_refit); colorbar;
title('Lambda'); set(gca,'dataAspectRatio',[1 1 1]);
subplot(1,2,2); imagesc(Theta_refit); colorbar;
title('Theta'); set(gca,'dataAspectRatio',[1 1 1]);

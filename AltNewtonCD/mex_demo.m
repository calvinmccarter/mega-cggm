mex -largeArrayDims cggmfast.cc cggmfast_mex.cc -DNDEBUG -I/usr/include/ -Ieigen3 -output cggmfast

% Mac OSX:
%setenv('DYLD_LIBRARY_PATH','');
%addpath('.'); % mexopts.sh
%mex -largeArrayDims cggmfast.cc cggmfast_mex.cc -DNDEBUG -Dchar16_t=uint16_t -I/usr/include/ -Ieigen3 -output cggmfast

n = 50;
p = 20;
q = 10;
X = 2*rand(n, p) - 1;
Lambda = eye(q, q);
Lambda(1,q) = 0.5;
Lambda(q,1) = 0.5;
Theta = full(sprand(p, q, 0.1));
Y = -X*Theta*inv(Lambda) + mvnrnd(zeros(q,1),inv(Lambda),n);
[estLambda, estTheta, stats] = cggmfast(Y, X, 0.05, 0.01);


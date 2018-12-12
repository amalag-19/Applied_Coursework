% original model perfect
clear all;
clc;

n = 1000;
X_1 = ones(n,1);
X_2 = rand(n,1);
X = [X_1, X_2];

% transformation on X
% X_3 = (X_2.^2).*log(X_2.^2);
% X_t = [X_1, X_2_sq, X_3];

beta0 = 1;
beta1 = 2;

% normal error
er = normrnd(0, 0.1, n, 1);

% original model
Y = beta0 + beta1*(X_2) + er;

% regression on original data
[b,bint,r] = regress(Y,X);
Y_cap = b(1) + b(2)*(X_2);

% residual plot to check 
figure, scatter(Y_cap, r);
title ('Plot of Residuals on original data vs. fitted values of the dependent variable');
xlabel('Y_cap');
ylabel('Residuals');

% qq plot to check deviation in distribution of Y from normality
m = 100000;
er_d = normrnd(0, 0.1, m, 1);
Y_d = beta0 + beta1*(X_2(50)) + er_d;
figure, qqplot(Y_d);

% Box-Cox transformation on Y
GM = geomean(Y);
V = @(lambda) (((Y.^lambda)-1)./(lambda*(GM^(lambda-1))));
V_cap = @(lambda) (X*(inv(X'*X))*(X')*V(lambda));
RSS = @(lambda) ((V(lambda)-V_cap(lambda))'*(V(lambda)-V_cap(lambda)));
lambda_max = fminbnd(RSS, -10, 10);
W = V(lambda_max);  % transformed Y
W_d = (((Y_d.^lambda_max)-1)./(lambda_max*(GM^(lambda_max-1))));

% regression on transformed data
[b_t,bint_t,r_t] = regress(W,X);
W_cap = b_t(1) + b_t(2)*(X_2);
figure, scatter(W_cap,r_t);
title ('Plot of Residuals on original data vs. X');
xlabel('W_cap');
ylabel('Residuals');

% qq plot to check distribution of Y after transformation
figure, qqplot(W_d);
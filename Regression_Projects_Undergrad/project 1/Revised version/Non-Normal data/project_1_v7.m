% original model non-normal (chi-square, t, F and beta distributions effect tested)
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

% error from chi-square distribution
df = 20; % degrees of freedom 
r = df*ones(n,1);
cr = chi2rnd(r);
ec = 0.1*cr;

% % error from t-distribution
% df = 4; % degrees of freedom 
% r = df*ones(n,1);
% tr = trnd(r);
% et = 0.1*tr;

% % error from F-distribution
% df1 = 1;    % degrees of freedom
% df2 = 1;    % degrees of freedom
% r1 = df1*ones(n,1);
% r2 = df2*ones(n,1);
% fr = frnd(r1,r2);
% ef = 0.1*fr;

% % error from Beta-distribution
% ap = 1000;    % degrees of freedom
% bp = 1000;    % degrees of freedom
% r1 = ap*ones(n,1);
% r2 = bp*ones(n,1);
% br = betarnd(r1,r2);
% eb = 0.1*br;

% original model
Y = beta0 + beta1*(X_2) + ec;

% regression on original data
[b,bint,r,rint,stats] = regress(Y,X);
Y_cap = b(1) + b(2)*(X_2);

% residual plot to check deviation from linearity & homoscedasticity and to detect the presence of outliers
rs = r/(sqrt(stats(4)));    % semi-studentized residuals
figure, scatter(X_2, rs);
title ('Plot of Semi-studentized Residuals of original data vs. X');
xlabel('X');
ylabel('Semi-studentized Residuals');

% qq plot to check deviation in distribution of Y from normality
m = 100000;
r_d = df*ones(m,1);
cr_d = trnd(r_d);
ec_d = 0.1*cr_d;
Y_d = beta0 + beta1*(X_2(n/2)) + ec_d;
figure, qqplot(Y_d);

% % qq plot to check deviation in distribution of Y from normality
% m = 100000;
% r_d = df*ones(m,1);
% tr_d = trnd(r_d);
% et_d = 0.1*tr_d;
% Y_d = beta0 + beta1*(X_2(50)) + et_d;
% figure, qqplot(Y_d);

% % qq plot to check deviation in distribution of Y from normality
% m = 100000;
% r1_d = df1*ones(m,1);
% r2_d = df2*ones(m,1);
% fr_d = frnd(r1_d,r2_d);
% ef_d = 0.1*fr_d;
% Y_d = beta0 + beta1*(X_2(50)) + ef_d;
% figure, qqplot(Y_d);

% % qq plot to check deviation in distribution of Y from normality
% m = 100000;
% r1_d = ap*ones(m,1);
% r2_d = bp*ones(m,1);
% br_d = frnd(r1_d,r2_d);
% eb_d = 0.1*br_d;
% Y_d = beta0 + beta1*(X_2(50)) + eb_d;
% figure, qqplot(Y_d);

% Box-Cox transformation on Y
GM = geomean(Y);
% V = @(lambda) ((log(Y))./(GM^(lambda-1)));
V = @(lambda) (((Y.^lambda)-1)./(lambda*(GM^(lambda-1))));
V_cap = @(lambda) (X*(inv(X'*X))*(X')*V(lambda));
RSS = @(lambda) ((V(lambda)-V_cap(lambda))'*(V(lambda)-V_cap(lambda)));
lambda_max = fminbnd(RSS, -10, 10);
W = V(lambda_max);  % transformed Y
W_d = (((Y_d.^lambda_max)-1)./(lambda_max*(GM^(lambda_max-1))));

% regression on transformed data
[b_t,bint_t,r_t,rint_t,stats_t] = regress(W,X);
W_cap = b_t(1) + b_t(2)*(X_2);
figure, scatter(W_cap,r_t);
title ('Plot of Residuals on original data vs. X');
xlabel('W_cap');
ylabel('Residuals');

% qq plot to check distribution of Y after transformation
figure, qqplot(W_d);
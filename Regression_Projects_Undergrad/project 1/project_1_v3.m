% original model non-linear and non-normal
clear;
clc;

n = 100;
X_1 = ones(n,1);
X_2 = 2+20*rand(n,1);
X = [X_1, X_2, X_2.^2];

% transformation on X
% X_3 = (X_2.^2).*log(X_2.^2);
% X_t = [X_1, X_2_sq, X_3];

beta0 = 0;
beta1 = 0;
beta2 = 9;

% error from t-distribution
df = 5; % degrees of freedom 
r = df*ones(n,1);
tr = trnd(r);
et = 0.1*tr;

% normal error
% var = zeros(n,1);
% for i = 1:n
%     var(i) = 0.0001*(i^2);
% end
% er = normrnd(0, var);

Y = beta0 + beta1*(X_2) + beta2*(X_2).^2 + et; % original model

% regression on original data
[b,bint,r] = regress(Y,X);
Y_cap = b(1) + b(2)*(X_2) + b(3)*((X_2).^2);

% residual plot to check 
figure, scatter(Y_cap, r);
title ('Plot of Residuals on original data vs. fitted values of the dependent variable');
xlabel('Y_cap');
ylabel('Residuals');

% qq plot to check deviation in distribution of Y from normality
m = 100000;
r_d = df*ones(m,1);
tr_d = trnd(r_d);
et_d = 0.1*tr_d;
Y_d = beta0 + beta1*(X_2(50)) + beta2*(X_2(50)).^2 + et_d;
figure, qqplot(Y_d);

% Box-Cox transformation on Y
GM = 1;
for i = 1:n
    GM = GM*Y(i);
end
GM = GM^(1/n);
V = @(lambda) (((Y.^lambda)-1)./(lambda*(GM^(lambda-1))));
V_cap = @(lambda) (X*(inv(X'*X))*(X')*V(lambda));
RSS = @(lambda) ((V(lambda)-V_cap(lambda))'*(V(lambda)-V_cap(lambda)));
lambda_max = fminbnd(RSS, -100, 100);
W = V(lambda_max);  % transformed Y
W_d = (((Y_d.^lambda_max)-1)./(lambda_max*(GM^(lambda_max-1))));

% regression on transformed data
[b_t,bint_t,r_t] = regress(W,X);
W_cap = b_t(1) + b_t(2)*(X_2) + b_t(3)*((X_2).^2);
figure, scatter(W_cap,r_t);
title ('Plot of Residuals on original data vs. X');
xlabel('W_cap');
ylabel('Residuals');

% qq plot to check distribution of Y after transformation
figure, qqplot(W_d);
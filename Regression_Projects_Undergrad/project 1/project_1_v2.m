clear;
clc;

n=100;
X_1 = ones(n,1);
X_2 = 2+20*rand(n,1);
X_2_sq = X_2.^2;
X = [X_1, X_2];
X_3 = X_2_sq.*log(X_2_sq);
X_t = [X_1, X_2_sq, X_3];

beta0 = 5;
beta1 = 3;
beta2 = 2;

r = 10*ones(n,1);
tr = trnd(r);
et = 0.1*tr;
var = zeros(n,1);
for i = 1:n
    var(i) = 0.0001*(i^2);
end
er = normrnd(0, 0.1, n, 1);

Y = beta0 + beta1*(X_2) + er; % + beta2*(X_2).^2 + et; % original model

[b,bint,r] = regress(Y,X);
Y_cap = b(1) + b(2)*(X_2); % + b(3)*((X_2).^2);
figure, scatter(Y_cap, r);
title ('Plot of Residuals on original data vs. X');
xlabel('X');
ylabel('Residuals');
figure, qqplot(Y);

GM=1;
for i=1:n
    GM=GM*Y(i);
end
GM=GM^(1/n);
V = @(lambda) (((Y.^lambda)-1)./(lambda*(GM^(lambda-1))));
V_cap = @(lambda) (X*(inv(X'*X))*(X')*V(lambda));
RSS = @(lambda) ((V(lambda)-V_cap(lambda))'*(V(lambda)-V_cap(lambda)));
lambda_max = fminbnd(RSS, -100, 100);
W = V(lambda_max);
W_cap = V_cap(lambda_max);

[b_t,bint_t,r_t] = regress(W,X);
W_cap = b_t(1) + b_t(2)*(X_2); % + b_t(3)*((X_2).^2);
figure, scatter(W_cap,r_t);
title ('Plot of Residuals on original data vs. X');
xlabel('X');
ylabel('Residuals');
figure, qqplot(W);
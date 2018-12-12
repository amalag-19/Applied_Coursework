% original model non-linear and heteroscedastic
clear all;
clc;

n = 1000;
alpha = 0.05;   % level of significance
X_1 = ones(n,1);
X_2 = 1+2*rand(n,1);
X = [X_1, X_2];

% transformation on X
% X_3 = (X_2.^2).*log(X_2.^2);
% X_t = [X_1, X_2_sq, X_3];

beta0 = 1;
beta1 = 2;

% normal error with variance not constant over the runs
% var = 0.1*abs(sin(0 + 2*pi*X_2));
var = 0.1*exp(X_2);
er = normrnd(0, var, n, 1);

% original model
Y = beta0 + beta1*(X_2.^2) + er;

% regression on original data
[b,bint,r,rint,stats] = regress(Y,X);
Y_cap = b(1) + b(2)*(X_2);

% residual plot to check deviation from linearity & homoscedasticity and to detect the presence of outliers
rs = r/(sqrt(stats(4)));    % semi-studentized residuals
figure, scatter(X_2, rs);
title ('Plot of Semi-studentized Residuals of original data vs. X');
xlabel('X');
ylabel('Semi-studentized Residuals');

% Lack of fit test on original model
m1 = 10;
y = n+1;
er_lf = zeros(m1, (y-1));
Y_lf = zeros(m1, (y-1));
Ymean_lf = zeros(1,(y-1));
YSSPE_lf = zeros(m1, (y-1));
YSSE_lf = zeros(m1, (y-1));
for j = 1:(y-1)
    for i = 1:m1
        er_lf(i,j) = normrnd(0,var(j));
        Y_lf(i,j) = beta0 + beta1*(X_2(j)) + er_lf(i,j);
    end
    Ymean_lf(j) = mean(Y_lf(:,j));
    for i = 1:m1
        YSSPE_lf(i,j) = (Y_lf(i,j)-Ymean_lf(j))^2;
        YSSE_lf(i,j) = (Y_lf(i,j)-(Y_cap(j)))^2;
    end
end
SSPE = sum(sum(YSSPE_lf));
SSE = sum(sum(YSSE_lf));
SSLF = SSE - SSPE;
MSLF = SSLF/((y-1)-2);
MSPE = SSPE/((m1*(y-1))-(y-1));
F_star = MSLF/MSPE;
F = finv((1-alpha),((y-1)-2),((m1*(y-1))-(y-1)));
if F_star > F
    disp('Regression Function is not linear for original model');
else disp('Regression Function is linear for original model');
end

% Outlier Analysis: Removing the outliers beyond 3 standard deviations 
k = 0;
for i = 1:n
    if (abs(rs(i))>3)
        k = [i,k];
    end
end

s = length(k);
km = zeros(1, s-1);
for i = 1:s-1
    km(i) = k(i);
end

X_2_temp1 = X_2;
var_temp1 = var;
er_temp1 = er;
Y_lf_temp1 = Y_lf;
for i=1:n
    for j=1:s-1
        if (km(j)==i)
            X_2_temp1(i) = 0;
            var_temp1(i) = 0;
            er_temp1(i) = 0;
            Y_lf_temp1(:,i) = zeros(m1,1);
        end
    end
end

X_2_temp2 = 0;
var_temp2 = 0;
er_temp2 = 0;
Y_lf_temp2 = zeros(m1,1);
for p=1:n
    if (X_2_temp1(p)~=0)
        X_2_temp2 = [X_2_temp1(p);X_2_temp2];
    end
    if (var_temp1(p)~=0)
        var_temp2 = [var_temp1(p);var_temp2];
    end
    if (er_temp1(p)~=0)
        er_temp2 = [er_temp1(p);er_temp2];
    end
    if (Y_lf_temp1(:,p)~=zeros(m1,1))
        Y_lf_temp2 = [Y_lf_temp1(:,p),Y_lf_temp2];
    end
end

y = length(X_2_temp2);  % = n-s+2
X_2_m = zeros((n-s+1),1);
var_m = zeros((n-s+1),1);
er_m = zeros((n-s+1),1);
Y_lf_m = zeros(m1,(n-s+1));
for i = 1:y-1
    X_2_m(i) = X_2_temp2(i);
    var_m(i) = var_temp2(i);
    er_m(i) = er_temp2(i);
    Y_lf_m(:,i) = Y_lf_temp2(:,i);
end

X_1_m = ones(y-1,1);
X_m = [X_1_m, X_2_m];

% original model modified by removing the outliers
Y_m = beta0 + beta1*(X_2_m.^2) + er_m;

% qq plot to check deviation in distribution of Y from normality
m2 = 100000;
mid = round((y-1)/2);
er_d = normrnd(0, var(mid), m2, 1);
Y_d = beta0 + beta1*(X_2(mid)) + er_d;
figure, qqplot(Y_d);

% Box-Cox transformation on modified Y
GM = geomean(Y_m);
V = @(lambda) (((Y_m.^lambda)-1)./(lambda*(GM^(lambda-1))));
V_cap = @(lambda) (X_m*(inv(X_m'*X_m))*(X_m')*V(lambda));
RSS = @(lambda) ((V(lambda)-V_cap(lambda))'*(V(lambda)-V_cap(lambda)));
lambda_max = fminbnd(RSS, -100, 100);
W = V(lambda_max);  % transformed Y
W_d = (((Y_d.^lambda_max)-1)./(lambda_max*(GM^(lambda_max-1))));
W_lf = (((Y_lf_m.^lambda_max)-1)./(lambda_max*(GM^(lambda_max-1))));

% regression on transformed data
[b_t,bint_t,r_t,rint_t,stats_t] = regress(W,X_m);
W_cap = b_t(1) + b_t(2)*(X_2_m);
rs_t = r_t/(sqrt(stats_t(4)));

% residual plot after transformation to check deviation from linearity & homoscedasticity and to detect the presence of outliers
figure, scatter(X_2_m, rs_t);
title ('Plot of Semi-studentized Residuals of transformed data vs. X');
xlabel('X');
ylabel('Semi-studentized Residuals');

% qq plot to check distribution of Y after transformation
figure, qqplot(W_d);

% Lack of fit test after transformation
Wmean_lf = zeros(1,(y-1));
WSSPE_lf = zeros(m1, (y-1));
WSSE_lf = zeros(m1, (y-1));
for j = 1:(y-1)
    Wmean_lf(j) = mean(W_lf(:,j));
    for i = 1:m1
        WSSPE_lf(i,j) = (W_lf(i,j)-Wmean_lf(j))^2;
        WSSE_lf(i,j) = (W_lf(i,j)-(W_cap(j)))^2;
    end
end
SSPE_t = sum(sum(WSSPE_lf));
SSE_t = sum(sum(WSSE_lf));
SSLF_t = SSE_t - SSPE_t;
MSLF_t = SSLF_t/((y-1)-2);
MSPE_t = SSPE_t/((m1*(y-1))-(y-1));
F_star_t = MSLF_t/MSPE_t;
F_t = finv((1-alpha),((y-1)-2),((m1*(y-1))-(y-1)));
if F_star_t > F_t
    disp('Regression Function is not linear after transformation');
else disp('Regression Function is linear after transformation');
end
clear;
clc;
n=100;
r = 10*ones(n,1);
X_1 = ones(n,1);
X_2 = 2+20*rand(n,1);
X_2_sq = X_2.^2;
X = [X_1,X_2_sq];
X_3 = X_2_sq.*log(X_2_sq);
X_t = [X_1, X_2, X_3];
% mu=zeros(n,1);
% sig=0.1*eye(n);
tr = trnd(r);
% m = normrnd(0, 0.1, n, 1);
m=0.1*tr; %normrnd(0, 0.1, n, 1);
Y = (X_2).^2 + m; % original model
qqplot(Y);

% sigma=zeros(n,1);
% Y=zeros(n,1);
% for i=1:n
%     sigma(i)=0.1*(i^3);
%     Y(i)=normrnd(X_2(i), sigma(i));
% end

%residual plots on the original data
k_old = (X_2-mean(X_2))./(sum((X_2-mean(X_2)).^2));
b1_old = sum(k_old.*Y);
b0_old = mean(Y)-b1_old*mean(X_2);
Y_cap = b0_old + b1_old*X_2;
Res_org = (Y-Y_cap);
figure, scatter(X_2,Res_org);
% axis([0 1800 0 2.5]);
title ('Plot of Residuals on original data vs. X');
xlabel('X');
ylabel('Residuals');
% hold on;

%for j=1:100
% estimation of lambda_max
GM=1;
for i=1:n
    GM=GM*Y(i);
end
GM=GM^(1/n);
V = @(lambda) (((Y.^lambda)-1)./(lambda*(GM^(lambda-1))));
V_cap = @(lambda) (X_t*(inv(X_t'*X_t))*(X_t')*V(lambda));
RSS=@(lambda) ((V(lambda)-V_cap(lambda))'*(V(lambda)-V_cap(lambda)));
lambda_max = fminbnd(RSS, -100, 100);
W = V(lambda_max);
W_cap = V_cap(lambda_max);
qqplot(V(lambda_max));
%Y=V(lambda_max);
%end

%residual plots on the transformed data
Res_new = V(lambda_max)-V_cap(lambda_max);
figure, scatter(X_2,Res_new);
title ('Plot of Residuals on transformed data vs. X');
xlabel('X');
ylabel('Residuals');

% hypothesis tests for b0=0 and b1=0 vs. the alternatives after transformation
b1 = sum(k_old.*V(lambda_max));
SSE = RSS(lambda_max);
MSE = SSE/(n-2);
s_b1 = (MSE/(sum((X_2-mean(X_2)).^2)))^(0.5);
t_b1 = b1/s_b1;
alpha = 0.05;
t_b1_alpha = tinv((1-(alpha/2)),(n-2));

b0 = mean(V(lambda_max))-b1*mean(X_2);
s_b0 = MSE*((1/n)+(((mean(X_2))^(2))/(sum((X_2-mean(X_2)).^2))));
t_b0 = abs(b0/s_b0);
t_b0_alpha = tinv((1-(alpha/2)),(n-2));

% Outlier Analysis
k=0;
for i=1:n
    if (abs(Res_new(i))>0.3)
        k=[i,k];
    end
end

s=length(k);
for i=1:s-1
    km(i) = k(i);
end

X_2_temp1=X_2;
tr_temp1=tr;
for i=1:n
    for j=1:s-1
        if (km(j)==i)
            X_2_temp1(i) = 0;
            tr_temp1(i) = 0;
        end
    end
end

X_2_temp2=[0];
tr_temp2=[0];
for p=1:n
    if (X_2_temp1(p)~=0)
        X_2_temp2 = [X_2_temp1(p);X_2_temp2];
    end
    if (tr_temp1(p)~=0)
        tr_temp2 = [tr_temp1(p);tr_temp2];
    end
end

y=length(X_2_temp2);
X_2_m=zeros((n-s+1),1);
tr_m=zeros((n-s+1),1);
for i=1:n-s+1
    X_2_m(i) = X_2_temp2(i);
    tr_m(i) = tr_temp2(i);
end

X_1_m = ones(y-1,1);
X_m = [X_1_m, X_2_m];
r=zeros(n-s+1,1);
sig=0.1*eye(n);
m_m=0.1*tr_m; %normrnd(0, 0.1, n, 1);
Y_m =(X_2_m.^2)+m_m;

% sigma=zeros(n,1);
% Y=zeros(n,1);
% for i=1:n
%     sigma(i)=0.1*(i^3);
%     Y(i)=normrnd(X_2(i), sigma(i));
% end

%residual plots on the original data after removing the outliers
k_old_m = (X_2_m-mean(X_2_m))./(sum((X_2_m-mean(X_2_m)).^2));
b1_old_m = sum(k_old_m.*Y_m);
b0_old_m = mean(Y_m)-b1_old_m*mean(X_2_m);
Y_cap_m = b0_old_m + b1_old_m*X_2_m;
Res_org_m = (Y_m-Y_cap_m);
figure, scatter(X_2_m,Res_org_m);
% axis([0 1800 0 2.5]);
title ('Plot of Residuals on original data vs. X removing the outliers');
xlabel('X_m');
ylabel('Residuals_m');
% hold on;

%for j=1:100
% estimation of lambda_max after removing the outliers
GM_m=1;
for i=1:n-s+1
    GM_m=GM_m*Y(i);
end
GM_m = GM_m^(1/(n-s+1));
V_m = @(lambda_m) (((Y_m.^lambda_m)-1)./(lambda_m*(GM_m^(lambda_m-1))));
V_cap_m = @(lambda_m) (X_m*(inv(X_m'*X_m))*(X_m')*V_m(lambda_m));
RSS_m = @(lambda_m) ((V_m(lambda_m)-V_cap_m(lambda_m))'*(V_m(lambda_m)-V_cap_m(lambda_m)));
lambda_max_m = fminbnd(RSS_m, -100, 100);
W_m = V_m(lambda_max_m);
W_cap_m = V_cap_m(lambda_max_m);
%Y=V(lambda_max);
%end

%residual plots on the transformed data after removing the outliers
Res_new_m = V_m(lambda_max_m)-V_cap_m(lambda_max_m);
figure, scatter(X_2_m,Res_new_m);
title ('Plot of Residuals on transformed data vs. X removing the outliers');
xlabel('X');
ylabel('Residuals');

% hypothesis tests for b0=0 and b1=0 vs. the alternatives after
% transformation removing the outliers
b1_m = sum(k_old_m.*V_m(lambda_max_m));
SSE_m = RSS_m(lambda_max_m);
MSE_m = SSE_m/(n-2);


s_b1_m = (MSE_m/(sum((X_2_m-mean(X_2_m)).^2)))^(0.5);
t_b1_m = b1_m/s_b1_m;
alpha = 0.05;
t_b1_alpha = tinv((1-(alpha/2)),(n-2));

b0_m = mean(V_m(lambda_max_m))-b1_m*mean(X_2_m);
s_b0_m = MSE_m*((1/(n-s+1))+(((mean(X_2_m))^(2))/(sum((X_2_m-mean(X_2_m)).^2))));
t_b0_m = abs(b0_m/s_b0_m);
t_b0_alpha = tinv((1-(alpha/2)),(n-2));
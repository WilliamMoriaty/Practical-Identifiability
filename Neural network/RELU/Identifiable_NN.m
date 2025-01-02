fig=figure(1);
clf();
subplot(2,2,1)
alpha = load("Nalpha_10.txt");
beta = load("Nbeta_10.txt");
omega = load("Nomega_10.txt");
N0 = 500;
t = linspace(0,1,N0);
y = sin(2*pi*t);
N1 = 301;
t1 = linspace(0,1,N1);
phi = zeros(N1,1);
for i = 1:N1
phi(i) = NN(t1(i),alpha,beta,omega);
end

%plot(t,y,'ro')
hold on
plot(t1,phi,'k-')

% S1: 500 data point
S1 = Sens_NN(alpha,beta,omega,N0);
F1 = S1'*S1;
[U,Sigma1,~]=svd(F1);
r = min(find(diag(Sigma1)<1e-6));
%Ur = U;
Ur = U(:,r:end);
m = size(alpha,1);
sigma_th = 1e-3; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
S_alpha = zeros(1,m);
S_beta = zeros(1,m);
S_omega = zeros(1,m);
    for j= 1:m
        S_alpha (j) = omega(j) * Dsigma(alpha(j)*((i-1)/(N1-1))+beta(j)) * (i-1)/(N1-1);
        S_beta (j) = omega(j) * Dsigma(alpha(j)*((i-1)/(N1-1))+beta(j));
        S_omega (j) = sigma(alpha(j)*((i-1)/(N1-1))+beta(j));        
    end
    S = [S_alpha S_beta S_omega];
    Var_NN(i)= sigma_th * dot(S * Ur, S*Ur);
end

x = t1;
Y_fit = phi;
CI=1.96*sqrt(Var_NN');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 

%%
subplot(2,2,2)
alpha = load("Nalpha_20.txt");
beta = load("Nbeta_20.txt");
omega = load("Nomega_20.txt");
N0 = 500;
t = linspace(0,1,N0);
y = sin(2*pi*t);
N1 = 301;
t1 = linspace(0,1,N1);
phi = zeros(N1,1);
for i = 1:N1
phi(i) = NN(t1(i),alpha,beta,omega);
end

%plot(t,y,'ko')


% S1: 500 data point
S1 = Sens_NN(alpha,beta,omega,N0);
F1 = S1'*S1;
[U,Sigma1,~]=svd(F1);
r = min(find(diag(Sigma1)<1e-6));
Ur = U;
%Ur = U(:,r:end);
m = size(alpha,1);
sigma_th = 1e-3; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
S_alpha = zeros(1,m);
S_beta = zeros(1,m);
S_omega = zeros(1,m);
    for j= 1:m
        S_alpha (j) = omega(j) * Dsigma(alpha(j)*((i-1)/(N1-1))+beta(j)) * (i-1)/(N1-1);
        S_beta (j) = omega(j) * Dsigma(alpha(j)*((i-1)/(N1-1))+beta(j));
        S_omega (j) = sigma(alpha(j)*((i-1)/(N1-1))+beta(j));        
    end
    S = [S_alpha S_beta S_omega];
    Var_NN(i)= sigma_th * dot(S * Ur, S*Ur);
end

x = t1;
Y_fit = phi;
CI=1.96*sqrt(Var_NN');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(t1,phi,'r-')
% 
%%
subplot(2,2,3)
alpha = load("Nalpha_40.txt");
beta = load("Nbeta_40.txt");
omega = load("Nomega_40.txt");
N0 = 20;
t = linspace(0,1,N0);
y = sin(2*pi*t);
N1 = 1001;
t1 = linspace(0,1,N1);
phi = zeros(N1,1);
for i = 1:N1
phi(i) = NN(t1(i),alpha,beta,omega);
end


% S1: 500 data point
S1 = Sens_NN(alpha,beta,omega,N0);
F1 = S1'*S1;
[U,Sigma1,~]=svd(F1);
r = min(find(diag(Sigma1)<1e-6));
Ur = U;
%Ur = U(:,r:end);
m = size(alpha,1);
sigma_th = 1e-3; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
S_alpha = zeros(1,m);
S_beta = zeros(1,m);
S_omega = zeros(1,m);
    for j= 1:m
        S_alpha (j) = omega(j) * Dsigma(alpha(j)*((i-1)/(N1-1))+beta(j)) * (i-1)/(N1-1);
        S_beta (j) = omega(j) * Dsigma(alpha(j)*((i-1)/(N1-1))+beta(j));
        S_omega (j) = sigma(alpha(j)*((i-1)/(N1-1))+beta(j));        
    end
    S = [S_alpha S_beta S_omega];
    Var_NN(i)= sigma_th * dot(S * Ur, S*Ur);
end

x = t1;
Y_fit = phi;
CI=1.96*sqrt(Var_NN');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(t,y,'g-o')
hold on
plot(t1,phi,'r-')



alpha = load("alpha_40.txt");
beta = load("beta_40.txt");
omega = load("omega_40.txt");

phi = zeros(N1,1);
for i = 1:N1
phi(i) = NN(t1(i),alpha,beta,omega);
end
hold on
plot(t1,phi,'k-')
% 
%%
subplot(2,2,4)
alpha = load("Nalpha_60.txt");
beta = load("Nbeta_60.txt");
omega = load("Nomega_60.txt");
N0 = 20;
t = linspace(0,1,N0);
y = sin(2*pi*t);
N1 = 300;
t1 = linspace(0,1,N1);
phi = zeros(N1,1);
for i = 1:N1
phi(i) = NN(t1(i),alpha,beta,omega);
end

% S1: 500 data point
S1 = Sens_NN(alpha,beta,omega,N1);
F1 = S1'*S1;
[U,Sigma1,~]=svd(F1);
r = min(find(diag(Sigma1)<1e-6));
Ur = U;
%Ur = U(:,r:end);
m = size(alpha,1);
sigma_th = 1e-3; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
S_alpha = zeros(1,m);
S_beta = zeros(1,m);
S_omega = zeros(1,m);
    for j= 1:m
        S_alpha (j) = omega(j) * Dsigma(alpha(j)*((i-1)/(N1-1))+beta(j)) * (i-1)/(N1-1);
        S_beta (j) = omega(j) * Dsigma(alpha(j)*((i-1)/(N1-1))+beta(j));
        S_omega (j) = sigma(alpha(j)*((i-1)/(N1-1))+beta(j));        
    end
    S = [S_alpha S_beta S_omega];
    Var_NN(i)= sigma_th * dot(S * Ur, S*Ur);
end



x = t1;
Y_fit = phi;
CI=1.96*sqrt(Var_NN');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(t,y,'g-o')
hold on
plot(t1,phi,'r-')

alpha = load("alpha_60.txt");
beta = load("beta_60.txt");
omega = load("omega_60.txt");

phi = zeros(N1,1);
for i = 1:N1
phi(i) = NN(t1(i),alpha,beta,omega);
end

plot(t1,phi,'k-')

% 
%%
function Ds = Dsigma(x)
% Derivative of Relu(x)
if x > 0
    Ds = 1;
else
    Ds = 0;
end
end

function s = sigma(x)
if x > 0
    s = x;
else
    s = 0;
end
end

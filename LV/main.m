% File: lotka_volterra.m

% Parameters
para = zeros(4,1);
para(1) = 0.5;   % Prey growth rate: alpha = para(1)
para(2) = 0.1;    % Predation rate: beta = para(2)
para(3) = 0.06; % Reproduction rate of predator: delta = para(3)
para(4) = 2.0;   % Predator death rate: gamma = para(4)

% Initial conditions for prey (x0) and predator (y0)
x0 = 40;  % Initial prey population
y0 = 9;   % Initial predator population
N = 21;
% Time span
tspan = linspace(0,20,N);

% Solve the ODEs
[t, Y] = ode45(@(t,y)lotka_volterra_eq(t,y,para), tspan, [x0 y0],para);

% Extract the prey and predator populations
prey = Y(:, 1);
predator = Y(:, 2);
yexp1 = prey;
yexp2 = predator;
%%
% options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',100000,'MaxIterations',20000,...
%     'algorithm','levenberg-marquardt','StepTolerance',1e-8);
% 
% y0 = [x0;y0];
% lb = [0 0 0 0];
% ub = [2 2 0.2 4];
% para_new0 = para;
% 
% [para_new,resnorm,residual,exitflag,output,lambda,jacobian] = ...
%     lsqnonlin(@Objective_lv,para_new0,lb,ub,options,tspan,y0,1,yexp1,yexp2);       
 para_new = para;
 y0 = [x0;y0];
%% Fisher information matrix
tspan1=tspan;
x0 = [0;0;y0];
[~, Xi4] = ode45(@(t,x)partial_theta4(t, x,para_new), tspan,x0);
[~, Xi3] = ode45(@(t,x)partial_theta3(t, x,para_new), tspan,x0);
[~, Xi2] = ode45(@(t,x)partial_theta2(t, x,para_new), tspan,x0);
[~, Xi1] = ode45(@(t,x)partial_theta1(t, x,para_new), tspan,x0);

Patial_phi1_theta = [Xi1(:,1) Xi2(:,1) Xi3(:,1) Xi4(:,1)];
Patial_phi2_theta = [Xi1(:,2) Xi2(:,2) Xi3(:,2) Xi4(:,2)];
S = [Patial_phi1_theta;Patial_phi2_theta]; %Sensitive Matrix
%S = Patial_phi1_theta;
%S = Patial_phi2_theta;
F = S'* S;
[U,Sigma,~]=svd(F);
%% RAU
nn = size(S,1);
s1=S(:,1);
A = S(:,2:4);
ss1 = (eye(nn)-A*pinv(A))*s1;

s2=S(:,2);
A = [S(:,1) S(:,3:4)];
ss2 = (eye(nn)-A*pinv(A))*s2;

s3=S(:,3);
A = [S(:,1:2) S(:,4)];
ss3 = (eye(nn)-A*pinv(A))*s3;

s4=S(:,4);
A = S(:,1:3);
ss4 = (eye(nn)-A*pinv(A))*s4;

SSS=[norm(ss1,"inf") norm(ss2,"inf") norm(ss3,"inf") norm(ss4,"inf")];

% %% UQ
fig1=figure(1);
clf()
subplot(1,2,1)
N1=201;
tspan = linspace(0,20,N1);
tspan1=tspan;
x0 = [0;0;y0];
[~, Xi4] = ode45(@(t,x)partial_theta4(t, x,para_new),tspan,x0);
[~, Xi3] = ode45(@(t,x)partial_theta3(t, x,para_new),tspan,x0);
[~, Xi2] = ode45(@(t,x)partial_theta2(t, x,para_new),tspan,x0);
[~, Xi1] = ode45(@(t,x)partial_theta1(t, x,para_new),tspan,x0);

Patial_phi1_theta = [Xi1(:,1) Xi2(:,1) Xi3(:,1) Xi4(:,1)];
Patial_phi2_theta = [Xi1(:,2) Xi2(:,2) Xi3(:,2) Xi4(:,2)];
%S = [Patial_phi1_theta;Patial_phi2_theta]; %Sensitive Matrix
[t, phi] = ode45(@(t,y)lotka_volterra_eq(t,y,para), tspan, y0,para);
% F = S'* S;
% [U,Sigma,~]=svd(F);
 r = min(find(diag(Sigma)<1e2));
 Ur = U(:,r:end);
%Ur = U;
S = Patial_phi1_theta;
sigma_th = 5e-1; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (S(i,:) * Ur)*(S(i,:) * Ur)';
end

x = tspan;
Y_fit = phi(:,1);
CI=1.96*sqrt(Var_NN');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

hold on

S = Patial_phi2_theta;
%Ur = U;
%sigma_th = 2e5; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (S(i,:) * Ur)*(S(i,:) * Ur)';
end
x = tspan;
Y_fit = phi(:,2);
CI=1.96*sqrt(Var_NN');
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% 
% Plot results

plot(t, phi(:,1), 'r', 'LineWidth', 2); 
hold on;
plot(t, phi(:,2), 'b', 'LineWidth', 2);
hold on;
ylim([0,60])
N3=N;
tspan2=linspace(0,20,N3);
[t2, phi1] = ode45(@(t,y)lotka_volterra_eq(t,y,para), tspan2, y0,para);

plot(t2, phi1(:,1), 'ro', 'LineWidth', 2); 
hold on;
plot(t2, phi1(:,2), 'b*', 'LineWidth', 2);

xlabel('Time');
ylabel('Population');
lgd = legend("95% CI prey","95% CI predator","Prey","Predator","Prey Data","Predator Data");
lgd.Location = 'best';
lgd.ItemTokenSize = [10,6];
lgd.FontWeight = 'bold';
lgd.Box='off';
lgd.ItemTokenSize = [10,6];
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

subplot(1,2,2)
Ur = U;
S = Patial_phi1_theta;
sigma_th = 1e-5; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (S(i,:) * Ur)*(S(i,:) * Ur)';
end

x = tspan;
Y_fit = phi(:,1);
CI=1.96*sqrt(Var_NN');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

hold on

S = Patial_phi2_theta;
Ur = U;
%sigma_th = 1e-4; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (S(i,:) * Ur)*(S(i,:) * Ur)';
end
x = tspan;
Y_fit = phi(:,2);
CI=1.96*sqrt(Var_NN');
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% 
% Plot results

plot(t, phi(:,1), 'r', 'LineWidth', 2); hold on;
plot(t, phi(:,2), 'b', 'LineWidth', 2);
N3=N;
tspan2=linspace(0,20,N3);
[t2, phi1] = ode45(@(t,y)lotka_volterra_eq(t,y,para), tspan2, y0,para);

plot(t2, phi1(:,1), 'ro', 'LineWidth', 2); 
hold on;
plot(t2, phi1(:,2), 'b*', 'LineWidth', 2);
xlabel('Time');
ylabel('Population');

fig2=figure(2);
clf();
subplot(2,2,3)
x=1:4;
bar(x,SSS)

ylabel('$$\|(I_n - AA^{\dagger}) \mathbf{s}_i\|_{\infty}$$', 'Interpreter', 'latex');
set(gca,'xticklabel',{'\alpha','\beta','\delta','\gamma'})

subplot(2,2,1)
x=1:4;
bar(x,[(Sigma(1,1)),(Sigma(2,2)),(Sigma(3,3)),(Sigma(4,4))])
hold on
plot([0,6],[100,100],'LineStyle','--','Color','k','LineWidth',1.5)
ylim([1,1e8])
xlim([0,5])
ylabel('Singular Value of FIM');
set(gca,'xticklabel',{'U_1^T\theta','U_2^T\theta','U_3^T\theta','U_4^T\theta'},'YScale','log')

subplot(2,2,2)
alphaData = ones(4,4);  % 初始化为全不透明
alphaData(:, r:4) = 0.2;  % 设置右半部分透明度为 0.2

imagesc(1:4,1:4,abs(U),'AlphaData',alphaData);

% 设置 colormap
cmap = othercolor('BuDRd_12');
colormap(cmap);  % 可以选择其他 colormap 例如 'jet', 'hot', 'cool' 等
clim([0,1.1])
% 添加 colorbar 并设置标签
c = colorbar;
c.Label.String = '|\partial U_i^T\theta/\partial \theta_j|';  % 设置 colorbar 的标签
c.FontSize = 12;  % 调整字体大小
c.Label.FontWeight = 'bold';  % 设置字体加粗
set(gca,'xticklabel',{'U_1','U_2','U_3','U_4'},...
    'YTick',1:4,'yticklabel',{'\alpha','\beta','\delta','\gamma'})
title('\theta=[\alpha,\beta,\delta,\gamma]','FontSize',14,'FontWeight','bold')
subplot(2,2,4)
NN=[20,50,100,200,400,800,1000,2000,4000,8000,1e4];
minsvd=zeros(1,size(NN,2));
maxsvd=zeros(1,size(NN,2));
%% Fisher information matrix
for i=1:size(NN,2)
tspan = linspace(0,20,NN(i));
tspan1=tspan;
x0 = [0;0;y0];
[~, Xi4] = ode45(@(t,x)partial_theta4(t, x,para_new), tspan,x0);
[~, Xi3] = ode45(@(t,x)partial_theta3(t, x,para_new), tspan,x0);
[~, Xi2] = ode45(@(t,x)partial_theta2(t, x,para_new), tspan,x0);
[~, Xi1] = ode45(@(t,x)partial_theta1(t, x,para_new), tspan,x0);

Patial_phi1_theta = [Xi1(:,1) Xi2(:,1) Xi3(:,1) Xi4(:,1)];
Patial_phi2_theta = [Xi1(:,2) Xi2(:,2) Xi3(:,2) Xi4(:,2)];
S = [Patial_phi1_theta;Patial_phi2_theta]; %Sensitive Matrix
Patial_phi1_theta = [Xi1(:,1) Xi2(:,1) Xi3(:,1) Xi4(:,1)];
Patial_phi2_theta = [Xi1(:,2) Xi2(:,2) Xi3(:,2) Xi4(:,2)];
S = [Patial_phi1_theta;Patial_phi2_theta]; %Sensitive Matrix
%S = Patial_phi1_theta;
%S = Patial_phi2_theta;
F = S'* S;
minsvd(i)=min(svd(F));
maxsvd(i)=max(svd(F));
end

% plot(NN,minsvd,'o-')
% hold on
plot(NN,minsvd./maxsvd,'*-')
xlabel('Number of Sample')
ylabel('\sigma_{min}/\sigma_{max}')
set(gca,'XScale','log')

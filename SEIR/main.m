% File: SEIR.m

% Parameters
para = zeros(3,1);
para(1) = 0.2;   % beta = para(1)
para(2) = 0.1;    % sigma = para(2)
para(3) = 0.06;    % gamma = para(3)

% Initial conditions for prey (x0) and predator (y0)
y10 = 0.99;  % Initial y1
y20 = 0.01;  % Initial y2
y30 = 0.0;  % Initial y3
y40 = 0.0;  % Initial y4
yzero = [y10 y20 y30 y40];
N = 21;
% Time span
tspan = linspace(0,150,N);

% Solve the ODEs
[t, Y] = ode45(@(t,y)SEIR_eq(t,y,para), tspan, yzero,para);

% % Extract the prey and predator populations
% y2 = Y(:, 2);
% y3 = Y(:, 3);
% yexp1 = y2;
% yexp2 = y3;
%%
% options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',100000,'MaxIterations',20000,...
%     'algorithm','levenberg-marquardt','StepTolerance',1e-8);
% 
 y0 = yzero';
% lb = [0 0 0];
% ub = [10.0 10.0 10.0];
para_new = para;

% [para_new,resnorm,residual,exitflag,output,lambda,jacobian] = ...
%     lsqnonlin(@Objective_seir,para_new0,lb,ub,options,tspan,y0,1,yexp1,yexp2);       
%% Fisher information matrix
tspan1=tspan;
x0 = zeros(4,1);
x0 = [x0;y0];
[~, Xi3] = ode45(@(t,x)partial_theta3(t, x,para_new), tspan1,x0);
[~, Xi2] = ode45(@(t,x)partial_theta2(t, x,para_new), tspan1,x0);
[~, Xi1] = ode45(@(t,x)partial_theta1(t, x,para_new), tspan1,x0);

Patial_phi1_theta = [Xi1(:,1) Xi2(:,1) Xi3(:,1)];
Patial_phi2_theta = [Xi1(:,2) Xi2(:,2) Xi3(:,2)];
Patial_phi3_theta = [Xi1(:,3) Xi2(:,3) Xi3(:,3)];
Patial_phi4_theta = [Xi1(:,4) Xi2(:,4) Xi3(:,4)];
%S = [Patial_phi1_theta;Patial_phi2_theta;Patial_phi3_theta;Patial_phi4_theta]; %Sensitive Matrix
%S = Patial_phi1_theta;
%S = Patial_phi2_theta;
%S = [Patial_phi2_theta;Patial_phi3_theta]; %Sensitive Matrix
S = Patial_phi3_theta;
%S = Patial_phi4_theta
F = S'* S;
[U,Sigma,~]=svd(F);
S123 = [Patial_phi1_theta;Patial_phi2_theta;Patial_phi3_theta]; %Sensitive Matrix
F123 = S123'* S123;
[U123,Sigma123,~]=svd(F123);
S13 = [Patial_phi1_theta;Patial_phi3_theta]; %Sensitive Matrix
F13 = S13'* S13;
[U13,Sigma13,~]=svd(F13);
S23 = [Patial_phi2_theta;Patial_phi3_theta]; %Sensitive Matrix
F23 = S23'* S23;
[U23,Sigma23,~]=svd(F23);
%% RAU
nn = size(S,1);
s1=S(:,1);
A = [S(:,3) S(:,2)];
ss1 = (eye(nn)-A*pinv(A))*s1;

s2=S(:,2);
A = [S(:,1) S(:,3)];
ss2 = (eye(nn)-A*pinv(A))*s2;

s3=S(:,3);
A = [S(:,1) S(:,2)];
ss3 = (eye(nn)-A*pinv(A))*s3;
SSS=[norm(ss1,"inf") norm(ss2,"inf") norm(ss3,"inf")];

nn = size(S123,1);
s1=S123(:,1);
A = [S123(:,3) S123(:,2)];
ss1 = (eye(nn)-A*pinv(A))*s1;

s2=S123(:,2);
A = [S123(:,1) S123(:,3)];
ss2 = (eye(nn)-A*pinv(A))*s2;

s3=S123(:,3);
A = [S123(:,1) S123(:,2)];
ss3 = (eye(nn)-A*pinv(A))*s3;
SSS123=[norm(ss1,"inf") norm(ss2,"inf") norm(ss3,"inf")];

nn = size(S13,1);
s1=S13(:,1);
A = [S13(:,3) S13(:,2)];
ss1 = (eye(nn)-A*pinv(A))*s1;

s2=S13(:,2);
A = [S13(:,1) S13(:,3)];
ss2 = (eye(nn)-A*pinv(A))*s2;

s3=S13(:,3);
A = [S13(:,1) S13(:,2)];
ss3 = (eye(nn)-A*pinv(A))*s3;
SSS13=[norm(ss1,"inf") norm(ss2,"inf") norm(ss3,"inf")];

nn = size(S23,1);
s1=S23(:,1);
A = [S23(:,3) S23(:,2)];
ss1 = (eye(nn)-A*pinv(A))*s1;

s2=S23(:,2);
A = [S23(:,1) S23(:,3)];
ss2 = (eye(nn)-A*pinv(A))*s2;

s3=S23(:,3);
A = [S23(:,1) S23(:,2)];
ss3 = (eye(nn)-A*pinv(A))*s3;
SSS23=[norm(ss1,"inf") norm(ss2,"inf") norm(ss3,"inf")];

%% UQ

fig1=figure(1);
clf()
N1=201;
tspan = linspace(0,150,N1);
tspan1=tspan;
x0 = zeros(4,1);
x0 = [x0;y0];
[~, Xi3] = ode45(@(t,x)partial_theta3(t, x,para_new), tspan,x0);
[~, Xi2] = ode45(@(t,x)partial_theta2(t, x,para_new), tspan,x0);
[~, Xi1] = ode45(@(t,x)partial_theta1(t, x,para_new), tspan,x0);

Patial_phi1_theta = [Xi1(:,1) Xi2(:,1) Xi3(:,1)];
Patial_phi2_theta = [Xi1(:,2) Xi2(:,2) Xi3(:,2)];
Patial_phi3_theta = [Xi1(:,3) Xi2(:,3) Xi3(:,3)];
Patial_phi4_theta = [Xi1(:,4) Xi2(:,4) Xi3(:,4)];
%S = [Patial_phi1_theta;Patial_phi2_theta;Patial_phi3_theta;Patial_phi4_theta]; %Sensitive Matrix
%S = [Patial_phi1_theta;Patial_phi2_theta]; %Sensitive Matrix
[t, phi] = ode45(@(t,y)SEIR_eq(t,y,para), tspan, y0,para);

Ur = U;
S = Patial_phi3_theta;
sigma_th = 1e-5; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (S(i,:) * Ur)*(S(i,:) * Ur)';
end

x = tspan;
Y_fit = phi(:,3);
CI=1.96*sqrt(Var_NN');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on

% S = Patial_phi4_theta;
% Ur = U;
% %sigma_th = 1e4; % parameter estimation
% Var_NN = zeros(1,N1);
% for i = 1:N1
% Var_NN(i)= sigma_th * (S(i,:) * Ur)*(S(i,:) * Ur)';
% end
% x = tspan;
% Y_fit = phi(:,4);
% CI=1.96*sqrt(Var_NN');
% fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% 
% Plot results

 plot(t, phi(:,3), 'r', 'LineWidth', 2); 
% hold on
% plot(t, phi(:,4), 'b', 'LineWidth', 2);

 plot(t(1:10:end), phi(1:10:end,3), 'ro', 'LineWidth', 2); 
% hold on;
% plot(t(1:10:end), phi(1:10:end,4), 'b*', 'LineWidth', 2);

xlabel('Time');
ylabel('Ratio');
lgd = legend("95% CI","I","data I");
lgd.ItemTokenSize = [10,6];
% lgd.Position = [0.746 0.301 0.139 0.326];
lgd.Box = 'off';
lgd.FontWeight = 'bold';
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

%%
fig2=figure(2);
clf();
set(gcf,'Position',[-32,463,1620,299])

subplot(1,4,1)
N1=201;
tspan = linspace(0,150,N1);
tspan1=tspan;
x0 = zeros(4,1);
x0 = [x0;y0];
[~, Xi3] = ode45(@(t,x)partial_theta3(t, x,para_new), tspan,x0);
[~, Xi2] = ode45(@(t,x)partial_theta2(t, x,para_new), tspan,x0);
[~, Xi1] = ode45(@(t,x)partial_theta1(t, x,para_new), tspan,x0);

Patial_phi1_theta = [Xi1(:,1) Xi2(:,1) Xi3(:,1)];
Patial_phi2_theta = [Xi1(:,2) Xi2(:,2) Xi3(:,2)];
Patial_phi3_theta = [Xi1(:,3) Xi2(:,3) Xi3(:,3)];
Patial_phi4_theta = [Xi1(:,4) Xi2(:,4) Xi3(:,4)];
%S = [Patial_phi1_theta;Patial_phi2_theta;Patial_phi3_theta;Patial_phi4_theta]; %Sensitive Matrix
%S = [Patial_phi1_theta;Patial_phi2_theta]; %Sensitive Matrix
[t, phi] = ode45(@(t,y)SEIR_eq(t,y,para), tspan, y0,para);
% F = S'* S;
% [U,Sigma,~]=svd(F);
r = min(find(diag(Sigma)<1e-1));
Ur = U(:,r:end);
%Ur = U;
S = Patial_phi3_theta;
sigma_th = 1e-2; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (S(i,:) * Ur)*(S(i,:) * Ur)';
end

x = tspan;
Y_fit = phi(:,3);
CI=1.96*sqrt(Var_NN');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

hold on

 plot(t, phi(:,3), 'r', 'LineWidth', 2); 
 hold on
%  plot(t, phi(:,4), 'b', 'LineWidth', 2);

 plot(t(1:10:end), phi(1:10:end,3), 'ro', 'LineWidth', 2); 
 hold on;
% plot(t(1:10:end), phi(1:10:end,4), 'b*', 'LineWidth', 2);
ylim([0 0.25])
xlabel('Time');
ylabel('Ratio');
% lg = legend("95% CI I","95% CI R","I","R","data I","data R");
lgd = legend("95% CI","I","data I");
lgd.ItemTokenSize = [10,6];
% lgd.Position = [0.746 0.301 0.139 0.326];
lgd.Box = 'off';
lgd.FontWeight = 'bold';
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

subplot(1,4,3)
x=1:3;
bar(x,[SSS;SSS13;SSS23;SSS123])
ylim([0.01,10])
ylabel('$$\|(I - A_iA^{\dagger}_i) \mathbf{s}_i\|_{\infty}$$', 'Interpreter', 'latex');
set(gca,'xticklabel',{'\beta','\sigma','\gamma'},'YScale','log')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off

subplot(1,4,2)
x=1:3;
bar(x,[Sigma(3,3),Sigma(2,2),Sigma(1,1);...
    Sigma13(3,3),Sigma13(2,2),Sigma13(1,1);...
    Sigma23(3,3),Sigma23(2,2),Sigma23(1,1);...
    Sigma123(3,3),Sigma123(2,2),Sigma123(1,1)])
hold on
plot([0,4],[0.1,0.1],'LineStyle','--','LineWidth',1.5,'Color','k')
xlim([0.5 3.5])
ylabel('Eigenvalue of FIM');
set(gca,'xticklabel',{'U_1^T\theta','U_2^T\theta','U_3^T\theta'},'YScale','log','FontSize',18,'FontWeight','bold')
lgd = legend("data I","data I & S","data I & E","data I & S & E");
lgd.Box = 'off';
lgd.Location = 'best';
lgd.FontWeight = 'bold';
lgd.ItemTokenSize = [10,6];
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

subplot(1,4,4)
NN=[10,15,20,40,60,80,100]+1;
minsvd3=zeros(1,size(NN,2));
minsvd13=zeros(1,size(NN,2));
minsvd23=zeros(1,size(NN,2));
minsvd123=zeros(1,size(NN,2));

maxsvd3=zeros(1,size(NN,2));
maxsvd13=zeros(1,size(NN,2));
maxsvd23=zeros(1,size(NN,2));
maxsvd123=zeros(1,size(NN,2));
%% Fisher information matrix
for i=1:size(NN,2)
tspan = linspace(0,150,NN(i));
tspan1=tspan;
x0 = zeros(4,1);
x0 = [x0;y0];
[~, Xi3] = ode15s(@(t,x)partial_theta3(t, x,para_new), tspan,x0);
[~, Xi2] = ode15s(@(t,x)partial_theta2(t, x,para_new), tspan,x0);
[~, Xi1] = ode15s(@(t,x)partial_theta1(t, x,para_new), tspan,x0);

Patial_phi1_theta = [Xi1(:,1) Xi2(:,1) Xi3(:,1)];
Patial_phi2_theta = [Xi1(:,2) Xi2(:,2) Xi3(:,2)];
Patial_phi3_theta = [Xi1(:,3) Xi2(:,3) Xi3(:,3)];
Patial_phi4_theta = [Xi1(:,4) Xi2(:,4) Xi3(:,4)];
S = Patial_phi3_theta;
F = S'* S;
minsvd3(i)=min(svd(F));
maxsvd3(i)=max(svd(F));
S = [Patial_phi1_theta;Patial_phi3_theta]; %Sensitive Matrix
F = S'* S;
minsvd13(i)=min(svd(F));
maxsvd13(i)=max(svd(F));
S = [Patial_phi2_theta;Patial_phi3_theta]; %Sensitive Matrix
F = S'* S;
minsvd23(i)=min(svd(F));
maxsvd23(i)=max(svd(F));
S = [Patial_phi1_theta;Patial_phi3_theta;Patial_phi3_theta]; %Sensitive Matrix
F = S'* S;
minsvd123(i)=min(svd(F));
maxsvd123(i)=max(svd(F));
end

plot(NN,minsvd3./maxsvd3,'o-','LineWidth',1.2,'MarkerSize',8)
hold on
plot(NN,minsvd13./maxsvd13,'*-','LineWidth',1.2,'MarkerSize',8)
hold on
plot(NN,minsvd23./maxsvd23,'LineStyle','-','Marker','square','LineWidth',1.2,'MarkerSize',8)
hold on
plot(NN,minsvd123./maxsvd123,'LineStyle','-','Marker','pentagram','LineWidth',1.2,'MarkerSize',8)
xlabel('Number of time points')
ylabel('\xi','FontSize',18,'FontWeight','bold')
ylim([0.0001,0.03])
xlim([10,101])
set(gca,'YScale','log')
lgd = legend("data I","data I & S","data I & E","data I & S & E");
lgd.Box = 'off';
lgd.Position = [0.757 0.622 0.068 0.239];
lgd.FontWeight = 'bold';
lgd.ItemTokenSize = [10,6];
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

% 
fig3=figure(3);
clf();
set(gcf,'Position',[-67,463,1620,249])
subplot(1,4,1)
alphadata = ones(3,3);  % 初始化为全不透明
alphadata(:, r:3) = 0.2;  % 设置右半部分透明度为 0.2

imagesc(1:3,1:3,abs(U),'Alphadata',alphadata);

% 设置 colormap
cmap = othercolor('BuDRd_12');
colormap(cmap);  % 可以选择其他 colormap 例如 'jet', 'hot', 'cool' 等
clim([0,1.0])
% 添加 colorbar 并设置标签
c = colorbar;
c.Label.String = '|\partial U_i^T\theta/\partial \theta_j|';  % 设置 colorbar 的标签
c.FontSize = 15;  % 调整字体大小
c.Label.FontWeight = 'bold';  % 设置字体加粗
set(gca,'xticklabel',{'U_1','U_2','U_3'},...
    'YTick',1:3,'yticklabel',{'\beta','\sigma','\gamma'})
title('data I','FontSize',14,'FontWeight','bold')

set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

subplot(1,4,2)

imagesc(1:3,1:3,abs(U13));

% 设置 colormap
cmap = othercolor('BuDRd_12');
colormap(cmap);  % 可以选择其他 colormap 例如 'jet', 'hot', 'cool' 等
clim([0,1.0])
% 添加 colorbar 并设置标签
c = colorbar;
c.Label.String = '|\partial U_i^T\theta/\partial \theta_j|';  % 设置 colorbar 的标签
c.FontSize = 15;  % 调整字体大小
c.Label.FontWeight = 'bold';  % 设置字体加粗
set(gca,'xticklabel',{'U_1','U_2','U_3'},...
    'YTick',1:3,'yticklabel',{'\beta','\sigma','\gamma'})
title('data I & S','FontSize',14,'FontWeight','bold')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

subplot(1,4,3)
imagesc(1:3,1:3,abs(U23));

% 设置 colormap
cmap = othercolor('BuDRd_12');
colormap(cmap);  % 可以选择其他 colormap 例如 'jet', 'hot', 'cool' 等
clim([0,1.0])
% 添加 colorbar 并设置标签
c = colorbar;
c.Label.String = '|\partial U_i^T\theta/\partial \theta_j|';  % 设置 colorbar 的标签
c.FontSize = 15;  % 调整字体大小
c.Label.FontWeight = 'bold';  % 设置字体加粗
set(gca,'xticklabel',{'U_1','U_2','U_3'},...
    'YTick',1:3,'yticklabel',{'\beta','\sigma','\gamma'})
title('data I & E','FontSize',14,'FontWeight','bold')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

subplot(1,4,4)
imagesc(1:3,1:3,abs(U123));

% 设置 colormap
cmap = othercolor('BuDRd_12');
colormap(cmap);  % 可以选择其他 colormap 例如 'jet', 'hot', 'cool' 等
clim([0,1.0])
% 添加 colorbar 并设置标签
c = colorbar;
c.Label.String = '|\partial U_i^T\theta/\partial \theta_j|';  % 设置 colorbar 的标签
c.FontSize = 15;  % 调整字体大小
c.Label.FontWeight = 'bold';  % 设置字体加粗
set(gca,'xticklabel',{'U_1','U_2','U_3'},...
    'YTick',1:3,'yticklabel',{'\beta','\sigma','\gamma'})
title('data I & E & S','FontSize',14,'FontWeight','bold')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

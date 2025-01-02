% File: Michaelis_Menten.m

% Parameters
para = zeros(3,1);
para(1) = 1e6;   % reaction rate: k1 = para(1)
para(2) = 1e-4;    % reaction rate: k2 = para(2)
para(3) = 0.1;    % reaction rate: k3 = para(3)

% Initial conditions for prey (x0) and predator (y0)
y10 = 5e-7;  % Initial y1
y20 = 2e-7;  % Initial y2
y30 = 0.0;  % Initial y3
y40 = 0.0;  % Initial y4
yzero = [y10 y20 y30 y40];
N = 6;
C = y10+y20+2*y30+y40;
para(1) = para(1)*C;
yzero = yzero/C;
% Time span
tspan = linspace(0,60,N);

% Solve the ODEs
[t, Y] = ode15s(@(t,y)Michaelis_Menten_eq(t,y,para), tspan, yzero,para);

% Extract data
y1 = Y(:, 1);
y4 = Y(:, 4);
yexp1 = y1;
yexp2 = y4;
 %%
% options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',100000,'MaxIterations',20000,...
%     'algorithm','levenberg-marquardt','StepTolerance',1e-8);
% 

% lb = [0 0 0];
% ub = [1e7*C 1e-2 1.0];
% para_new0 = para;
% 
% [para_new,resnorm,residual,exitflag,output,lambda,jacobian] = ...
%     lsqnonlin(@Objective_mm,para_new0,lb,ub,options,tspan,y0,1,yexp1,yexp2);       
%% Fisher information matrix
% tspan1=tspan;
para_new = para;
y0 = yzero';
x0 = zeros(4,1);
x0 =[x0;y0];
[~, Xi3] = ode15s(@(t,x)partial_theta3(t, x,para_new), tspan,x0);
[~, Xi2] = ode15s(@(t,x)partial_theta2(t, x,para_new), tspan,x0);
[~, Xi1] = ode15s(@(t,x)partial_theta1(t, x,para_new), tspan,x0);

Patial_phi1_theta = [Xi1(:,1) Xi2(:,1) Xi3(:,1)];
% Patial_phi2_theta = [Xi1(:,2) Xi2(:,2) Xi3(:,2)];
% Patial_phi3_theta = [Xi1(:,3) Xi2(:,3) Xi3(:,3)];
Patial_phi4_theta = [Xi1(:,4) Xi2(:,4) Xi3(:,4)];
S = [Patial_phi1_theta;Patial_phi4_theta]; %Sensitive Matrix
%S = Patial_phi1_theta;
%S = Patial_phi2_theta;
%S = Patial_phi3_theta;
%S = Patial_phi4_theta;
 F = S'* S;
 [U,Sigma,~]=svd(F);
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

%% UQ

fig1=figure(1);
clf()
set(gcf,"Position",[650,481,578,323])

N1=101;
tspan = linspace(0,60,N1);
tspan1=tspan;
x0 = zeros(4,1);
x0 =[x0;y0];
[~, Xi3] = ode15s(@(t,x)partial_theta3(t, x,para_new), tspan,x0);
[~, Xi2] = ode15s(@(t,x)partial_theta2(t, x,para_new), tspan,x0);
[~, Xi1] = ode15s(@(t,x)partial_theta1(t, x,para_new), tspan,x0);

Patial_phi1_theta = [Xi1(:,1) Xi2(:,1) Xi3(:,1)];
Patial_phi2_theta = [Xi1(:,2) Xi2(:,2) Xi3(:,2)];
Patial_phi3_theta = [Xi1(:,3) Xi2(:,3) Xi3(:,3)];
Patial_phi4_theta = [Xi1(:,4) Xi2(:,4) Xi3(:,4)];
%S = [Patial_phi1_theta;Patial_phi2_theta;Patial_phi3_theta;Patial_phi4_theta]; %Sensitive Matrix
%S = [Patial_phi1_theta;Patial_phi2_theta]; %Sensitive Matrix
[t, phi] = ode15s(@(t,y)Michaelis_Menten_eq(t,y,para), tspan, y0,para);
% F = S'* S;
% [U,Sigma,~]=svd(F);
r = min(find(diag(Sigma)<1e-6));
Ur = U(:,r:end);
%Ur = U;
S = Patial_phi1_theta;
sigma_th = 5e-1; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (S(i,:) * Ur)*(S(i,:) * Ur)';
end

% x = tspan;
% Y_fit = phi(:,1);
% CI=1.96*sqrt(Var_NN');
% % 绘制置信区域（拟合曲线上下的置信区间）
% fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 
% hold on

S = Patial_phi4_theta;
%Ur = U;
%sigma_th = 2e5; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (S(i,:) * Ur)*(S(i,:) * Ur)';
end
% x = tspan;
% Y_fit = phi(:,4);
% CI=1.96*sqrt(Var_NN');
% fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% 
% Plot results

plot(t, phi(:,1), 'r', 'LineWidth', 2,'LineStyle','-'); 
hold on
plot(t, phi(:,4), 'b', 'LineWidth', 2);
hold on
tspan = [0,10,20,50];%[0,13.75,30.75,48];
[t, phi] = ode15s(@(t,y)Michaelis_Menten_eq(t,y,para), tspan, y0,para);

plot(t(2:end), phi(2:end,4), 'marker','square','MarkerSize',12,'color','k', 'LineStyle','none','LineWidth',2);
hold on
t = linspace(0,60,N);
plot(t(2:end), y1(2:end), 'ko', 'LineWidth', 2,'MarkerSize',10); 
hold on;
plot(t(2:end),y4(2:end), 'ko', 'LineWidth', 2,'MarkerSize',10);



xlabel('Time');
ylabel('Concentration');

lgd = legend("S simu.","P simu.",'critical data',"synthetic data");
lgd.ItemTokenSize = [10,6];
lgd.Position = [0.746 0.301 0.139 0.326];
lgd.Box = 'off';
lgd.FontWeight = 'bold';
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)



fig3=figure(3);
clf()
set(gcf,"Position",[650,481,578,323])


Ur = U;
S = Patial_phi1_theta;
sigma_th = 5e-4; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (S(i,:) * Ur)*(S(i,:) * Ur)';
end
tspan = linspace(0,60,N1);
[t, phi] = ode15s(@(t,y)Michaelis_Menten_eq(t,y,para), tspan, y0,para);

x = tspan;
Y_fit = phi(:,1);
CI=1.96*sqrt(Var_NN');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

hold on

S = Patial_phi4_theta;
Ur = U;
%sigma_th = 1e-3; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (S(i,:) * Ur)*(S(i,:) * Ur)';
end
x = tspan;
Y_fit = phi(:,4);
CI=1.96*sqrt(Var_NN');
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% 
% Plot results
t = linspace(0,60,N1);
plot(t, phi(:,1), 'r', 'LineWidth', 2); 
hold on
plot(t, phi(:,4), 'b', 'LineWidth', 2);
t = linspace(0,60,N);
plot(t(2:end), y1(2:end), 'ko', 'LineWidth', 2,'MarkerSize',10); 
hold on;
plot(t(2:end),y4(2:end), 'ko', 'LineWidth', 2,'MarkerSize',10);

xlabel('Time');
ylabel('Concentration');

lgd = legend("95% CI S","95% CI P","S simu.","P simu.","synthetic data");
lgd.ItemTokenSize = [10,6];
lgd.Position = [0.746 0.301 0.139 0.326];
lgd.Box = 'off';
lgd.FontWeight = 'bold';
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)


%%
fig2=figure(2);
clf();
set(gcf,'Position',[479,231,745,525])


subplot(2,2,3)
alphaData = ones(3,3);  % 初始化为全不透明
alphaData(:, r:3) = 0.2;  % 设置右半部分透明度为 0.2

imagesc(1:3,1:3,abs(U),'AlphaData',alphaData);

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
    'YTick',1:3,'XTick',1:3,'yticklabel',{'k_1','k_2','k_3'})
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
title('synthetic data','FontSize',15,'FontWeight','bold')
box off
% title('\theta=[k_1,k_2,k_3]','FontSize',14,'FontWeight','bold')





%% Fisher information matrix
tspan2 = [0,10,20,50];
%tspan2 =[0,20,40,60];
tspan = tspan2;
x0 = zeros(4,1);
x0 =[x0;y0];
[~, Xi3] = ode15s(@(t,x)partial_theta3(t, x,para_new), tspan,x0);
[~, Xi2] = ode15s(@(t,x)partial_theta2(t, x,para_new), tspan,x0);
[~, Xi1] = ode15s(@(t,x)partial_theta1(t, x,para_new), tspan,x0);

Patial_phi4_theta = [Xi1(:,4) Xi2(:,4) Xi3(:,4)];
S1 = Patial_phi4_theta; %Sensitive Matrix
%S = Patial_phi1_theta;
%S = Patial_phi2_theta;
F = S1'* S1;

% F = S'* S;
[U1,Sigma1,~]=svd(F);
r = min(find(diag(Sigma)<1e-6));

%% RAU
nn = size(S1,1);
s1=S1(:,1);
A = S1(:,2:3);
ss1 = (eye(nn)-A*pinv(A))*s1;

s2=S1(:,2);
A = [S1(:,1) S1(:,3)];
ss2 = (eye(nn)-A*pinv(A))*s2;

s3=S1(:,3);
A = [S1(:,1) S1(:,2)];
ss3 = (eye(nn)-A*pinv(A))*s3;
SSS=[SSS;norm(ss1,"inf") norm(ss2,"inf") norm(ss3,"inf")];
%Ur = U;
subplot(2,2,2)
x=1:3;
bar(x,SSS)
ylim([1e-4,10])
ylabel('$$\|(I - AA^{\dagger}) \mathbf{s}_i\|_{\infty}$$', 'Interpreter', 'latex');
set(gca,'xticklabel',{'k_1','k_2','k_3'},'YScale','log')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
subplot(2,2,4)
alphaData = ones(3,3);  % 初始化为全不透明
alphaData(:, r:3) = 0.2;  % 设置右半部分透明度为 0.2

imagesc(1:3,1:3,abs(U1),'AlphaData',alphaData);

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
    'YTick',1:3,'XTick',1:3,'yticklabel',{'k_1','k_2','k_3'})
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
title('critical data','FontSize',15,'FontWeight','bold')
box off

subplot(2,2,1)
x=1:3;
bar(x,[(Sigma(1,1)),(Sigma(2,2)),(Sigma(3,3));...
    (Sigma1(1,1)),(Sigma1(2,2)),(Sigma1(3,3))])
hold on
plot([0,5],[1e-6,1e-6],'LineStyle','--','LineWidth',1.5,'Color','k')
ylim([1e-8,1e3])
xlim([0,4])
ylabel('Eigenvalue of FIM');
set(gca,'xticklabel',{'U_1^T\theta','U_2^T\theta','U_3^T\theta'},'YScale','log',...
    'ytick',[1e-6,1e-4,1e-2,1,1e2])
lgd = legend("synthetic data",'critical data');
lgd.ItemTokenSize = [10,6];

lgd.Box = 'off';
lgd.FontWeight = 'bold';
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

box off
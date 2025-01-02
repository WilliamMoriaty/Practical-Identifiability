% File: lotka_volterra.m

% Parameters
para = zeros(3,1);
para(1) = 5.4;   % beta = para(1)
para(2) = 5.8;    % sigma = para(2)
para(3) = 4.54;    % gamma = para(3)

data = readmatrix('Influenza_data.csv');

% Initial conditions for prey (x0) and predator (y0)
y10 = 157756/157759;  % Initial y1
y20 = 0.0;  % Initial y2
y30 = 3.0/157759;  % Initial y3
y40 = 0.0;  % Initial y4
yzero = [y10 y20 y30 y40];
N = 21;
% Time span
tspan = data(:,1);

% Solve the ODEs
[t, Y] = ode45(@(t,y)SEIR_eq(t,y,para), tspan, yzero,para);

% % Extract the prey and predator populations
% y2 = Y(:, 2);
% y3 = Y(:, 3);
yexp1 = data(:,2)/157759;
yexp2 = [];
%% Fit parameter
options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',1e8,'MaxIterations',20000,...
    'algorithm','levenberg-marquardt','StepTolerance',1e-8);

y0 = yzero';
lb = [0 0 0];
ub = [10.0 10.0 10.0];
para_new0 = para;

[para_new,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqnonlin(@Objective_seir,para_new0,lb,ub,options,tspan,y0,1,yexp1,yexp2,[],[],[]);       
% N=501;
% tspan1=linspace(0,48,N);
% [t, phi] = ode45(@(t,y)SEIR_eq(t,y,para_new), tspan1, y0,para);
% plot(t,phi(:,3),'r-')
% hold on
% plot(data(:,1),data(:,2),'ro')
% hold on
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
%% RAU
nn = size(S,1);
s1=S(:,1);
A = [S(:,2) S(:,3)];
ss1 = (eye(nn)-A*pinv(A))*s1;

s2=S(:,2);
A = [S(:,1) S(:,3)];
ss2 = (eye(nn)-A*pinv(A))*s2;

s3=S(:,3);
A = [S(:,1) S(:,2)];
ss3 = (eye(nn)-A*pinv(A))*s3;
SSS=[norm(ss1,"inf") norm(ss2,"inf") norm(ss3,"inf")];
r = min(find(diag(Sigma)<1e-6));
Ur = U(:,r:end);
%% normalization
[para_new1,resnorm1,residual1,exitflag1,output1,lambda1,jacobian1] = ...
    lsqnonlin(@Objective_seir,para_new,lb,ub,options,tspan,y0,1,yexp1,yexp2,0.1,Ur,para_new);       
%%
fig1=figure(1);
clf()
set(gcf,"Position",[-67,463,1620,249])
subplot(1,4,1)
N1=201;
tspan = linspace(0,48,N1);
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

%Ur = U;
S = Patial_phi3_theta;
sigma_th = 1e2; % parameter estimation
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

plot(data(:,1),data(:,2)/157759,'ro')
hold on
% plot(t(1:10:end), phi(1:10:end,4), 'b*', 'LineWidth', 2);

xlabel('Weeks');
ylabel('Ratio');
xlim([0,50])
ylim([0,0.01])
% lg = legend("95% CI I","95% CI R","I","R","data I","data R");
lgd = legend("95% CI","I","Influenza data");
lgd.ItemTokenSize = [10,6];
% lgd.Position = [0.746 0.301 0.139 0.326];
lgd.Box = 'off';
lgd.FontWeight = 'bold';
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)


%%

subplot(1,4,2)
x=1:3;
bar(x,SSS,'FaceColor',[0,0,0])
ylim([1e-6,1e-3])
ylabel('$$\|(I - AA^{\dagger}) \mathbf{s}_i\|_{\infty}$$', 'Interpreter', 'latex');
set(gca,'xticklabel',{'\beta','\sigma','\gamma'},'YScale','log')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
subplot(1,4,3)
x=1:3;
bar(x,[Sigma(1,1),Sigma(2,2),Sigma(3,3)],'FaceColor',[0,0,0])
hold on
plot([0,4],[1e-6,1e-6],'LineStyle','--','LineWidth',1.5,'Color','k')

ylabel('Eigenvalue of FIM');
set(gca,'xticklabel',{'U_1^T\theta','U_2^T\theta','U_3^T\theta'},'YScale','log')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
subplot(1,4,4)
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
c.FontSize = 12;  % 调整字体大小
c.Label.FontWeight = 'bold';  % 设置字体加粗
set(gca,'xticklabel',{'U_1','U_2','U_3'},...
    'YTick',1:3,'yticklabel',{'\beta','\sigma','\gamma'})
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
%% 
fig2=figure(2);
clf();


Ur = U;
S = Patial_phi3_theta;
sigma_th = 1e-3; % parameter estimation
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

% Plot results

 plot(t, phi(:,3), 'r', 'LineWidth', 2); 
% hold on
% plot(t, phi(:,4), 'b', 'LineWidth', 2);
plot(data(:,1),data(:,2)/157759,'ro')
hold on
% plot(t(1:10:end), phi(1:10:end,4), 'b*', 'LineWidth', 2);

xlabel('Weeks');
ylabel('Population');
lgd = legend("95% CI","I","Influenza data");
lgd.ItemTokenSize = [10,6];
% lgd.Position = [0.746 0.301 0.139 0.326];
lgd.Box = 'off';
lgd.FontWeight = 'bold';
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

% File: lotka_volterra.m

% Parameters
para = zeros(4,1);
para(1) = 0.545;   % Prey growth rate: alpha = para(1)
para(2) = 0.028;    % Predation rate: beta = para(2)
para(3) = 0.024; % Reproduction rate of predator: delta = para(3)
para(4) = 0.803;   % Predator death rate: gamma = para(4)

data = readmatrix('hudson-bay-lynx-hare.csv');

% Initial conditions for prey (x0) and predator (y0)
x0 = 33.956;  % Initial prey population
y0 = 5.933;   % Initial predator population
%N = 21;
% Time span
tspan = data(:,1);

% Solve the ODEs
%[t, Y] = ode45(@(t,y)lotka_volterra_eq(t,y,para), tspan, [x0 y0],para);

% Extract the prey and predator populations
% prey = Y(:, 1);
% predator = Y(:, 2);
yexp1 = data(:,2);
yexp2 = data(:,1);

%%
options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',1e8,'MaxIterations',20000,...
    'algorithm','levenberg-marquardt','StepTolerance',1e-6);

y0 = [x0;y0];
lb = [0 0 0 0];
ub = [1 0.05 0.05 1];
para_new0 = para;

[para_new,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqnonlin(@Objective_lv,para_new0,lb,ub,options,tspan,y0,1,yexp1,yexp2,[],[],[]);       
% N=501;
% tspan1=linspace(1900,1920,N);
% [t, phi] = ode45(@(t,y)lotka_volterra_eq(t,y,para_new), tspan1, y0,para);
% plot(t,phi(:,1),'r-')
% hold on
% plot(t,phi(:,2),'b-')
% hold on

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
r = min(find(diag(Sigma)<1e0));
Ur = U(:,r:end);
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
% File: lotka_volterra.m

%% normalization
[para_new1,resnorm1,residual1,exitflag1,output1,lambda1,jacobian1] = ...
    lsqnonlin(@Objective_lv,para_new,lb,ub,options,tspan,y0,2,yexp1,yexp2,0.1,Ur,para_new);  
N=501;
% tspan1=linspace(1900,1920,N);
% [t1, phi1] = ode45(@(t,y)lotka_volterra_eq(t,y,para_new1), tspan1, y0,para);

%% UQ
fig1=figure(1);
clf()
set(gcf,"Position",[650,481,578,323])
N1=1001;
tspan = linspace(1900,1920,N1);
tspan1=tspan;
x0 = [0;0;y0];
[~, Xi4] = ode45(@(t,x)partial_theta4(t, x,para_new),tspan,x0);
[~, Xi3] = ode45(@(t,x)partial_theta3(t, x,para_new),tspan,x0);
[~, Xi2] = ode45(@(t,x)partial_theta2(t, x,para_new),tspan,x0);
[~, Xi1] = ode45(@(t,x)partial_theta1(t, x,para_new),tspan,x0);

Patial_phi1_theta = [Xi1(:,1) Xi2(:,1) Xi3(:,1) Xi4(:,1)];
Patial_phi2_theta = [Xi1(:,2) Xi2(:,2) Xi3(:,2) Xi4(:,2)];
%S = [Patial_phi1_theta;Patial_phi2_theta]; %Sensitive Matrix
[t, phi] = ode45(@(t,y)lotka_volterra_eq(t,y,para_new1), tspan1, y0,para_new1);
% F = S'* S;
% [U,Sigma,~]=svd(F);


S = Patial_phi1_theta;
sigma_th = 1e-4; % parameter estimation
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
hold on
plot(data(:,1),data(:,3),'ro','MarkerSize',6,'LineWidth',1.2)
hold on
plot(data(:,1),data(:,2),'bo','MarkerSize',6,'LineWidth',1.2)
hold on
plot(t,phi(:,1),'r-','LineWidth',1.2)
hold on
plot(t,phi(:,2),'b-','LineWidth',1.2)
hold on

xlabel('Time');
ylabel('Population');
lgd = legend("95% CI hare","95% CI lynx","hare data","lynx data","hare","lynx");
lgd.ItemTokenSize = [10,6];
lgd.Box = 'off';
lgd.FontWeight = 'bold';
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

% 
fig3=figure(3);
clf()
set(gcf,"Position",[650,481,578,323])
% subplot(1,2,2)
Ur = U;
S = Patial_phi1_theta;
sigma_th = 1e-6; % parameter estimation
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

hold on
plot(data(:,1),data(:,3),'ro','MarkerSize',6,'LineWidth',1.2)
hold on
plot(data(:,1),data(:,2),'bo','MarkerSize',6,'LineWidth',1.2)
hold on
plot(t,phi(:,1),'r-','LineWidth',1.2)
hold on
plot(t,phi(:,2),'b-','LineWidth',1.2)
hold on

xlabel('Time');
ylabel('Population');
% lgd = legend("95% CI hare","95% CI lynx","hare data","lynx data","hare","lynx");
% lgd.ItemTokenSize = [10,6];
% lgd.Box = 'off';
% lgd.FontWeight = 'bold';
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)


fig2=figure(2);
clf();
set(gcf,'Position',[228,445,768,280])
% subplot(1,3,1)
% x=1:4;
% bar(x,SSS,'FaceColor',[0,0,0])
% ylim([10 15000])
% ylabel('$$\|(I - AA^{\dagger}) \mathbf{s}_i\|_{\infty}$$', 'Interpreter', 'latex');
% set(gca,'xticklabel',{'\alpha','\beta','\delta','\gamma'},'YScale','log')
% set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
% box off
subplot(1,2,1)
x=1:4;
bar(x,[(Sigma(1,1)),(Sigma(2,2)),(Sigma(3,3)),(Sigma(4,4))],'FaceColor',[0,0,0])
hold on
plot([0,5],[1e6,1e6],'LineStyle','--','Color','k','LineWidth',1.2)
ylim([1e3,1e9])
xlim([0,5])
ylabel('Eigenvalue of FIM');
set(gca,'xticklabel',{'U_1^T\theta','U_2^T\theta','U_3^T\theta','U_4^T\theta'},'YScale','log')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
subplot(1,2,2)
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
c.FontSize = 15;  % 调整字体大小
c.Label.FontWeight = 'bold';  % 设置字体加粗
set(gca,'xticklabel',{'U_1','U_2','U_3','U_4'},...
    'YTick',1:4,'yticklabel',{'\alpha','\beta','\delta','\gamma'})
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
% title('\theta=[\alpha,\beta,\delta,\gamma]','FontSize',14,'FontWeight','bold')


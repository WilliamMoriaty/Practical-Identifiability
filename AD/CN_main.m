% File: CN.m

% Parameters
para = zeros(9,1);
para(1) = 6.35*10^(-2);%lambda_Ab = para(1);
para(2) = 0.15;%lambda_tau = para(2);
para(3) = 2e-2;%lambda_Ntaup = para(3);
para(4) = 1.67e-2;%lambda_CN = para(4);
para(5) = 0.83e-5;%lambda_Ctau = para(5);
para(6) = 139.94;% K_Ab = para(6);
para(7) = 123.35;% K_taup = para(7);
para(8) = 1.00;% K_N = para(8);
para(9) = 50.48;% K_C = para(9);

% Initial conditions 

y10 = 44.92;  % Initial Ab
y20 = 3.68;  % Initial tau_p
y30 = 0.52;  % Initial N
y40 = 1.68; % Initial C
yzero = [y10 y20 y30 y40];
% N = 51;
% % Time span

% % Extract the prey and predator populations
load('CN_data/CN_677_13.mat')
load('CN_data/CN_677_45.mat')
[NAge, idx] = sort(NAge);
 N_PID = N_PID(idx);
 C_PID = C_PID(idx);
[NAge1, ~, idx] = unique(NAge);
maxnum=max(idx);
N_PID1=NAge1;
C_PID1=NAge1;

for i=1:maxnum
aa = find(idx == i);
N_PID1(i)=sum(N_PID(aa))/size(aa,1);
C_PID1(i)=sum(C_PID(aa))/size(aa,1);
end



tspan1 = [50;AAge];
tspan2 = [50;NAge1];
% y2 = Y(:, 2);
% y3 = Y(:, 3);
yexp1 = [yzero(1);Abeta_PID];
yexp2 = [yzero(2);ptau_PID];
yexp3 = [yzero(3);N_PID1];
yexp4 = [yzero(4);C_PID1];


%%
options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',1e6,'MaxIterations',20000,...
    'algorithm','levenberg-marquardt','StepTolerance',1e-8);

  y0 = yzero';
 lb = [0 0 0 0 0 0 0 0 0];
 ub = [10.0 1.0 10.0 10 10 500 500 2 500];
para_new0 = para;

[para_new,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqnonlin(@Objective_AD,para_new0,lb,ub,options,tspan1,tspan2,y0,yexp1,yexp2,yexp3,yexp4,0,zeros(9,1),para);       
% N=101;
% tspan = linspace(50,100,N);
% 
% % Solve the ODEs
% [t, Y] = ode45(@(t,y)AD_eq(t,y,para_new), tspan, yzero);
% fig1=figure(1);
% clf();
% subplot(1,5,1)
% plot(t,Y(:,1))
% hold on
% plot(tspan1(2:end),yexp1(2:end),'ro')
% xlim([50,100])
% xlabel('Age')
% ylabel('A_\beta')
% subplot(1,5,2)
% plot(t,Y(:,2))
% hold on
% plot(tspan1(2:end),yexp2(2:end),'ro')
% xlim([50,100])
% xlabel('Age')
% ylabel('\tau_p')
% subplot(1,5,3)
% plot(t,Y(:,3))
% hold on
% plot(tspan1(2:end),yexp3(2:end),'ro')
% xlim([50,100])
% xlabel('Age')
% ylabel('\tau_0')
% subplot(1,5,4)
% plot(t,Y(:,4))
% hold on
% plot(tspan2(2:end),yexp4(2:end),'ro')
% xlim([50,100])
% xlabel('Age')
% ylabel('N')
% subplot(1,5,5)
% plot(t,Y(:,5))
% hold on
% plot(tspan2(2:end),yexp5(2:end),'ro')
% xlim([50,100])
% xlabel('Age')
% ylabel('C')

%% Fisher informatino matrix

x0 = zeros(4,1);
x0 = [x0;y0];
tspan = tspan1;

[~, Xi9] = ode45(@(t,x)partial_theta9(t, x,para_new), tspan,x0);
[~, Xi8] = ode45(@(t,x)partial_theta8(t, x,para_new), tspan,x0);
[~, Xi7] = ode45(@(t,x)partial_theta7(t, x,para_new), tspan,x0);
[~, Xi6] = ode45(@(t,x)partial_theta6(t, x,para_new), tspan,x0);
[~, Xi5] = ode45(@(t,x)partial_theta5(t, x,para_new), tspan,x0);
[~, Xi4] = ode45(@(t,x)partial_theta4(t, x,para_new), tspan,x0);
[~, Xi3] = ode45(@(t,x)partial_theta3(t, x,para_new), tspan,x0);
[~, Xi2] = ode45(@(t,x)partial_theta2(t, x,para_new), tspan,x0);
[~, Xi1] = ode45(@(t,x)partial_theta1(t, x,para_new), tspan,x0);

Patial_phi1_theta = [Xi1(:,1) Xi2(:,1) Xi3(:,1) Xi4(:,1) Xi5(:,1) Xi6(:,1) Xi7(:,1) Xi8(:,1) Xi9(:,1)];
Patial_phi2_theta = [Xi1(:,2) Xi2(:,2) Xi3(:,2) Xi4(:,2) Xi5(:,2) Xi6(:,2) Xi7(:,2) Xi8(:,2) Xi9(:,2)];

S = [Patial_phi1_theta;Patial_phi2_theta]; %Sensitive Matrix

tspan = tspan2;

[~, Xi9] = ode45(@(t,x)partial_theta9(t, x,para_new), tspan,x0);
[~, Xi8] = ode45(@(t,x)partial_theta8(t, x,para_new), tspan,x0);
[~, Xi7] = ode45(@(t,x)partial_theta7(t, x,para_new), tspan,x0);
[~, Xi6] = ode45(@(t,x)partial_theta6(t, x,para_new), tspan,x0);
[~, Xi5] = ode45(@(t,x)partial_theta5(t, x,para_new), tspan,x0);
[~, Xi4] = ode45(@(t,x)partial_theta4(t, x,para_new), tspan,x0);
[~, Xi3] = ode45(@(t,x)partial_theta3(t, x,para_new), tspan,x0);
[~, Xi2] = ode45(@(t,x)partial_theta2(t, x,para_new), tspan,x0);
[~, Xi1] = ode45(@(t,x)partial_theta1(t, x,para_new), tspan,x0);

Patial_phi3_theta = [Xi1(:,3) Xi2(:,3) Xi3(:,3) Xi4(:,3) Xi5(:,3) Xi6(:,3) Xi7(:,3) Xi8(:,3) Xi9(:,3)];
Patial_phi4_theta = [Xi1(:,4) Xi2(:,4) Xi3(:,4) Xi4(:,4) Xi5(:,4) Xi6(:,4) Xi7(:,4) Xi8(:,4) Xi9(:,4)];

S = [S;Patial_phi3_theta;Patial_phi4_theta]; %Sensitive Matrix


 F = S'* S;
 [U,Sigma,~]=svd(F);
 thresh = 1e-3;
r = min(find(diag(Sigma)<thresh));
Ur = U(:,r:end);
 %% RAU
nn = size(S,1);
s1=S(:,1);
A = [S(:,9) S(:,2:8)];
ss1 = (eye(nn)-A*pinv(A))*s1;

s2=S(:,2);
A = [S(:,1) S(:,9) S(:,3:8)];
ss2 = (eye(nn)-A*pinv(A))*s2;

s3=S(:,3);
A = [S(:,1:2) S(:,9) S(:,4:8)];
ss3 = (eye(nn)-A*pinv(A))*s3;

s4=S(:,4);
A = [S(:,1:3) S(:,9) S(:,5:8)];
ss4 = (eye(nn)-A*pinv(A))*s4;

s5=S(:,5);
A = [S(:,1:4) S(:,9) S(:,6:8)];
ss5 = (eye(nn)-A*pinv(A))*s5;

s6=S(:,6);
A = [S(:,1:5) S(:,9) S(:,7:8)];
ss6 = (eye(nn)-A*pinv(A))*s6;

s7=S(:,7);
A = [S(:,1:6) S(:,9) S(:,8:8)];
ss7 = (eye(nn)-A*pinv(A))*s7;

s8=S(:,8);
A = [S(:,1:7) S(:,9:9)];
ss8 = (eye(nn)-A*pinv(A))*s8;

s9=S(:,9);
A = S(:,1:8);
ss9 = (eye(nn)-A*pinv(A))*s9;



SSS=[norm(ss1,"inf") norm(ss2,"inf") norm(ss3,"inf"),...
    norm(ss4,"inf") norm(ss5,"inf") norm(ss6,"inf"),...
    norm(ss7,"inf") norm(ss8,"inf") norm(ss9,"inf")];


%% normalization



[para_new1,resnorm1,residual1,exitflag1,output1,lambda1,jacobian1] = ...
    lsqnonlin(@Objective_AD,para_new0,lb,ub,options,tspan1,tspan2,y0,yexp1,yexp2,yexp3,yexp4,100,Ur,para_new);       
% 
% N=101;
% tspan = linspace(50,100,N);
% 
% % Solve the ODEs
% [t, Y1] = ode45(@(t,y)AD_eq(t,y,para_new), tspan, yzero);
% [~, Y2] = ode45(@(t,y)AD_eq(t,y,para_new1), tspan, yzero);
% fig1=figure(1);
% clf();
% subplot(1,5,1)
% plot(t,Y1(:,1))
% hold on
% plot(t,Y2(:,1),'r-')
% hold on
% plot(tspan1(2:end),yexp1(2:end),'ro')
% xlim([50,100])
% xlabel('Age')
% ylabel('A_\beta')
% subplot(1,5,2)
% plot(t,Y1(:,2))
% hold on
% plot(t,Y2(:,2),'r-')
% hold on
% plot(tspan1(2:end),yexp2(2:end),'ro')
% xlim([50,100])
% xlabel('Age')
% ylabel('\tau_p')
% subplot(1,5,3)
% plot(t,Y1(:,3))
% hold on
% plot(t,Y2(:,3),'r-')
% hold on
% plot(tspan1(2:end),yexp3(2:end),'ro')
% xlim([50,100])
% xlabel('Age')
% ylabel('\tau_0')
% subplot(1,5,4)
% plot(t,Y1(:,4))
% hold on
% plot(t,Y2(:,4),'r-')
% hold on
% plot(tspan2(2:end),yexp4(2:end),'ro')
% xlim([50,100])
% xlabel('Age')
% ylabel('N')
% subplot(1,5,5)
% plot(t,Y1(:,5))
% hold on
% plot(t,Y2(:,5),'r-')
% hold on
% plot(tspan2(2:end),yexp5(2:end),'ro')
% xlim([50,100])
% xlabel('Age')
% ylabel('C')

%% UQ
fig3=figure(3);
clf()
set(gcf,'Position',[112,270,918,200])

N1=201;
tspan = linspace(50,100,N1);

x0 = zeros(4,1);
x0 = [x0;y0];

[~, Xi9] = ode45(@(t,x)partial_theta9(t, x,para_new1), tspan,x0);
[~, Xi8] = ode45(@(t,x)partial_theta8(t, x,para_new1), tspan,x0);
[~, Xi7] = ode45(@(t,x)partial_theta7(t, x,para_new1), tspan,x0);
[~, Xi6] = ode45(@(t,x)partial_theta6(t, x,para_new1), tspan,x0);
[~, Xi5] = ode45(@(t,x)partial_theta5(t, x,para_new1), tspan,x0);
[~, Xi4] = ode45(@(t,x)partial_theta4(t, x,para_new1), tspan,x0);
[~, Xi3] = ode45(@(t,x)partial_theta3(t, x,para_new1), tspan,x0);
[~, Xi2] = ode45(@(t,x)partial_theta2(t, x,para_new1), tspan,x0);
[~, Xi1] = ode45(@(t,x)partial_theta1(t, x,para_new1), tspan,x0);

Patial_phi1_theta = [Xi1(:,1) Xi2(:,1) Xi3(:,1) Xi4(:,1) Xi5(:,1) Xi6(:,1) Xi7(:,1) Xi8(:,1) Xi9(:,1)];
Patial_phi2_theta = [Xi1(:,2) Xi2(:,2) Xi3(:,2) Xi4(:,2) Xi5(:,2) Xi6(:,2) Xi7(:,2) Xi8(:,2) Xi9(:,2)];
Patial_phi3_theta = [Xi1(:,3) Xi2(:,3) Xi3(:,3) Xi4(:,3) Xi5(:,3) Xi6(:,3) Xi7(:,3) Xi8(:,3) Xi9(:,3)];
Patial_phi4_theta = [Xi1(:,4) Xi2(:,4) Xi3(:,4) Xi4(:,4) Xi5(:,4) Xi6(:,4) Xi7(:,4) Xi8(:,4) Xi9(:,4)];

SS = [Patial_phi1_theta;Patial_phi2_theta;Patial_phi3_theta;Patial_phi4_theta]; %Sensitive Matrix
[t, phi] = ode45(@(t,y)AD_eq(t,y,para_new1), tspan, yzero);


%Ur = U;
subplot(1,3,1)
SS1 = Patial_phi1_theta;
sigma_th = 1e-4; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (SS1(i,:) * Ur)*(SS1(i,:) * Ur)';
end

x = tspan;
Y_fit = phi(:,1);
CI=1.96*sqrt(Var_NN');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on

plot(t,phi(:,1),'k-','LineWidth',1.5)
hold on
plot(tspan1(2:end),yexp1(2:end),'ko','LineWidth',1.5,'MarkerSize',8)
xlim([50,100])
xlabel('Age')
ylabel('A_\beta')

lgd = legend("95% CI","simu","data");
lgd.ItemTokenSize = [10,6];
% lgd.Position = [0.746 0.301 0.139 0.326];
lgd.Location = 'southeast';
lgd.Box = 'off';
lgd.FontWeight = 'bold';
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

subplot(1,3,2)
SS2 = Patial_phi2_theta;
% sigma_th = 1e1; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= 1* (SS2(i,:) * Ur)*(SS2(i,:) * Ur)';
end

x = tspan;
Y_fit = phi(:,2);
CI=1.96*sqrt(Var_NN');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on

plot(t,phi(:,2),'k-','LineWidth',1.5)
hold on
plot(tspan1(2:end),yexp2(2:end),'ko','LineWidth',1.5,'MarkerSize',8)
xlim([50,100])
xlabel('Age')
ylabel('\tau_p')

box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

subplot(1,3,3)

SS4 = Patial_phi4_theta;
% sigma_th = 1e1; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (SS4(i,:) * Ur)*(SS4(i,:) * Ur)';
end

x = tspan;
Y_fit = phi(:,4);
CI=1.96*sqrt(Var_NN');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(t,phi(:,4),'k-','LineWidth',1.5,'MarkerSize',8)
hold on
plot(tspan2(2:end),yexp4(2:end),'ko','LineWidth',1.5,'MarkerSize',8)
xlim([50,100])
xlabel('Age')
ylabel('C')

box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
%%
fig2=figure(2);
clf();
set(gcf,"Position",[340,97,1097,673])


subplot(2,2,3)
x=1:9;

bar(x,[(Sigma(1,1)),(Sigma(2,2)),(Sigma(3,3)),...
    (Sigma(4,4)),(Sigma(5,5)),(Sigma(6,6)),...
    (Sigma(7,7)),(Sigma(8,8)),(Sigma(9,9))],'FaceColor',[0 0 0])
hold on
plot([0 10],[thresh thresh],'k--','LineWidth',1.5)

ylim([1e-8,1e10])
ylabel('Eigenvalue of FIM');
set(gca,'xticklabel',{'U_1^T\theta','U_2^T\theta','U_3^T\theta',...
    'U_4^T\theta','U_5^T\theta','U_6^T\theta','U_7^T\theta','U_8^T\theta','U_9^T\theta'},'YScale','log')
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

subplot(2,2,4)
alphaData = ones(9,9);  % 初始化为全不透明
alphaData(:, r:9) = 0.2;  % 设置右半部分透明度为 0.2

imagesc(1:9,1:9,abs(U),'AlphaData', alphaData);

% 设置 colormap
cmap = othercolor('BuDRd_12');
colormap(cmap);  % 可以选择其他 colormap 例如 'jet', 'hot', 'cool' 等
clim([0,1.0])
% 添加 colorbar 并设置标签
c = colorbar;
c.Label.String = '|\partial U_i^T\theta/\partial \theta_j|';  % 设置 colorbar 的标签
c.FontSize = 15;  % 调整字体大小
c.Label.FontWeight = 'bold';  % 设置字体加粗
set(gca,'XTick',1:9,'xticklabel',{'U_1','U_2','U_3','U_4','U_5','U_6','U_7','U_8','U_9'},...
    'YTick',1:9,'yticklabel',{'\lambda_{A\beta}','\lambda_\tau',...
    '\lambda_{N\tau_p}','\lambda_{CN}','\lambda_{C\tau}',...
    'K_{A\beta}','K_{\tau_p}','K_N','K_C'})
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)


subplot(2,2,2)
x=1:9;
bar(x,SSS,'FaceColor',[0,0,0])
ylim([0.0001,1000])
ylabel('$$\|(I - AA^{\dagger}) \mathbf{s}_i\|_{\infty}$$', 'Interpreter', 'latex');
set(gca,'XTick',1:9,'xticklabel',{'\lambda_{A\beta}','\lambda_\tau'...
    ,'\lambda_{N\tau_p}','\lambda_{CN}','\lambda_{C\tau}',...
    'K_{A\beta}','K_{\tau_p}','K_N','K_C'},'YScale','log')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off

subplot(2,2,1)
SS3 = Patial_phi3_theta;
% sigma_th = 1e1; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (SS3(i,:) * Ur)*(SS3(i,:) * Ur)';
end

x = tspan;
Y_fit = phi(:,3);
CI=1.96*sqrt(Var_NN');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(t,phi(:,3),'k-','LineWidth',1.5)
hold on
plot(tspan2(2:end),yexp3(2:end),'ko','LineWidth',1.5,'MarkerSize',8)
xlim([50,100])
xlabel('Age')
ylabel('N')

box off
lgd = legend("95% CI","simu","data");
lgd.ItemTokenSize = [10,6];
% lgd.Position = [0.746 0.301 0.139 0.326];
lgd.Location = 'southeast';
lgd.Box = 'off';
lgd.FontWeight = 'bold';
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

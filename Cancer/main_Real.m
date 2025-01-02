% File: PDE.m

% Parameters
para = zeros(10,1);
s = 1.36e4*10;   % reaction rate: s = para(1)
f = 0.2998e8*5;    % reaction rate: f = para(2)
g = 2.02e7*3;    % reaction rate: g = para(3)
d1 = 0.0412;   % reaction rate: d1 = para(4)
k1 = 1.3e-7;    % reaction rate: k1 = para(5)
k_1 = 24.0;    % reaction rate: k_1 = para(6)
k2 = 7.2;   % reaction rate: k2 = para(7)
p = 0.9997;    % reaction rate: p = para(8)
b1 = 0.18;    % reaction rate: b1 = para(9)
b2 = 2.5e-7*(1/5);    % reaction rate: b2 = para(10)
% Nondimensional parameters
 E0=3.3e5*6;  % Effector T cells
 T0=4.0e6*5;   % Tumor cells
 C0=3.3e5*6; % Cell-complex 
 D1 = 1e-6;t0=1/D1;
para(1) = s*t0/E0;   % reaction rate: sigma = s*t0/E0
para(2) = f*t0*C0/(E0*T0);    % reaction rate: rho = f*t0*(C0/(E0*T0))
para(3) = g/T0;    % reaction rate: eta = g/T0
para(4) = k1*T0*t0;   % reaction rate: mu = k1*T0*t0
para(5) = t0*C0*(k_1+k2*p)/E0;    % reaction rate: eps = t0*C0*(k_1+k2*p)/E0
para(6) = b1*t0;    % reaction rate: beta1 = b1*t0
para(7) = b2*T0;   % reaction rate: beta2 = b2*T0
para(8) = k1*t0*E0;    % reaction rate: phi = k1*t0*E0
para(9) = t0*C0*(k_1+k2*(1-p))/T0;    % reaction rate: lambda = t0*C0*(k_1+k2*(1-p))/T0
para(10) = t0*(k_1+k2);    % reaction rate: psi = t0*(k_1+k2)


% Initial conditions
nn=128;
disx = linspace(0,1,nn+1);
l = 0.2*1;
y10=zeros(nn+1,1);
y20=zeros(nn+1,1);
y30=zeros(nn+1,1);
for i=1:nn+1
y10(i) = 0.001*hvd(disx(i)).*(1-exp(-1000*(disx(i)-l).^2));  % Initial E
y20(i) = 0.01*(1-hvd(disx(i)).*(1-exp(-1000*(disx(i)-l).^2)));  % Initial T
y30(i) = 0.001*(exp(-1000*(disx(i)-l).^2));  % Initial C

end

yzero = [y10;y20;y30];

% Time span
%tspan = linspace(0,1000*D1,101);

% Solve the ODEs
% [t, Y] = ode15s(@(t,y)PDE_eq(t,y,para,nn), tspan, yzero,para);
% sol = ode15s(@(t,y)PDE_eq(t,y,para,nn),tspan,y0);
% solpts = deval(sol,tspan);
% Extract the data of E and T
TCell=load('Glioma_tcell.txt');
TumorCell=load('Glioma_tumor.txt');

tspan2 = [0;TumorCell(:,1)]*D1;
tspan1 = [0;TCell(:,1)]*D1;

 yexp2=[sum(y20)/nn;TumorCell(:,2)/T0];
 yexp1=[sum(y10)/nn;TCell(:,2)/E0];
 
%% 
 options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',1e6,'MaxIterations',20000,...
    'algorithm','levenberg-marquardt','StepTolerance',1e-8);

  y0 = yzero;
 lb = [0 0 0 0 0 0 0 0 0 0];
 ub = 2*para';
para_new0 = para;

[para_new,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqnonlin(@Objective_PDE,para_new0,lb,ub,options,tspan1,tspan2,y0,yexp1,yexp2,0,zeros(10,1),para,nn);       


%% Fisher information
%tspan = tspan1;
% y0 = yzero;
% x0 = zeros((nn+1)*3,1);
% x0 = [x0;y0];
% [~, Xi10] = ode15s(@(t,x)partial_theta10(t, x,para_new,nn), tspan,x0);
% [~, Xi9] = ode15s(@(t,x)partial_theta9(t, x,para_new,nn), tspan,x0);
% [~, Xi8] = ode15s(@(t,x)partial_theta8(t, x,para_new,nn), tspan,x0);
% [~, Xi7] = ode15s(@(t,x)partial_theta7(t, x,para_new,nn), tspan,x0);
% [~, Xi6] = ode15s(@(t,x)partial_theta6(t, x,para_new,nn), tspan,x0);
% [~, Xi5] = ode15s(@(t,x)partial_theta5(t, x,para_new,nn), tspan,x0);
% [~, Xi4] = ode15s(@(t,x)partial_theta4(t, x,para_new,nn), tspan,x0);
% [~, Xi3] = ode15s(@(t,x)partial_theta3(t, x,para_new,nn), tspan,x0);
% [~, Xi2] = ode15s(@(t,x)partial_theta2(t, x,para_new,nn), tspan,x0);
% [~, Xi1] = ode15s(@(t,x)partial_theta1(t, x,para_new,nn), tspan,x0);
load('S_Matrix_Real_Tcell.mat')
mm = size(tspan1,1)-1;
Patial_phi1_theta = [reshape(Xi1(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi2(2:end,1:(nn+1))',[mm*(nn+1),1])...
    reshape(Xi3(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi4(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi5(2:end,1:(nn+1))',[mm*(nn+1),1])...
    reshape(Xi6(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi7(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi8(2:end,1:(nn+1))',[mm*(nn+1),1])...
    reshape(Xi9(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi10(2:end,1:(nn+1))',[mm*(nn+1),1])];
% 
% tspan = tspan2;
% y0 = yzero;
% x0 = zeros((nn+1)*3,1);
% x0 = [x0;y0];
% [~, Xi10] = ode15s(@(t,x)partial_theta10(t, x,para_new,nn), tspan,x0);
% [~, Xi9] = ode15s(@(t,x)partial_theta9(t, x,para_new,nn), tspan,x0);
% [~, Xi8] = ode15s(@(t,x)partial_theta8(t, x,para_new,nn), tspan,x0);
% [~, Xi7] = ode15s(@(t,x)partial_theta7(t, x,para_new,nn), tspan,x0);
% [~, Xi6] = ode15s(@(t,x)partial_theta6(t, x,para_new,nn), tspan,x0);
% [~, Xi5] = ode15s(@(t,x)partial_theta5(t, x,para_new,nn), tspan,x0);
% [~, Xi4] = ode15s(@(t,x)partial_theta4(t, x,para_new,nn), tspan,x0);
% [~, Xi3] = ode15s(@(t,x)partial_theta3(t, x,para_new,nn), tspan,x0);
% [~, Xi2] = ode15s(@(t,x)partial_theta2(t, x,para_new,nn), tspan,x0);
% [~, Xi1] = ode15s(@(t,x)partial_theta1(t, x,para_new,nn), tspan,x0);
load('S_Matrix_Real_Tumor.mat')
mm = size(tspan2,1)-1;
Patial_phi2_theta = [reshape(Xi1(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi2(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])...
    reshape(Xi3(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi4(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi5(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])...
    reshape(Xi6(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi7(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi8(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])...
    reshape(Xi9(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi10(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])];

% % 
 S = [Patial_phi1_theta;Patial_phi2_theta]; %Sensitive Matrix
% % S = Patial_phi1_theta;
% % S = Patial_phi2_theta;
% % S = Patial_phi3_theta;
% % S = Patial_phi4_theta;
  F = S'* S;
[U,Sigma,~]=svd(F);
r = min(find(diag(Sigma)<1e-5));
Ur = U(:,r:end);

%
fig3=figure(3);
clf();
% set(gcf,'Position',[169,223,987,541])

x=1:10;
bar(x,[(Sigma(1,1)),(Sigma(2,2)),(Sigma(3,3)),...
    (Sigma(4,4)),(Sigma(5,5)),(Sigma(6,6)),...
    (Sigma(7,7)),(Sigma(8,8)),(Sigma(9,9)),(Sigma(10,10))],'FaceColor',[0,0,0])

hold on
plot([0 11],[1e-5 1e-5],'LineStyle','--','LineWidth',1.5,'Color','k')
ylim([1e-24,1e2])
xlim([0 11])
ylabel('Eigenvalue of FIM');
set(gca,'xticklabel',{'U_1^T\theta','U_2^T\theta','U_3^T\theta'...
    'U_4^T\theta','U_5^T\theta','U_6^T\theta'...
    'U_7^T\theta','U_8^T\theta','U_9^T\theta','U_{10}^T\theta'},'YScale','log')

box off
set(gca,'FontName','Helvetica','FontSize',12,'FontWeight','bold','linewidth',1.2)


%% plot 


fig2=figure(2);
clf();
set(gcf,'Position',[169,223,987,541])
subplot(2,2,1)
 %% RAU
nnn = size(S,1);
s1=S(:,1);
A = [S(:,10) S(:,2:9)];
ss1 = (eye(nnn)-A*pinv(A))*s1;

s2=S(:,2);
A = [S(:,1) S(:,10) S(:,3:9)];
ss2 = (eye(nnn)-A*pinv(A))*s2;

s3=S(:,3);
A = [S(:,1:2) S(:,10) S(:,4:9)];
ss3 = (eye(nnn)-A*pinv(A))*s3;

s4=S(:,4);
A = [S(:,1:3) S(:,10) S(:,5:9)];
ss4 = (eye(nnn)-A*pinv(A))*s4;

s5=S(:,5);
A = [S(:,1:4) S(:,10) S(:,6:9)];
ss5 = (eye(nnn)-A*pinv(A))*s5;

s6=S(:,6);
A = [S(:,1:5) S(:,10) S(:,7:9)];
ss6 = (eye(nnn)-A*pinv(A))*s6;

s7=S(:,7);
A = [S(:,1:6) S(:,10) S(:,8:9)];
ss7 = (eye(nnn)-A*pinv(A))*s7;

s8=S(:,8);
A = [S(:,1:7) S(:,10) S(:,9:9)];
ss8 = (eye(nnn)-A*pinv(A))*s8;

s9=S(:,9);
A = [S(:,1:8) S(:,10:10)];
ss9 = (eye(nnn)-A*pinv(A))*s9;

s10=S(:,10);
A = S(:,1:9) ;
ss10 = (eye(nnn)-A*pinv(A))*s10;


SSS=[norm(ss1,"inf") norm(ss2,"inf") norm(ss3,"inf"),...
    norm(ss4,"inf") norm(ss5,"inf") norm(ss6,"inf"),...
    norm(ss7,"inf") norm(ss8,"inf") norm(ss9,"inf") norm(ss10,"inf")];

x=1:10;
bar(x,SSS,'FaceColor',[0,0,0])
ylim([1e-12,1])
ylabel('$$\|(I - AA^{\dagger}) \mathbf{s}_i\|_{\infty}$$', 'Interpreter', 'latex');
set(gca,'XTick',1:10,'xticklabel',{'\sigma','\rho','\eta','\mu','\epsilon','\beta_1','\beta_2','\phi','\lambda','\psi'}...
    ,'YScale','log')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
subplot(2,2,2)
alphaData = ones(10,10);  % 初始化为全不透明
alphaData(:, r:10) = 0.2;  % 设置右半部分透明度为 0.2

imagesc(1:10,1:10,abs(U),'AlphaData',alphaData);

% 设置 colormap
cmap = othercolor('BuDRd_12');
colormap(cmap);  % 可以选择其他 colormap 例如 'jet', 'hot', 'cool' 等
clim([0,1.0])
% 添加 colorbar 并设置标签
c = colorbar;
c.Label.String = '|\partial U_i^T\theta/\partial \theta_j|';  % 设置 colorbar 的标签
c.FontSize = 15;  % 调整字体大小
c.Label.FontWeight = 'bold';  % 设置字体加粗
set(gca,'XTick',1:10,'xticklabel',{'U_1','U_2','U_3','U_4','U_5','U_6','U_7','U_8','U_9','U_{10}'},...
    'YTick',1:10,'yticklabel',{'\sigma','\rho','\eta','\mu','\epsilon','\beta_1','\beta_2','\phi','\lambda','\psi'})
%title('\theta=[\sigma,\rho,\eta,\mu,\epsilon,\beta_1,\beta_2,\phi,\lambda,\psi]','FontSize',14,'FontWeight','bold')
%title('T & Tumor cells')

box off
set(gca,'FontName','Helvetica','FontSize',12,'FontWeight','bold','linewidth',1.2)


%% normalization

[para_new1,resnorm1,residual1,exitflag1,output1,lambda1,jacobian1] = ...
    lsqnonlin(@Objective_PDE,para_new,lb,ub,options,tspan1,tspan2,y0,yexp1,yexp2,100,Ur,para_new,nn);       

tspan = linspace(0,60*D1,101);
sol = ode45(@(t,y)PDE_eq(t,y,para_new1,nn),tspan,y0);
solpts = deval(sol,tspan);
sss1 = zeros(size(tspan,2),1);
for i = 1:size(tspan,2)
 sss1 (i) = sum(solpts(1:nn+1,i))/nn;
end
sss2 = zeros(size(tspan,2),1);
for i = 1:size(tspan,2)
sss2 (i) = sum(solpts(nn+2:2*(nn+1),i))/nn;
end
% subplot(1,2,1)
% plot(tspan*t0,ss1*E0,'k-')
% hold on
% plot(tspan1(2:end)*t0,yexp1(2:end)*E0,'o')
% subplot(1,2,2)
% plot(tspan*t0,ss2*T0,'k-')
% hold on
% plot(tspan2(2:end)*t0,yexp2(2:end)*T0,'o')
%% UQ
N1=101;
tspan = linspace(0,60*D1,101);
% y0 = yzero;
% x0 = zeros((nn+1)*3,1);
% x0 = [x0;y0];
% [~, Xi10] = ode15s(@(t,x)partial_theta10(t, x,para_new1,nn), tspan,x0);
% [~, Xi9] = ode15s(@(t,x)partial_theta9(t, x,para_new1,nn), tspan,x0);
% [~, Xi8] = ode15s(@(t,x)partial_theta8(t, x,para_new1,nn), tspan,x0);
% [~, Xi7] = ode15s(@(t,x)partial_theta7(t, x,para_new1,nn), tspan,x0);
% [~, Xi6] = ode15s(@(t,x)partial_theta6(t, x,para_new1,nn), tspan,x0);
% [~, Xi5] = ode15s(@(t,x)partial_theta5(t, x,para_new1,nn), tspan,x0);
% [~, Xi4] = ode15s(@(t,x)partial_theta4(t, x,para_new1,nn), tspan,x0);
% [~, Xi3] = ode15s(@(t,x)partial_theta3(t, x,para_new1,nn), tspan,x0);
% [~, Xi2] = ode15s(@(t,x)partial_theta2(t, x,para_new1,nn), tspan,x0);
% [~, Xi1] = ode15s(@(t,x)partial_theta1(t, x,para_new1,nn), tspan,x0);
load('S_Matrix_UQ_Real.mat')
mm = N1-1;
Patial_phi1_theta = [reshape(Xi1(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi2(2:end,1:(nn+1))',[mm*(nn+1),1])...
    reshape(Xi3(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi4(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi5(2:end,1:(nn+1))',[mm*(nn+1),1])...
    reshape(Xi6(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi7(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi8(2:end,1:(nn+1))',[mm*(nn+1),1])...
    reshape(Xi9(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi10(2:end,1:(nn+1))',[mm*(nn+1),1])];
Patial_phi2_theta = [reshape(Xi1(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi2(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])...
    reshape(Xi3(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi4(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi5(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])...
    reshape(Xi6(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi7(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi8(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])...
    reshape(Xi9(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi10(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])];

subplot(2,2,3)
%Ur = U(:,r:end);
%Ur = U;
S = Patial_phi1_theta;
SS = zeros(N1,10);
for i=1:N1-1
SS(i+1,:)=sum(S((i-1)*129+1:i*129,1:end),1)*(1/(nn));
end
sigma_th = 5e8; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (SS(i,:) * Ur)*(SS(i,:) * Ur)';
end
% 

x = tspan*t0;

Y_fit = sss1*E0;

CI=1.96*sqrt(Var_NN')*E0;
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(x, Y_fit, 'b', 'LineWidth', 2,'LineStyle','-')
hold on
plot(tspan1(2:end)*t0,yexp1(2:end)*E0,'bo','Markersize',8,'LineWidth',1.5)
ylabel('T cells')
xlabel('Time (days)')
lgd=legend('95% CI T cells','T cells','Real data');
 lgd.ItemTokenSize = [10,6];
lgd.Location = 'northwest';
lgd.Box = 'off';
lgd.FontWeight = 'bold';
box off
set(gca,'FontName','Helvetica','FontSize',12,'FontWeight','bold','linewidth',1.2)



 subplot(2,2,4)
% Ur = U(:,r:end);
%Ur = U;

S = Patial_phi2_theta;
SS = zeros(N1,10);
for i=1:N1-1
SS(i+1,:)=sum(S((i-1)*129+1:i*129,1:end),1)*(1/(nn));
end
sigma_th = 5e7; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (SS(i,:) * Ur)*(SS(i,:) * Ur)';
end
% 

x = tspan*t0;
Y_fit = sss2*T0;
CI=1.96*sqrt(Var_NN')*T0;
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(x, Y_fit, 'r', 'LineWidth', 2,'LineStyle','-')
hold on
plot(tspan2(2:end)*t0,yexp2(2:end)*T0,'ro','Markersize',8,'LineWidth',1.5)

ylabel("Tumor cells")
xlabel('Time (days)')
lgd=legend('95% CI Tumor','Tumor cells','Real data');

 lgd.ItemTokenSize = [10,6];
lgd.Location = 'northwest';
lgd.Box = 'off';
lgd.FontWeight = 'bold';
box off
set(gca,'FontName','Helvetica','FontSize',12,'FontWeight','bold','linewidth',1.2)

%%

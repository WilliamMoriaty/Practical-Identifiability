% File: PDE.m

% Parameters
para = zeros(10,1);
s = 1.36e4;   % reaction rate: s = para(1)
f = 0.2998e8;    % reaction rate: f = para(2)
g = 2.02e7;    % reaction rate: g = para(3)
d1 = 0.0412;   % reaction rate: d1 = para(4)
k1 = 1.3e-7;    % reaction rate: k1 = para(5)
k_1 = 24.0;    % reaction rate: k_1 = para(6)
k2 = 7.2;   % reaction rate: k2 = para(7)
p = 0.9997;    % reaction rate: p = para(8)
b1 = 0.18;    % reaction rate: b1 = para(9)
b2 = 2.0e-9;    % reaction rate: b2 = para(10)
% Nondimensional parameters
 E0=3.3e5;  % Effector T cells
 T0=0.5e9;   % Tumor cells
 C0=3.3e5; % Cell-complex 
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
y10(i) = hvd(disx(i)).*(1-exp(-1000*(disx(i)-l).^2));  % Initial E
y20(i) = (1-hvd(disx(i)).*(1-exp(-1000*(disx(i)-l).^2)));  % Initial T
y30(i) = (exp(-1000*(disx(i)-l).^2));  % Initial C

end

yzero = [y10;y20;y30];

% Time span
tspan = linspace(0,1000*D1,101);

% Solve the ODEs
[t, Y] = ode15s(@(t,y)PDE_eq(t,y,para,nn), tspan, yzero,para);
sol = ode15s(@(t,y)PDE_eq(t,y,para,nn),tspan,yzero);
solpts = deval(sol,tspan);
% Extract the data of E and T
y1 = Y(:, 1:nn+1);
y2 = Y(:, nn+1+1:2*(nn+1));
yexp1 = y1;
yexp2 = y2;
figure(1)
clf();
subplot(2,2,1)
plot(disx,yexp1(11,:))
xlabel('distance (x)')
ylabel('T cells(E)')
title('t=100 days')
subplot(2,2,2)
plot(disx,yexp1(41,:))
xlabel('distance (x)')
ylabel('T cells(E)')
title('t=400 days')
subplot(2,2,3)
plot(disx,yexp1(71,:))
xlabel('distance (x)')
ylabel('T cells(E)')
title('t=700 days')
subplot(2,2,4)
plot(disx,yexp1(101,:))
xlabel('distance (x)')
ylabel('T cells(E)')
title('t=1000 days')
figure(2)
clf();
subplot(2,2,1)
plot(disx,yexp2(11,:))
xlabel('distance (x)')
ylabel('Tumor cells(T)')
title('t=100 days')
subplot(2,2,2)
plot(disx,yexp2(41,:))
xlabel('distance (x)')
ylabel('Tumor cells(T)')
title('t=400 days')
subplot(2,2,3)
plot(disx,yexp2(71,:))
xlabel('distance (x)')
ylabel('Tumor cells(T)')
title('t=700 days')
subplot(2,2,4)
plot(disx,yexp2(101,:))
xlabel('distance (x)')
ylabel('Tumor cells(T)')
title('t=1000 days')

%% Fisher information matrix
tspan1=D1*[0 100 400 700 1000];
tspan = tspan1;
para_new = para;
y0 = yzero;
x0 = zeros((nn+1)*3,1);
x0 = [x0;y0];
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
load('S_Matrix.mat')
mm = size(tspan1,2)-1;
Patial_phi1_theta = [reshape(Xi1(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi2(2:end,1:(nn+1))',[mm*(nn+1),1])...
    reshape(Xi3(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi4(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi5(2:end,1:(nn+1))',[mm*(nn+1),1])...
    reshape(Xi6(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi7(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi8(2:end,1:(nn+1))',[mm*(nn+1),1])...
    reshape(Xi9(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi10(2:end,1:(nn+1))',[mm*(nn+1),1])];
Patial_phi2_theta = [reshape(Xi1(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi2(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])...
    reshape(Xi3(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi4(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi5(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])...
    reshape(Xi6(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi7(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi8(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])...
    reshape(Xi9(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi10(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])];
% 
 S = [Patial_phi1_theta;Patial_phi2_theta]; %Sensitive Matrix
% S = Patial_phi1_theta;
% S = Patial_phi2_theta;
% S = Patial_phi3_theta;
% S = Patial_phi4_theta;
  F = S'* S;
[U,Sigma,~]=svd(F);
i=mm-2;
Patial_phi1_theta_1 = [reshape(Xi1(2:i+1,1:(nn+1))',[i*(nn+1),1]) reshape(Xi2(2:i+1,1:(nn+1))',[i*(nn+1),1])...
    reshape(Xi3(2:i+1,1:(nn+1))',[i*(nn+1),1]) reshape(Xi4(2:i+1,1:(nn+1))',[i*(nn+1),1]) reshape(Xi5(2:i+1,1:(nn+1))',[i*(nn+1),1])...
    reshape(Xi6(2:i+1,1:(nn+1))',[i*(nn+1),1]) reshape(Xi7(2:i+1,1:(nn+1))',[i*(nn+1),1]) reshape(Xi8(2:i+1,1:(nn+1))',[i*(nn+1),1])...
    reshape(Xi9(2:i+1,1:(nn+1))',[i*(nn+1),1]) reshape(Xi10(2:i+1,1:(nn+1))',[i*(nn+1),1])];
Patial_phi2_theta_2 = [reshape(Xi1(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1]) reshape(Xi2(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1])...
    reshape(Xi3(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1]) reshape(Xi4(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1]) reshape(Xi5(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1])...
    reshape(Xi6(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1]) reshape(Xi7(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1]) reshape(Xi8(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1])...
    reshape(Xi9(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1]) reshape(Xi10(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1])];

S1 = Patial_phi1_theta_1; %Sensitive Matrix
F1 = S1'* S1;
[U1,Sigma1,~]=svd(F1);
S2 = Patial_phi2_theta_2; %Sensitive Matrix
F2 = S2'* S2;
[U2,Sigma2,~]=svd(F2);
r = min(find(diag(Sigma)<1e-5));
r1 = min(find(diag(Sigma1)<1e-5));
r2 = min(find(diag(Sigma2)<1e-5));

fig4=figure(4);
clf();
set(gcf,'Position',[18,477,1456,264])
subplot(1,3,1)
alphaData = ones(10,10);  
alphaData(:, r:10) = 0.2;  

imagesc(1:10,1:10,abs(U),'AlphaData',alphaData);

% 
cmap = othercolor('BuDRd_12');
colormap(cmap);  % 
clim([0,1.0])
% 
c = colorbar;
c.Label.String = '|\partial U_i^T\theta/\partial \theta_j|';  % 
c.FontSize = 15;  % 
c.Label.FontWeight = 'bold';  % 
set(gca,'XTick',1:10,'xticklabel',{'U_1','U_2','U_3','U_4','U_5','U_6','U_7','U_8','U_9','U_{10}'},...
    'YTick',1:10,'yticklabel',{'\sigma','\rho','\eta','\mu','\epsilon','\beta_1','\beta_2','\phi','\lambda','\psi'})
%title('\theta=[\sigma,\rho,\eta,\mu,\epsilon,\beta_1,\beta_2,\phi,\lambda,\psi]','FontSize',14,'FontWeight','bold')
title('T & Tumor cells','FontSize',12,'FontWeight','bold')
box off
set(gca,'FontName','Helvetica','FontSize',12,'FontWeight','bold','linewidth',1.2)

subplot(1,3,2)
alphaData = ones(10,10);  % 
alphaData(:, r1:10) = 0.2;  % 

imagesc(1:10,1:10,abs(U1),'AlphaData',alphaData);

% 
cmap = othercolor('BuDRd_12');
colormap(cmap);  % 
clim([0,1.0])
% 
c = colorbar;
c.Label.String = '|\partial U_i^T\theta/\partial \theta_j|';  % 
c.FontSize = 15;  % 
c.Label.FontWeight = 'bold';  % 
set(gca,'XTick',1:10,'xticklabel',{'U_1','U_2','U_3','U_4','U_5','U_6','U_7','U_8','U_9','U_{10}'},...
    'YTick',1:10,'yticklabel',{'\sigma','\rho','\eta','\mu','\epsilon','\beta_1','\beta_2','\phi','\lambda','\psi'})
%title('\theta=[\sigma,\rho,\eta,\mu,\epsilon,\beta_1,\beta_2,\phi,\lambda,\psi]','FontSize',14,'FontWeight','bold')
title('T cells','FontSize',12,'FontWeight','bold')
box off
set(gca,'FontName','Helvetica','FontSize',12,'FontWeight','bold','linewidth',1.2)

subplot(1,3,3)
alphaData = ones(10,10);  % 
alphaData(:, r2:10) = 0.2;  % 

imagesc(1:10,1:10,abs(U2),'AlphaData',alphaData);

% 
cmap = othercolor('BuDRd_12');
colormap(cmap);  % 
clim([0,1.0])
% 
c = colorbar;
c.Label.String = '|\partial U_i^T\theta/\partial \theta_j|';  % 
c.FontSize = 15;  % 
c.Label.FontWeight = 'bold';  % 
set(gca,'XTick',1:10,'xticklabel',{'U_1','U_2','U_3','U_4','U_5','U_6','U_7','U_8','U_9','U_{10}'},...
    'YTick',1:10,'yticklabel',{'\sigma','\rho','\eta','\mu','\epsilon','\beta_1','\beta_2','\phi','\lambda','\psi'})
%title('\theta=[\sigma,\rho,\eta,\mu,\epsilon,\beta_1,\beta_2,\phi,\lambda,\psi]','FontSize',14,'FontWeight','bold')
title('Tumor cells','FontSize',12,'FontWeight','bold')
box off
set(gca,'FontName','Helvetica','FontSize',12,'FontWeight','bold','linewidth',1.2)

%% plot 
% 

fig3=figure(3);
clf();
% set(gcf,'Position',[124,85,1291,650])
% subplot(2,2,3)
x=1:10;
bar(x,[(Sigma(1,1)),(Sigma(2,2)),(Sigma(3,3)),...
    (Sigma(4,4)),(Sigma(5,5)),(Sigma(6,6)),...
    (Sigma(7,7)),(Sigma(8,8)),(Sigma(9,9)),(Sigma(10,10));...
    (Sigma1(1,1)),(Sigma1(2,2)),(Sigma1(3,3)),...
    (Sigma1(4,4)),(Sigma1(5,5)),(Sigma1(6,6)),...
    (Sigma1(7,7)),(Sigma1(8,8)),(Sigma1(9,9)),(Sigma1(10,10));...
    (Sigma2(1,1)),(Sigma2(2,2)),(Sigma2(3,3)),...
    (Sigma2(4,4)),(Sigma2(5,5)),(Sigma2(6,6)),...
    (Sigma2(7,7)),(Sigma2(8,8)),(Sigma2(9,9)),(Sigma2(10,10))])

hold on
plot([0 11],[1e-5 1e-5],'LineStyle','--','LineWidth',1.5,'Color','k')
ylim([1e-18,1e10])
xlim([0 11])
ylabel('Eigenvalue of FIM');
set(gca,'xticklabel',{'U_1^T\theta','U_2^T\theta','U_3^T\theta'...
    'U_4^T\theta','U_5^T\theta','U_6^T\theta'...
    'U_7^T\theta','U_8^T\theta','U_9^T\theta','U_{10}^T\theta'},'YScale','log')
lgd=legend('T  & Tumor cells','T cells','Tumor cells');
lgd.ItemTokenSize = [10,6];
%lgd.Position = [0.746 0.301 0.139 0.326];
lgd.Box = 'off';
lgd.FontWeight = 'bold';
lgd.FontSize = 15;
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

fig4=figure(4);
clf();
% subplot(2,2,4)
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
A = [S(:,1:8) S(:,10)];
ss9 = (eye(nnn)-A*pinv(A))*s9;

s10=S(:,10);
A = S(:,1:9) ;
ss10 = (eye(nnn)-A*pinv(A))*s10;


SSS1=[norm(ss1,"inf") norm(ss2,"inf") norm(ss3,"inf"),...
    norm(ss4,"inf") norm(ss5,"inf") norm(ss6,"inf"),...
    norm(ss7,"inf") norm(ss8,"inf") norm(ss9,"inf") norm(ss10,"inf")];
%% S1
nnn = size(S1,1);
s1=S1(:,1);
A = [S1(:,10) S1(:,2:9)];
ss1 = (eye(nnn)-A*pinv(A))*s1;

s2=S1(:,2);
A = [S1(:,1) S1(:,10) S1(:,3:9)];
ss2 = (eye(nnn)-A*pinv(A))*s2;

s3=S1(:,3);
A = [S1(:,1:2) S1(:,10) S1(:,4:9)];
ss3 = (eye(nnn)-A*pinv(A))*s3;

s4=S1(:,4);
A = [S1(:,1:3) S1(:,10) S1(:,5:9)];
ss4 = (eye(nnn)-A*pinv(A))*s4;

s5=S1(:,5);
A = [S1(:,1:4) S1(:,10) S1(:,6:9)];
ss5 = (eye(nnn)-A*pinv(A))*s5;

s6=S1(:,6);
A = [S1(:,1:5) S1(:,10) S1(:,7:9)];
ss6 = (eye(nnn)-A*pinv(A))*s6;

s7=S1(:,7);
A = [S1(:,1:6) S1(:,10) S1(:,8:9)];
ss7 = (eye(nnn)-A*pinv(A))*s7;

s8=S1(:,8);
A = [S1(:,1:7) S1(:,10) S1(:,9:9)];
ss8 = (eye(nnn)-A*pinv(A))*s8;

s9=S1(:,9);
A = [S1(:,1:8) S1(:,10)];
ss9 = (eye(nnn)-A*pinv(A))*s9;

s10=S1(:,10);
A = S1(:,1:9) ;
ss10 = (eye(nnn)-A*pinv(A))*s10;


SSS2=[norm(ss1,"inf") norm(ss2,"inf") norm(ss3,"inf"),...
    norm(ss4,"inf") norm(ss5,"inf") norm(ss6,"inf"),...
    norm(ss7,"inf") norm(ss8,"inf") norm(ss9,"inf") norm(ss10,"inf")];
%% S2
nnn = size(S2,1);
s1=S2(:,1);
A = [S2(:,10) S2(:,2:9)];
ss1 = (eye(nnn)-A*pinv(A))*s1;

s2=S2(:,2);
A = [S2(:,1) S2(:,10) S2(:,3:9)];
ss2 = (eye(nnn)-A*pinv(A))*s2;

s3=S2(:,3);
A = [S2(:,1:2) S2(:,10) S2(:,4:9)];
ss3 = (eye(nnn)-A*pinv(A))*s3;

s4=S2(:,4);
A = [S2(:,1:3) S2(:,10) S2(:,5:9)];
ss4 = (eye(nnn)-A*pinv(A))*s4;

s5=S2(:,5);
A = [S2(:,1:4) S2(:,10) S2(:,6:9)];
ss5 = (eye(nnn)-A*pinv(A))*s5;

s6=S2(:,6);
A = [S2(:,1:5) S2(:,10) S2(:,7:9)];
ss6 = (eye(nnn)-A*pinv(A))*s6;

s7=S2(:,7);
A = [S2(:,1:6) S2(:,10) S2(:,8:9)];
ss7 = (eye(nnn)-A*pinv(A))*s7;

s8=S2(:,8);
A = [S2(:,1:7) S2(:,10) S2(:,9:9)];
ss8 = (eye(nnn)-A*pinv(A))*s8;

s9=S2(:,9);
A = [S2(:,1:8) S2(:,10)];
ss9 = (eye(nnn)-A*pinv(A))*s9;

s10=S2(:,10);
A = S2(:,1:9) ;
ss10 = (eye(nnn)-A*pinv(A))*s10;


SSS3=[norm(ss1,"inf") norm(ss2,"inf") norm(ss3,"inf"),...
    norm(ss4,"inf") norm(ss5,"inf") norm(ss6,"inf"),...
    norm(ss7,"inf") norm(ss8,"inf") norm(ss9,"inf") norm(ss10,"inf")];

SSS = [SSS1;SSS2;SSS3];
x=1:10;
bar(x,SSS)
% ylim([1e-12,1])
ylabel('$$\|(I - AA^{\dagger}) \mathbf{s}_i\|_{\infty}$$', 'Interpreter', 'latex');
set(gca,'XTick',1:10,'xticklabel',{'\sigma','\rho','\eta','\mu','\epsilon','\beta_1','\beta_2','\phi','\lambda','\psi'}...
    ,'YScale','log')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)


set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

box off

fig5=figure(5);
clf();
% subplot(2,2,1)
tspan1=D1*[0 100 400 700 1000];
mm = size(tspan1,2);
minsvd1=zeros(1,mm-1);
minsvd2=zeros(1,mm-1);
 minsvd12=zeros(1,mm-1);


maxsvd1=zeros(1,mm-1);
maxsvd2=zeros(1,mm-1);
 maxsvd12=zeros(1,mm-1);

for i=1:mm-1

Patial_phi1_theta_1 = [reshape(Xi1(2:i+1,1:(nn+1))',[i*(nn+1),1]) reshape(Xi2(2:i+1,1:(nn+1))',[i*(nn+1),1])...
    reshape(Xi3(2:i+1,1:(nn+1))',[i*(nn+1),1]) reshape(Xi4(2:i+1,1:(nn+1))',[i*(nn+1),1]) reshape(Xi5(2:i+1,1:(nn+1))',[i*(nn+1),1])...
    reshape(Xi6(2:i+1,1:(nn+1))',[i*(nn+1),1]) reshape(Xi7(2:i+1,1:(nn+1))',[i*(nn+1),1]) reshape(Xi8(2:i+1,1:(nn+1))',[i*(nn+1),1])...
    reshape(Xi9(2:i+1,1:(nn+1))',[i*(nn+1),1]) reshape(Xi10(2:i+1,1:(nn+1))',[i*(nn+1),1])];
Patial_phi2_theta_2 = [reshape(Xi1(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1]) reshape(Xi2(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1])...
    reshape(Xi3(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1]) reshape(Xi4(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1]) reshape(Xi5(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1])...
    reshape(Xi6(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1]) reshape(Xi7(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1]) reshape(Xi8(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1])...
    reshape(Xi9(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1]) reshape(Xi10(2:i+1,(nn+1)+1:2*(nn+1))',[i*(nn+1),1])];

S = Patial_phi1_theta_1;
F = S'* S;
minsvd1(i)=min(svd(F));
maxsvd1(i)=max(svd(F));
S = Patial_phi2_theta_2; %Sensitive Matrix
F = S'* S;
minsvd2(i)=min(svd(F));
maxsvd2(i)=max(svd(F));
S = [Patial_phi1_theta;Patial_phi2_theta]; %Sensitive Matrix
F = S'* S;
minsvd12(i)=min(svd(F));
maxsvd12(i)=max(svd(F));

end

plot(1:mm-1,minsvd1./maxsvd1*1e20,'o-','LineWidth',1.2,'MarkerSize',8)
hold on
plot(1:mm-1,minsvd2./maxsvd2*1e20,'*-','LineWidth',1.2,'MarkerSize',8)
hold on
plot(1:mm-1,minsvd12./maxsvd12*1e20,'LineStyle','-','Marker','square','LineWidth',1.2,'MarkerSize',8)
ylim([0 2.5])
% hold on
xlabel('Number of time points')
ylabel('\xi \times10^{-20}')
 set(gca,'Xtick',1:4)
lgd = legend("T cells","Tumor cells","T & Tumor cells");
lgd.ItemTokenSize = [10,6];
%lgd.Position = [0.746 0.301 0.139 0.326];
lgd.Box = 'off';
lgd.FontWeight = 'bold';
lgd.FontSize = 15;
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)


%% UQ
fig6=figure(6);
clf();
% subplot(2,2,2)
N1=101;
tspan = linspace(0,1000*D1,N1);
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
load('S_Matrix_UQ.mat')
mm = N1-1;
Patial_phi1_theta = [reshape(Xi1(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi2(2:end,1:(nn+1))',[mm*(nn+1),1])...
    reshape(Xi3(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi4(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi5(2:end,1:(nn+1))',[mm*(nn+1),1])...
    reshape(Xi6(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi7(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi8(2:end,1:(nn+1))',[mm*(nn+1),1])...
    reshape(Xi9(2:end,1:(nn+1))',[mm*(nn+1),1]) reshape(Xi10(2:end,1:(nn+1))',[mm*(nn+1),1])];
Patial_phi2_theta = [reshape(Xi1(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi2(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])...
    reshape(Xi3(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi4(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi5(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])...
    reshape(Xi6(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi7(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi8(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])...
    reshape(Xi9(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1]) reshape(Xi10(2:end,(nn+1)+1:2*(nn+1))',[mm*(nn+1),1])];


Ur = U(:,r:end);
%Ur = U;
S = Patial_phi1_theta;
SS = zeros(N1,10);
for i=1:N1-1
SS(i+1,:)=sum(S((i-1)*129+1:i*129,1:end),1)*(1/(nn));
end
sigma_th = 3e8; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (SS(i,:) * Ur)*(SS(i,:) * Ur)';
end
% 
yyaxis left
x = tspan*t0;

Y_fit = sum(yexp1,2)/nn;

CI=1.96*sqrt(Var_NN');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(x, Y_fit, 'Color',[0.07,0.62,1.00], 'LineWidth', 2,'LineStyle','-'); 
ylabel('T cells')

% subplot(1,2,2)
Ur = U(:,r:end);
%Ur = U;
S = Patial_phi2_theta;
SS = zeros(N1,10);
for i=1:N1-1
SS(i+1,:)=sum(S((i-1)*129+1:i*129,1:end),1)*(1/(nn));
end
% sigma_th = 1e5; % parameter estimation
Var_NN = zeros(1,N1);
for i = 1:N1
Var_NN(i)= sigma_th * (SS(i,:) * Ur)*(SS(i,:) * Ur)';
end
% 
yyaxis right
x = tspan*t0;
Y_fit = sum(yexp2,2)/nn;
CI=1.96*sqrt(Var_NN');
% 
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(x, Y_fit, 'r', 'LineWidth', 2,'LineStyle','-'); 
ylabel("Tumor cells")
xlabel('Time (days)')
lgd=legend('95% CI Tcell','T cells','95% CI Tumor','Tumor cells');
lgd.Box='off';
lgd.FontSize = 15;
lgd.ItemTokenSize = [10,6];
lgd.Position = [0.654,0.801,0.097,0.122];

box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

% 
% % 
% % 
% % 
% % 
% % 

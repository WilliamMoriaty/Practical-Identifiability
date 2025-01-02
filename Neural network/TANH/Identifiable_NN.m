N0 = 20;
N1 = 200;

%%
fig2=figure(2);
clf();
set(gcf,"Position",[211,517,1010,274])

%%
subplot(1,3,3)
alpha = load("Nalpha_30.txt");
beta = load("Nbeta_30.txt");
omega = load("Nomega_30.txt");

t = linspace(0,1,N0);
y = sin(2*pi*t);

t1 = linspace(0,1,N1);
phi = zeros(N1,1);
for i = 1:N1
phi(i) = NN(t1(i),alpha,beta,omega);
end

% S0: 20 data point
S0 = Sens_NN(alpha,beta,omega,N0);
F1 = S0'*S0;
[U,Sigma1,~]=svd(F1);
thresh = 1e-12;
r = min(find(diag(Sigma1)<thresh));
Ur = U(:,r:end);
Sigma13 = Sigma1;
% Ur = U;
m = size(alpha,1);
sigma_th = 5e9; % parameter estimation
Var_NN3 = zeros(1,N1);
for i = 1:N1
S_alpha = zeros(1,m);
S_beta = zeros(1,m);
S_omega = zeros(1,m);
    for j= 1:m
        S_alpha (j) = omega(j) * (1-(tanh(alpha(j)*((i-1)/(N1-1))+beta(j)))^2) * (i-1)/(N1-1);
        S_beta (j) = omega(j) * (1-(tanh(alpha(j)*((i-1)/(N1-1))+beta(j)))^2);
        S_omega (j) = tanh(alpha(j)*((i-1)/(N1-1))+beta(j));        
    end
    S = [S_alpha S_beta S_omega];
    Var_NN3(i)= sigma_th * (S * Ur)*(S * Ur)';
end

x = t1;
Y_fit = phi;
CI=1.96*sqrt(Var_NN3');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 

hold on
plot(t1,phi,'k-','LineWidth',1.2)
hold on
plot(t,y,'ko','MarkerSize',6,'LineWidth',1.2)
xlabel('t')
ylabel('h')
ylim([-1 1.5])
title("M=30",'FontWeight','bold')

box off

set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

box off
subplot(1,3,1)

alpha = load("Nalpha_10.txt");
beta = load("Nbeta_10.txt");
omega = load("Nomega_10.txt");

t = linspace(0,1,N0);
y = sin(2*pi*t);

t1 = linspace(0,1,N1);
phi = zeros(N1,1);
for i = 1:N1
phi(i) = NN(t1(i),alpha,beta,omega);
end

% S0: 20 data point
S0 = Sens_NN(alpha,beta,omega,N0);
F1 = S0'*S0;
[U,Sigma1,~]=svd(F1);
Sigma11 = Sigma1;
thresh = 1e-12;
r = min(find(diag(Sigma1)<thresh));
Ur = U(:,r:end);
%Ur = U;
m = size(alpha,1);
sigma_th = 1e11; % parameter estimation
Var_NN1 = zeros(1,N1);
for i = 1:N1
S_alpha = zeros(1,m);
S_beta = zeros(1,m);
S_omega = zeros(1,m);
    for j= 1:m
        S_alpha (j) = omega(j) * (1-(tanh(alpha(j)*((i-1)/(N1-1))+beta(j)))^2) * (i-1)/(N1-1);
        S_beta (j) = omega(j) * (1-(tanh(alpha(j)*((i-1)/(N1-1))+beta(j)))^2);
        S_omega (j) = tanh(alpha(j)*((i-1)/(N1-1))+beta(j));        
    end
    S = [S_alpha S_beta S_omega];
    Var_NN1(i)= sigma_th * (S * Ur)*(S * Ur)';
end

x = t1;
Y_fit = phi;
CI=1.96*sqrt(Var_NN1');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 

hold on
plot(t1,phi,'k-','LineWidth',1.2)
hold on
plot(t,y,'ko','MarkerSize',6,'LineWidth',1.2)
xlabel('t')
ylabel('h')
title("M=10",'FontWeight','bold')
lgd = legend("95% CI","Neural Network","Data");
lgd.FontWeight = 'bold';
lgd.Box='off';
lgd.ItemTokenSize = [10,6];
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)


subplot(1,3,2)
alpha = load("Nalpha_20.txt");
beta = load("Nbeta_20.txt");
omega = load("Nomega_20.txt");

t = linspace(0,1,N0);
y = sin(2*pi*t);

t1 = linspace(0,1,N1);
phi = zeros(N1,1);
for i = 1:N1
phi(i) = NN(t1(i),alpha,beta,omega);
end

% S0: 20 data point
S0 = Sens_NN(alpha,beta,omega,N0);
F1 = S0'*S0;
[U,Sigma2,~]=svd(F1);
Sigma12 = Sigma2;
r3 = min(find(diag(Sigma2)<thresh));
Ur = U(:,r3:end);
m = size(alpha,1);
%sigma_th = 5e6; % parameter estimation
Var_NN2 = zeros(1,N1);
for i = 1:N1
S_alpha = zeros(1,m);
S_beta = zeros(1,m);
S_omega = zeros(1,m);
    for j= 1:m
        S_alpha (j) = omega(j) * (1-(tanh(alpha(j)*((i-1)/(N1-1))+beta(j)))^2) * (i-1)/(N1-1);
        S_beta (j) = omega(j) * (1-(tanh(alpha(j)*((i-1)/(N1-1))+beta(j)))^2);
        S_omega (j) = tanh(alpha(j)*((i-1)/(N1-1))+beta(j));        
    end
    S = [S_alpha S_beta S_omega];
    Var_NN2(i)= sigma_th * (S * Ur)*(S * Ur)';
end

x = t1;
Y_fit = phi;
CI=1.96*sqrt(Var_NN2');
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 

hold on
plot(t1,phi,'k-','LineWidth',1.2)
hold on
plot(t,y,'ko','MarkerSize',6,'LineWidth',1.2)
xlabel('t')
ylabel('h')
title("M=20",'FontWeight','bold')
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

fig3=figure(3);
clf();
subplot(1,3,1)
Sigma = diag(Sigma11);
bar(Sigma,'FaceColor',[0,0,0])
hold on
plot([0,30.5],[thresh,thresh],'LineStyle','--','LineWidth',1.5,'Color','k')
ylabel('Eigenvalue of FIM');
set(gca,'YScale','log')
title('M=10','FontSize',15,'FontWeight','bold')

set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

box off
subplot(1,3,2)
Sigma = diag(Sigma12);
bar(Sigma,'FaceColor',[0,0,0])
hold on
plot([0,60.5],[thresh,thresh],'LineStyle','--','LineWidth',1.5,'Color','k')
ylabel('Eigenvalue of FIM');
set(gca,'YScale','log')
title('M=20','FontSize',15,'FontWeight','bold')

set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

box off
subplot(1,3,3)
Sigma = diag(Sigma13);
bar(Sigma,'FaceColor',[0,0,0])
hold on
plot([0,90.5],[thresh,thresh],'LineStyle','--','LineWidth',1.5,'Color','k')
ylabel('Eigenvalue of FIM');
set(gca,'YScale','log')

set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
title('M=30','FontSize',15,'FontWeight','bold')
box off
%%
%fig=figure(1);
% clf();
% sigma_th = 1e4; % parameter estimation
% N0 = 20;
% N1 = 100;
% subplot(2,2,1)
% alpha = load("Nalpha_5.txt");
% beta = load("Nbeta_5.txt");
% omega = load("Nomega_5.txt");
% 
% t = linspace(0,1,N0);
% y = sin(2*pi*t);
% 
% t1 = linspace(0,1,N1);
% phi = zeros(N1,1);
% for i = 1:N1
% phi(i) = NN(t1(i),alpha,beta,omega);
% end
% 
% % S1: 500 data point
% S1 = Sens_NN(alpha,beta,omega,N1);
% F1 = S1'*S1;
% [U,Sigma1,~]=svd(F1);
% r = min(find(diag(Sigma1)<1e-6));
% Ur = U(:,r:end);
% m = size(alpha,1);
% %sigma_th = 5e6; % parameter estimation
% Var_NN = zeros(1,N1);
% for i = 1:N1
% S_alpha = zeros(1,m);
% S_beta = zeros(1,m);
% S_omega = zeros(1,m);
%     for j= 1:m
%         S_alpha (j) = omega(j) * (1-tanh(alpha(j)*((i-1)/(N1-1))+beta(j))^2) * (i-1)/(N1-1);
%         S_beta (j) = omega(j) * (1-tanh(alpha(j)*((i-1)/(N1-1))+beta(j))^2);
%         S_omega (j) = tanh(alpha(j)*((i-1)/(N1-1))+beta(j));        
%     end
%     S = [S_alpha S_beta S_omega];
%     Var_NN(i)= sigma_th * (S * Ur)*(S * Ur)';
% end
% 
% x = t1;
% Y_fit = phi;
% CI=1.96*sqrt(Var_NN');
% % 绘制置信区域（拟合曲线上下的置信区间）
% fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% % 
% 
% hold on
% plot(t1,phi,'k-')
% 
% xlabel('t')
% ylabel('y')
% title("M=5")
% lg = legend("95% CI","Neural Network");
% 
% box off
% 
% 
% %%
% subplot(2,2,2)
% alpha = load("Nalpha_10.txt");
% beta = load("Nbeta_10.txt");
% omega = load("Nomega_10.txt");
% 
% t = linspace(0,1,N0);
% y = sin(2*pi*t);
% 
% t1 = linspace(0,1,N1);
% phi = zeros(N1,1);
% for i = 1:N1
% phi(i) = NN(t1(i),alpha,beta,omega);
% end
% 
% % S1: 500 data point
% S1 = Sens_NN(alpha,beta,omega,N1);
% F1 = S1'*S1;
% [U,Sigma1,~]=svd(F1);
% r = min(find(diag(Sigma1)<1e-6));
% Ur = U(:,r:end);
% m = size(alpha,1);
% %sigma_th = 5e6; % parameter estimation
% Var_NN = zeros(1,N1);
% for i = 1:N1
% S_alpha = zeros(1,m);
% S_beta = zeros(1,m);
% S_omega = zeros(1,m);
%     for j= 1:m
%         S_alpha (j) = omega(j) * (1-tanh(alpha(j)*((i-1)/(N1-1))+beta(j))^2) * (i-1)/(N1-1);
%         S_beta (j) = omega(j) * (1-tanh(alpha(j)*((i-1)/(N1-1))+beta(j))^2);
%         S_omega (j) = tanh(alpha(j)*((i-1)/(N1-1))+beta(j));        
%     end
%     S = [S_alpha S_beta S_omega];
%     Var_NN(i)= sigma_th * (S * Ur)*(S * Ur)';
% end
% 
% x = t1;
% Y_fit = phi;
% CI=1.96*sqrt(Var_NN');
% % 绘制置信区域（拟合曲线上下的置信区间）
% fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% % 
% 
% hold on
% plot(t1,phi,'k-')
% 
% xlabel('t')
% ylabel('y')
% title("M=10")
% box off
% 
% %%
% subplot(2,2,3)
% alpha = load("Nalpha_20.txt");
% beta = load("Nbeta_20.txt");
% omega = load("Nomega_20.txt");
% 
% t = linspace(0,1,N0);
% y = sin(2*pi*t);
% 
% t1 = linspace(0,1,N1);
% phi = zeros(N1,1);
% for i = 1:N1
% phi(i) = NN(t1(i),alpha,beta,omega);
% end
% 
% % S1: 500 data point
% S1 = Sens_NN(alpha,beta,omega,N1);
% F1 = S1'*S1;
% [U,Sigma1,~]=svd(F1);
% r = min(find(diag(Sigma1)<1e-6));
% Ur = U(:,r:end);
% m = size(alpha,1);
% %sigma_th = 5e6; % parameter estimation
% Var_NN = zeros(1,N1);
% for i = 1:N1
% S_alpha = zeros(1,m);
% S_beta = zeros(1,m);
% S_omega = zeros(1,m);
%     for j= 1:m
%         S_alpha (j) = omega(j) * (1-tanh(alpha(j)*((i-1)/(N1-1))+beta(j))^2) * (i-1)/(N1-1);
%         S_beta (j) = omega(j) * (1-tanh(alpha(j)*((i-1)/(N1-1))+beta(j))^2);
%         S_omega (j) = tanh(alpha(j)*((i-1)/(N1-1))+beta(j));        
%     end
%     S = [S_alpha S_beta S_omega];
%     Var_NN(i)= sigma_th * (S * Ur)*(S * Ur)';
% end
% 
% x = t1;
% Y_fit = phi;
% CI=1.96*sqrt(Var_NN');
% % 绘制置信区域（拟合曲线上下的置信区间）
% fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% % 
% 
% hold on
% plot(t1,phi,'k-')
% xlabel('t')
% ylabel('y')
% title("M=20")
% box off
% %%
% subplot(2,2,4)
% alpha = load("Nalpha_30.txt");
% beta = load("Nbeta_30.txt");
% omega = load("Nomega_30.txt");
% 
% t = linspace(0,1,N0);
% y = sin(2*pi*t);
% 
% t1 = linspace(0,1,N1);
% phi = zeros(N1,1);
% for i = 1:N1
% phi(i) = NN(t1(i),alpha,beta,omega);
% end
% 
% % S1: 500 data point
% S1 = Sens_NN(alpha,beta,omega,N1);
% F1 = S1'*S1;
% [U,Sigma1,~]=svd(F1);
% r = min(find(diag(Sigma1)<1e-6));
% Ur = U(:,r:end);
% %Ur = U;
% m = size(alpha,1);
% %sigma_th = 1e-2; % parameter estimation
% Var_NN = zeros(1,N1);
% for i = 1:N1
% S_alpha = zeros(1,m);
% S_beta = zeros(1,m);
% S_omega = zeros(1,m);
%     for j= 1:m
%         S_alpha (j) = omega(j) * (1-tanh(alpha(j)*((i-1)/(N1-1))+beta(j))^2) * (i-1)/(N1-1);
%         S_beta (j) = omega(j) * (1-tanh(alpha(j)*((i-1)/(N1-1))+beta(j))^2);
%         S_omega (j) = tanh(alpha(j)*((i-1)/(N1-1))+beta(j));        
%     end
%     S = [S_alpha S_beta S_omega];
%     Var_NN(i)= sigma_th * (S * Ur)*(S * Ur)';
% end
% 
% x = t1;
% Y_fit = phi;
% CI=1.96*sqrt(Var_NN');
% % 绘制置信区域（拟合曲线上下的置信区间）
% fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% % 
% 
% hold on
% plot(t1,phi,'k-')
% xlabel('t')
% ylabel('y')
% title("M=30")
% box off
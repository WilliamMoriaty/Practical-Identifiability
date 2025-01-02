fig1=figure(1);
clf();
set(gcf,'Position',[780,134,275,678])
% 设定参数
theta1 = 1.0; 
theta2 = 1.0;    
theta3 = 1.0;    
N1=3;
N2=50;
x = linspace(1, 3, N1);  % fit data
x1 = linspace(0,4,N2); %uncertainty data
subplot(3,1,2)
% 计算polynomial函数值 
y = poly_function(x, theta1,theta2,theta3);
y1 = poly_function(x1, theta1,theta2,theta3);
yy=zeros(N1,3);
for i=1:N1
yy(i,:)=poly_para(x(i));
end
% 绘图
plot(x(1:N1), y(1:N1),'ko','MarkerSize',8,'LineWidth',1.2);

hold on
plot(x1,y1,'k-','LineWidth',1.2)

x = linspace(1, 3, N1);
y = poly_function(x1, theta1,theta2,theta3);
for i=1:N1
yy(i,:)=poly_para(x(i));
end
[U,Sigma,~]=svd(yy'*yy);
Var = zeros(N2,1);
thresh = 1.0;
r = min(find(diag(Sigma)<thresh));
for i=1:N2
yy1 = poly_para(x1(i));

Var(i) = yy1*U(:,r:end)*U(:,r:end)'*yy1';
end

sigma_th=10;
Y_fit=y1(1:N2)';
CI=1.96*sqrt(sigma_th*Var);
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x1'; flipud(x1')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

xlabel('t');
ylabel('h(t,\theta)')
%ylabel({'$\varphi(t,\theta)$'},'Interpreter','latex');
lgd = legend("Data","Polynomial","95% CI");
lgd.FontWeight = 'bold';
lgd.Location = 'best';
lgd.Box='off';
lgd.ItemTokenSize = [10,6];
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off


subplot(3,1,3)
x = linspace(1, 3, N1);  % fit data
x1 = linspace(0,4,N2); %uncertainty data
% 计算polynomial函数值 
y = poly_function(x, theta1,theta2,theta3);
y1 = poly_function(x1, theta1,theta2,theta3);
[U,Sigma,V]=svd(yy'*yy);
Var = zeros(N2,1);
thresh = 1.0;
%r = min(find(diag(Sigma)<thresh));
r =1;
for i=1:N2
yy1 = poly_para(x1(i));

Var(i) = yy1*U(:,r:end)*U(:,r:end)'*yy1';
end

sigma_th=1e-1;
Y_fit=y1(1:N2)';
CI=1.96*sqrt(sigma_th*Var);
% 绘制置信区域（拟合曲线上下的置信区间）
fill([x1'; flipud(x1')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(x(1:N1), y(1:N1),'ko','MarkerSize',8,'LineWidth',1.2);

hold on
plot(x1,y1,'k-','LineWidth',1.2)

xlabel('t');
ylabel('h(t,\theta)')
%ylabel({'$\bm{\varphi(t,\theta)}$'},'Interpreter','latex');
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off

subplot(3,1,1)

alphaData = ones(3,3);  % 初始化为全不透明
alphaData(:, 3:3) = 0.05;  % 设置右半部分透明度为 0.2

imagesc(1:3,1:3,abs(U),'AlphaData',alphaData);

% 设置 colormap
cmap = othercolor('BuDRd_12');
colormap(cmap);  % 可以选择其他 colormap 例如 'jet', 'hot', 'cool' 等
clim([0,1.1])
% 添加 colorbar 并设置标签
c = colorbar;
c.Label.String = '|\partial U_i^T\theta/\partial \theta_j|';  % 设置 colorbar 的标签
c.FontSize = 12;  % 调整字体大小
c.Label.FontWeight = 'bold';  % 设置字体加粗
set(gca,'FontWeight','bold','FontSize',12)
set(gca,'xticklabel',{'U_1','U_2','U_3'},...
    'YTick',1:3,'yticklabel',{'\theta_1','\theta_2','\theta_3'})
% title('\theta=[\phi_1,\alpha,n]','FontSize',14,'FontWeight','bold')

figure(2)
clf();
set(gcf,'Position',[314,399,818,233])
subplot(1,2,1)

S = yy;
nn = size(S,1);
s1=S(:,1);
A = [S(:,3) S(:,2)];
ss1 = (eye(nn)-A*pinv(A))*s1;

s2=S(:,2);
A = [S(:,1) S(:,3)];
ss2 = (eye(nn)-A*pinv(A))*s2;

s3=S(:,3);
A = [S(:,1:2)];
ss3 = (eye(nn)-A*pinv(A))*s3;


SSS=[norm(ss1,"inf") norm(ss2,"inf") norm(ss3,"inf") ];
bar(1:3,SSS,'EdgeColor','none','FaceColor',[0,0,0])
hold on
ylim([1e-3,2])
ylabel('$$\|(I - AA^{\dagger}) \mathbf{s}_i\|_{\infty}$$', 'Interpreter', 'latex');
set(gca,'xticklabel',{'\theta_1','\theta_2','\theta_3'},'Yscale','log','ytick',[1e-3,1e-2,1e-1,1])
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off

subplot(1,2,2)
x=1:3;
bar(x,[(Sigma(1,1)),(Sigma(2,2)),(Sigma(3,3))],'FaceColor',[0,0,0])
hold on
plot([0,4],[1e-4,1e-4],'LineStyle','--','LineWidth',1.5,'Color','k')
ylim([1e-6,1e3])
xlim([0,4])
ylabel('Eigenvalue of FIM');
set(gca,'xticklabel',{'U_1^T\theta','U_2^T\theta','U_3^T\theta'},'YScale','log',...
    'ytick',[1e-4,1e-2,1e0,1e2])
lgd.ItemTokenSize = [10,6];

lgd.Box = 'off';
lgd.FontWeight = 'bold';
box off
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

box off
function y = poly_function(x, theta1,theta2,theta3)
   
    y = theta3*((x-1).*(x-2).*(x-3)+2)+theta2*x.^2+theta1;
end
function yy=poly_para(x)

yy=zeros(1,3);
yy(1)=1;
yy(2)=x.^2;
yy(3)=(x-1).*(x-2).*(x-3)+2;
end


Fig1=figure(1);
clf();
set(gcf,'Position',[246,186,798,579])
% 

Vmax = 1.0;  % 
Kd = 3.0;    % 
n = 8;       % 

N2=51;
VM = 10;      % maximum concentration

subplot(2,2,4)
x = 1.6:1.6:8.0;
x3 = [1.6,3.2,5.8];  % 
N1 = size(x,2);
x1 = linspace(0,VM,N2);
x2 = [x x3];
% 
y = hill_function(x, Vmax, Kd, n);
y1 = hill_function(x1, Vmax, Kd, n);
y3 = hill_function(x3, Vmax, Kd, n);
yy=zeros(N1,3);
for i=1:N1
yy(i,:)=hill_para(x(i),Vmax,Kd,n);
end
yy3=yy;
[U,Sigma,~]=svd(yy3'*yy3);
SSSigma = Sigma(3,3);
% plot
plot(x, y,'ko','markersize',8,'LineWidth',1.2);

hold on
plot(x1,y1,'k-','LineWidth',1.2)
hold on
plot(x3, y3,'r*','markersize',8,'LineWidth',1.2);
hold on
x = linspace(0, VM, N2);
y = hill_function(x, Vmax, Kd, n);
Var1 = zeros(1,N2-1);
for i=1:N2-1
yy=hill_para(x(i+1),Vmax,Kd,n);
Var1(i) = yy*U(:,3:3)*U(:,3:3)'*yy';
end

sigma=1e1;
Y_fit=[0;y(2:N2)'];
CI=[0;1.96*sqrt(sigma*Var1)'];
% 
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

xlabel(' ligand concentration (x)');
ylabel(['fraction of the receptor' sprintf('\n') 'bound by ligand (h)']);
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
lgd = legend("data","hill","critical data","95% CI");
lgd.FontWeight = 'bold';
lgd.Location = 'best';
lgd.Box='off';
lgd.ItemTokenSize = [10,6];


subplot(2,2,2)

alphaData = ones(3,3);  % 
alphaData(:, 3:3) = 0.2;  % 

imagesc(1:3,1:3,abs(U),'AlphaData',alphaData);

% set colormap
cmap = othercolor('BuDRd_12');
colormap(cmap);  % 
clim([0,1.1])
%xlim([1 3])
% add colorbar 
c = colorbar;
c.Label.String = '|\partial U_i^T\theta/\partial \theta_j|';  % 
c.FontSize = 15;  % 
c.Label.FontWeight = 'bold';  % 
set(gca,'FontWeight','bold','FontSize',15)
set(gca,'xticklabel',{'U_1','U_2','U_3'},...
    'YTick',1:3,'yticklabel',{'V_{max}','K_d','n'})
% title('\theta=[\phi_1,\alpha,n]','FontSize',14,'FontWeight','bold')

%% RAU
subplot(2,2,3)
S = yy3;
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
% plot([0,4],[1e-2,1e-2],'LineStyle','--','LineWidth',1.5,'Color','k')
ylim([1e-3,2])
ylabel('$$\|(I - AA^{\dagger}) \mathbf{s}_i\|_{\infty}$$', 'Interpreter', 'latex');
set(gca,'xticklabel',{'V_{max}','K_d','n'},'Yscale','log','ytick',[1e-3,1e-2,1e-1,1])
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off

subplot(2,2,1)
x=1:3;
bar(x,[(Sigma(1,1)),(Sigma(2,2)),(Sigma(3,3))],'FaceColor',[0,0,0])
hold on
plot([0,5],[1e-4,1e-4],'LineStyle','--','LineWidth',1.5,'Color','k')
ylim([1e-8,1e3])
xlim([0,4])
ylabel('Eigenvalue of FIM');
set(gca,'xticklabel',{'U_1^T\theta','U_2^T\theta','U_3^T\theta'},'YScale','log',...
    'ytick',[1e-6,1e-4,1e-2,1,1e2])

set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

box off

%%
fig3=figure(3);
clf();
set(gcf,'Position',[298,525,319,243])
x = 1.6:1.6:8.0;
x3 = [3.2,5.8,0.4];  % 
N1 = size(x,2);
x1 = linspace(0,VM,N2);
x2 = [x x3];
% 
y = hill_function(x, Vmax, Kd, n);
y1 = hill_function(x1, Vmax, Kd, n);

% 
plot(x, y,'ko','markersize',8,'LineWidth',1.2);

hold on

plot(x1,y1,'k-','LineWidth',1.2)

x = linspace(0, VM, N2);
y = hill_function(x, Vmax, Kd, n);
Var2 = zeros(1,N2-1);
for i=1:N2-1
yy=hill_para(x(i+1),Vmax,Kd,n);
Var2(i) = yy*U(:,1:3)*U(:,1:3)'*yy';
end

sigma=1e-3;
Y_fit=[0;y(2:N2)'];
CI=[0;1.96*sqrt(sigma*Var2)'];
% 
fill([x'; flipud(x')], [Y_fit + CI; flipud(Y_fit - CI)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
xlabel(' ligand concentration (x)');
ylabel(['fraction of the receptor' sprintf('\n') 'bound by ligand (h)']);
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
lgd = legend("data","hill","95% CI");
lgd.FontWeight = 'bold';
lgd.Location = 'best';
lgd.Box='off';
lgd.ItemTokenSize = [10,6];

set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off

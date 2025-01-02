fig1=figure(1);
clf();
set(gcf,'Position',[298,525,957,243])
Vmax = 1.0;  % 最大反应速率
Kd = 3.0;    % 解离常数
n = 8;       % Hill系数
N1 = 6;
VM = 8;      % maximum concentration

N=51;
X1=[];X2=[];X3=[];
l1=zeros(21,1);
t = 1.6:1.6:8.0;
y = hill_function(t, Vmax, Kd, n);
options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',10000,'MaxIterations',20000,...
    'algorithm','levenberg-marquardt','StepTolerance',1e-8);
x11=1;delta=0.2;
for i=1:N
x1=linspace(x11-delta,x11+delta,N);
lb = [0 0];
ub = [3 10];
x0 = [3.2 8];
fun = @(x)Hill_loss(t,x1(i),x(1),x(2),y);
x = lsqnonlin(fun,x0,lb,ub,options);
X1=[X1;x];
l1(i)=Hill_loss(t,x1(i),x(1),x(2),y);

end

%plot(x1,l2)

l2=zeros(N,1);
x21=3;delta=0.3;
for i=1:N
x2=linspace(x21-delta,x21+delta,N);
lb = [0 6];
ub = [2 10];

x0 = [1.0 8];
fun = @(x)Hill_loss(t,x(1),x2(i),x(2),y);
x = lsqnonlin(fun,x0,lb,ub,options);
X2=[X2;x];
l2(i)=Hill_loss(t,x(1),x2(i),x(2),y);

end

% % 
l3=zeros(N,1);
x31=8;delta=1;
for i=1:N
x3=linspace(x31-delta,x31+delta,N);
lb = [1 1];
ub = [15 15];
x0 = [1.2 3.2];
fun = @(x)Hill_loss(t,x(1),x(2),x3(i),y);
x = lsqnonlin(fun,x0,lb,ub,options);
X3=[X3;x];
l3(i)=Hill_loss(t,x(1),x(2),x3(i),y);

end
% %L=loss(0,1,1);
subplot(1,3,1)
plot(x1,l1,'k-','LineWidth',1.2)
%ylim([-1e-6,0.01])
xlabel('V_{max}')
ylabel('$l(\hat{h},\theta|V_{max})$','Interpreter','latex')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
hold on
subplot(1,3,2)
plot(x2,l2,'k-','LineWidth',1.2)
% ylim([-0.01,0.01])
xlabel('K_d')
ylabel('$l(\hat{h},\theta|K_d)$','Interpreter','latex')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
xlim([2.7,3.3])
hold on
subplot(1,3,3)
plot(x3,l3,'k-','LineWidth',1.2)
% ylim([-1e-4,1e-4])
xlabel('n')
ylabel('$l(\hat{h},\theta|n)$','Interpreter','latex')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off

%%
% %% normalization
% fig2=figure(2);
% X1=[];X2=[];X3=[];
% l1=zeros(21,1);
% t = [1.14,1.6,2.29,3.2,4.8,6.4];
% y = hill_function(t, Vmax, Kd, n);
% for i=1:5
% yy=[yy;hill_para(t(i),Vmax,Kd,n)];
% end
% F = yy'*yy;
% [U,Sigma,V]=svd(F);
% lambda=1;
% x11=1;delta=0.2;
% for i=1:N
% x1=linspace(x11-delta,x11+delta,N);
% lb = [0 0];
% ub = [3 10];
% x0 = [3.2 8];
% fun = @(x)NHill_loss(t,x1(i),x(1),x(2),y,lambda,U(:,3));
% x = lsqnonlin(fun,x0,lb,ub,options);
% X1=[X1;x];
% l1(i)=NHill_loss(t,x1(i),x(1),x(2),y,lambda,U(:,3));
% 
% end
% 
% l2=zeros(N,1);
% x21=3;delta=0.3;
% for i=1:N
% x2=linspace(x21-delta,x21+delta,N);
% lb = [0 6];
% ub = [2 10];
% 
% x0 = [1.0 8];
% fun = @(x)NHill_loss(t,x(1),x2(i),x(2),y,lambda,U(:,3));
% x = lsqnonlin(fun,x0,lb,ub,options);
% X2=[X2;x];
% l2(i)=NHill_loss(t,x(1),x2(i),x(2),y,lambda,U(:,3));
% 
% end
% 
% % % 
% l3=zeros(N,1);
% x31=8;delta=1;
% for i=1:N
% x3=linspace(x31-delta,x31+delta,N);
% lb = [1 1];
% ub = [15 15];
% x0 = [1.2 3.2];
% fun = @(x)NHill_loss(t,x(1),x(2),x3(i),y,lambda,U(:,3));
% x = lsqnonlin(fun,x0,lb,ub,options);
% X3=[X3;x];
% l3(i)=NHill_loss(t,x(1),x(2),x3(i),y,lambda,U(:,3));
% 
% end
% % %L=loss(0,1,1);
% subplot(1,3,1)
% plot(x1,l1)
% %ylim([-1e-3,1e-3])
% xlabel('$\tilde{V}_{max}$','Interpreter','latex')
% ylabel('$\tilde{l}(y,\theta|V_{max})$','Interpreter','latex')
% hold on
% subplot(1,3,2)
% plot(x2,l2)
% % ylim([-0.01,0.01])
% xlabel('$\tilde{K}$','Interpreter','latex')
% ylabel('$\tilde{l}(y,\theta|K)$','Interpreter','latex')
% hold on
% subplot(1,3,3)
% plot(x3,l3)
% %ylim([-1e-3,1e-3])
% xlabel('$\tilde{n}$','Interpreter','latex')
% ylabel('$\tilde{l}(y,\theta|n)$','Interpreter','latex')
% 

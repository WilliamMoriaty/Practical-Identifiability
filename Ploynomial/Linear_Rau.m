N=51;
X1=[];X2=[];X3=[];
l1=zeros(21,1);
t = [1 2 3];
y = [4 7 12];
options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',10000,'MaxIterations',20000,...
    'algorithm','levenberg-marquardt','StepTolerance',1e-8);
x11=1;delta=0.2;
for i=1:N
x1=linspace(x11-delta,x11+delta,N);
lb = [-3 -3];
ub = [3 3];
x0 = [1.2 1.2];
fun = @(x)Linear_loss(t,x1(i),x(1),x(2),y);
x = lsqnonlin(fun,x0,lb,ub,options);
X1=[X1;x];
l1(i)=Linear_loss(t,x1(i),x(1),x(2),y);

end

%plot(x1,l2)

l2=zeros(N,1);
x21=1;delta=0.2;
for i=1:N
x2=linspace(x21-delta,x21+delta,N);
lb = [-3 -3];
ub = [3 3];

x0 = [1.2 1.2];
fun = @(x)Linear_loss(t,x(1),x2(i),x(2),y);
x = lsqnonlin(fun,x0,lb,ub,options);
X2=[X2;x];
l2(i)=Linear_loss(t,x(1),x2(i),x(2),y);

end

% % 
l3=zeros(N,1);
x31=1;delta=0.2;
for i=1:N
x3=linspace(x31-delta,x31+delta,N);
lb = [-3 -3];
ub = [3 3];
x0 = [1.2 1.2];
fun = @(x)Linear_loss(t,x(1),x(2),x3(i),y);
x = lsqnonlin(fun,x0,lb,ub,options);
X3=[X3;x];
l3(i)=Linear_loss(t,x(1),x(2),x3(i),y);

end
fig1=figure(1);
clf();
set(gcf,'Position',[314,399,818,233])
% %L=loss(0,1,1);
subplot(1,3,1)
plot(x1,l1,'k-','LineWidth',1.2)
ylim([-1e-3,1e-3])
xlabel('\theta_1')
ylabel('$l(\hat{h},\theta|\theta_1)$','Interpreter','latex')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
hold on
subplot(1,3,2)
plot(x2,l2,'k-','LineWidth',1.2)
% ylim([-0.01,0.01])
xlabel('\theta_2')
ylabel('$l(\hat{h},\theta|\theta_2)$','Interpreter','latex')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
hold on
subplot(1,3,3)
plot(x3,l3,'k-','LineWidth',1.2)
ylim([-1e-3,1e-3])
xlabel('\theta_3')
ylabel('$l(\hat{h},\theta|\theta_3)$','Interpreter','latex')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off

%%
fig2=figure(2);
clf();
set(gcf,'Position',[314,399,818,233])
N=51;
X1=[];X2=[];X3=[];
S = [1 1 2;1 4 2;1 9 2];
F = S'*S;
[U,~,~]=svd(F);
xx=[1;1;1];
l1=zeros(21,1);
x11=0;delta=0.2;
for i=1:N
x1=linspace(x11-delta,x11+delta,N);
lb = [-3 -3];
ub = [3 3];
x0 = [0 0];
fun = @(x)ULinear_loss(t,x1(i),x(1),x(2),y,U,xx);
x = lsqnonlin(fun,x0,lb,ub,options);
X1=[X1;x];
l1(i)=ULinear_loss(t,x1(i),x(1),x(2),y,U,xx);

end

%plot(x1,l2)

l2=zeros(N,1);
x21=0;delta=0.2;
for i=1:N
x2=linspace(x21-delta,x21+delta,N);
lb = [-3 -3];
ub = [3 3];
x0 = [0,0];
fun = @(x)ULinear_loss(t,x(1),x2(i),x(2),y,U,xx);
x = lsqnonlin(fun,x0,lb,ub,options);
X2=[X2;x];
l2(i)=ULinear_loss(t,x(1),x2(i),x(2),y,U,xx);

end

% % 
l3=zeros(N,1);
x31=0;delta=0.2;
for i=1:N
x3=linspace(x31-delta,x31+delta,N);
lb = [-3 -3];
ub = [3 3];
x0 = [0 0];
fun = @(x)ULinear_loss(t,x(1),x(2),x3(i),y,U,xx);
x = lsqnonlin(fun,x0,lb,ub,options);
X3=[X3;x];
l3(i)=ULinear_loss(t,x(1),x(2),x3(i),y,U,xx);

end
% %L=loss(0,1,1);
subplot(1,3,1)
plot(x1,l1,'k-','LineWidth',1.2)
%ylim([-1e-3,1e-3])
xlabel('U_1^T(\theta-\theta^*)')
ylabel('$l(\hat{h},\theta|U_1^T\theta)$','Interpreter','latex')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
subplot(1,3,2)
plot(x2,l2,'k-','LineWidth',1.2)
% ylim([-0.01,0.01])
xlabel('U_2^T(\theta-\theta^*)')
ylabel('$l(\hat{h},\theta|U_2^T\theta)$','Interpreter','latex')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
subplot(1,3,3)
plot(x3,l3,'k-','LineWidth',1.2)
ylim([-1e-3,1e-3])
xlabel('U_3^T(\theta-\theta^*)')
ylabel('$l(\hat{h},\theta|U_3^T\theta)$','Interpreter','latex')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
%% normalization
fig3=figure(3);
clf();
set(gcf,'Position',[314,399,818,233])
N=51;
X1=[];X2=[];X3=[];
l1=zeros(21,1);
t = [1 2 3];
y = [4 7 12];
S = [1 1 2;1 4 2;1 9 2];
F = S'*S;
[U,Sigma,V]=svd(F);
lambda=1;
options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',10000,'MaxIterations',20000,...
    'algorithm','levenberg-marquardt','StepTolerance',1e-8);
x11=1;delta=0.2;
for i=1:N
x1=linspace(x11-delta,x11+delta,N);
lb = [-3 -3];
ub = [3 3];
x0 = [1.2 1.2];
fun = @(x)NLinear_loss(t,x1(i),x(1),x(2),y,lambda,U(:,3));
x = lsqnonlin(fun,x0,lb,ub,options);
X1=[X1;x];
l1(i)=NLinear_loss(t,x1(i),x(1),x(2),y,lambda,U(:,3));

end

%plot(x1,l2)

l2=zeros(N,1);
x21=1;delta=0.2;
for i=1:N
x2=linspace(x21-delta,x21+delta,N);
lb = [-3 -3];
ub = [3 3];

x0 = [1.2 1.2];
fun = @(x)NLinear_loss(t,x(1),x2(i),x(2),y,lambda,U(:,3));
x = lsqnonlin(fun,x0,lb,ub,options);
X2=[X2;x];
l2(i)=NLinear_loss(t,x(1),x2(i),x(2),y,lambda,U(:,3));

end

% % 
l3=zeros(N,1);
x31=1;delta=0.2;
for i=1:N
x3=linspace(x31-delta,x31+delta,N);
lb = [-3 -3];
ub = [3 3];
x0 = [1.2 1.2];
fun = @(x)NLinear_loss(t,x(1),x(2),x3(i),y,lambda,U(:,3));
x = lsqnonlin(fun,x0,lb,ub,options);
X3=[X3;x];
l3(i)=NLinear_loss(t,x(1),x(2),x3(i),y,lambda,U(:,3));

end
% %L=loss(0,1,1);
subplot(1,3,1)
plot(x1,l1,'k-','LineWidth',1.2)
%ylim([-1e-3,1e-3])
% xlabel('$\tilde{\theta}_1$','Interpreter','latex')
% ylabel('$\tilde{l}(y,\theta|\theta_1)$','Interpreter','latex')
xlabel('\theta_1')
ylabel('$\bf {\hat{l}(\hat{h},\theta|\theta_1)}$','Interpreter','latex')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
subplot(1,3,2)
plot(x2,l2,'k-','LineWidth',1.2)
% ylim([-0.01,0.01])
% xlabel('$\tilde{\theta}_2$','Interpreter','latex')
% ylabel('$\tilde{l}(y,\theta|\theta_2)$','Interpreter','latex')
xlabel('\theta_2')
ylabel('$\bf {\hat{l}(\hat{h},\theta|\theta_2)}$','Interpreter','latex')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
subplot(1,3,3)
plot(x3,l3,'k-','LineWidth',1.2)
%ylim([-1e-3,1e-3])
xlabel('\theta_3')
ylabel('$\bf {\hat{l}(\hat{h},\theta|\theta_3)}$','Interpreter','latex')
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)
box off
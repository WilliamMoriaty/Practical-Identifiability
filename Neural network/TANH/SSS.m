alpha = load("Nalpha_30.txt");
beta = load("Nbeta_30.txt");
omega = load("Nomega_30.txt");
N0 = 20;
t = linspace(0,1,N0);
y = sin(2*pi*t);
N1 = N0;

fig1=figure(1);
set(gcf,"Position",[1,430,1506,284])

% % S1: 500 data point
S1 = Sens_NN(alpha,beta,omega,N0);
F1 = S1'*S1;
[U,Sigma1,~]=svd(F1);
%Ur = U;
m = size(alpha,1);

NN = 30;
nn = size(S1,1);
SS = [];
for j=1:3*NN-1
s1=S1(:,j);
A = [S1(:,1:j-1) S1(:,3*NN)  S1(:,j+1:3*NN-1)];
ss1 = (eye(nn)-A*pinv(A))*s1;
SS = [SS norm(ss1,"inf")];
end
s1=S1(:,3*NN);
A = S1(:,1:3*NN-1);
ss1 = (eye(nn)-A*pinv(A))*s1;
SS = [SS norm(ss1,"inf")];

S= [SS(1:NN);SS(NN+1:2*NN);SS(2*NN+1:3*NN)];
bar(1:NN,S)
% ylim([1e-3,2])
xlabel("Neuron index (i)")
ylabel('$$\|(I - AA^{\dagger}) \mathbf{s}_i\|_{\infty}$$', 'Interpreter', 'latex');
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','yscale','log','linewidth',1.2)
lgd = legend("\alpha_i","\beta_i","\omega_i");
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
lgd.Box='off';
lgd.ItemTokenSize = [10,6];
set(gca,'FontName','Helvetica','FontSize',15,'FontWeight','bold','linewidth',1.2)

box off
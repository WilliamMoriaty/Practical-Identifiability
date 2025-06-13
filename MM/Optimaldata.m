%% output t
t = [0,10,20]; 
para = zeros(3,1);
para(1) = 1e6;   % reaction rate: k1 = para(1)
para(2) = 1e-4;    % reaction rate: k2 = para(2)
para(3) = 0.1;    % reaction rate: k3 = para(3)

y10 = 5e-7;  % Initial y1
y20 = 2e-7;  % Initial y2
y30 = 0.0;  % Initial y3
y40 = 0.0;  % Initial y4
yzero = [y10 y20 y30 y40];
N = 6;
C = y10+y20+2*y30+y40;
para(1) = para(1)*C;
yzero = yzero/C;
y0 = yzero';

tspan=t;
x0 = zeros(4,1);
x0 =[x0;y0];
[~, Xi3] = ode15s(@(t,x)partial_theta3(t, x,para), tspan,x0);
[~, Xi2] = ode15s(@(t,x)partial_theta2(t, x,para), tspan,x0);
[~, Xi1] = ode15s(@(t,x)partial_theta1(t, x,para), tspan,x0);
Patial_phi4_theta = [Xi1(:,4) Xi2(:,4) Xi3(:,4)];
S = Patial_phi4_theta;
F = S'* S;
[U,Sigma,~]=svd(F);


dt=30;
eps = 5e-5;
M = 5;
m1 = size(t,2);
alpha = zeros(1,3);
sigmak = Sigma(3,3);
for m=m1+1:M+m1
T = max(t)+dt;
t = [t,T];

tspan=t;
x0 = zeros(4,1);
x0 =[x0;y0];
[~, Xi3] = ode15s(@(t,x)partial_theta3(t, x,para), tspan,x0);
[~, Xi2] = ode15s(@(t,x)partial_theta2(t, x,para), tspan,x0);
[~, Xi1] = ode15s(@(t,x)partial_theta1(t, x,para), tspan,x0);
Patial_phi4_theta = [Xi1(:,4) Xi2(:,4) Xi3(:,4)];
S = Patial_phi4_theta;

x = diag((S*U(:,end))'*S*U(:,end));
non_zero_indices = x ~= 0;
num_non_zero = sum(non_zero_indices);

if num_non_zero > 0
    break;
end

end
m = m-1;
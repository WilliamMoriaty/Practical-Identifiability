function J = Jacobian(para,y1,y2,y3,y4)
% Jacobian matrix \partial f/\partial x
J = zeros(4);
lambda_Ab = para(1);
lambda_tau = para(2);
lambda_Ntaup = para(3);
lambda_CN = para(4);
lambda_Ctau = para(5);
K_Ab = para(6);
K_taup = para(7);
K_N = para(8);
K_C = para(9);

  Ab = y1;
  tau_p = y2;
  N = y3;
  C = y4;
J(1,1) = lambda_Ab * (1 - 2 * Ab / K_Ab);

J(2,1) = lambda_tau * (1 - tau_p/K_taup );
J(2,2) = (-lambda_tau * Ab) / K_taup;


J(3,2) = lambda_Ntaup * (1 - N / K_N);
J(3,3) = (-1) * (lambda_Ntaup * tau_p) / K_N;

J(4,2) = lambda_Ctau * (1 - C / K_C);
J(4,3) = lambda_CN * (1 - C / K_C);
J(4,4) = (-1) * (lambda_CN * N + lambda_Ctau * tau_p) / K_C;

end
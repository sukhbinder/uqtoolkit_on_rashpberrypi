function [Urlz Vrlz Wrlz] = runNISPRealizations(quaddata)
% function that solves the system at the NISP quadrature points

Xi = quaddata.nodes;
Nq = quaddata.nquad;

% random parameter
mu = 20.75;
sigma = mu * 0.005;

% time mesh
Tf = 1000;
tym = 0 : 1e-2 : Tf;
Nt = length(tym);

% state variables
Urlz = zeros(Nq, Nt); 
Vrlz = zeros(Nq, Nt); 
Wrlz = zeros(Nq, Nt); 

h = waitbar(0, 'Running simulations, please wait ...');
for j = 1 : Nq 
   
   waitbar(j / Nq);

   xi = Xi(j, :);
   paramb = mu + sigma * xi;
  
   [X T] = solve_dae_det(paramb, Tf);
   X = X(1:Nt,:);
   Urlz(j,:) = X(:,1)';
   Vrlz(j,:) = X(:,2)';
   Wrlz(j,:) = X(:,3)';

end

close(h);

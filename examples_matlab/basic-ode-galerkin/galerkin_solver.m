function [t Y pcdata] = galerkin_solver
% Solve the ODE 
%
% y' = lambda * y^2, y(0) = 1
%
% with a random parameter lambda via Intrusive Spectral Projection

% time mesh:
Tfin = 5;
dt = Tfin / 100;
t = [0 : dt : Tfin];
ntyms = length(t);

% the basic PC setup
nord = 3;
ndim = 1;
pc_type = 'LEGENDRE';
pcdata = uq_pcset(nord, ndim, pc_type);

% the pce for lambda
lambda = zeros(pcdata.nPCTerms, 1);
lambda(1) = -3;       
lambda(2) =  1;

% allocate memory for pce of solution vector:
Y = zeros(ntyms, pcdata.nPCTerms);

% impose the initial cond:
Y(1, 1) = 1;            % Y(0) = 1 deterministic initial cond

k = 1;
while k < ntyms

   tmp = uq_product(pcdata, Y(k, :), Y(k, :));

   dudt = uq_product(pcdata, tmp, lambda);

   Y(k+1, :) = Y(k, :) + dt * dudt;
 
   k = k + 1;

end

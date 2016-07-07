% file: driver.m
%
% purpose: 
%    This is a driver script for galerkin_solver_2D;
%    calls galerkin_solver_2D.m to solve the ODE
%    and then performs some postprocessing. 


[t,Y,pcdata] = galerkin_solver_2D;

% compute the first and second moments:
[mu sigma2] = uq_meanvar(pcdata, Y);

% plot the mean +/- 2 std deviations
stddev = sqrt(sigma2);
plot(t, mu, 'linewidth', 2);
hold on;
plot(t, mu + 2 * stddev, '--r', t, mu - 2 * stddev, '--r');




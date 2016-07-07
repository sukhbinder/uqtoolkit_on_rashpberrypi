function odetest
% Consider the ODE:
%     y' = lambda * y^2, u(0) = 1
% we know that the  
% analytic sol is 
%
%         y(t) = 1 / (1 - lambda * t)
%
% this will be used in couple of tests of 
% the UQ Toolbox.


% numerical sol with ode45



lambda = -1;
frhs  = @(t, y)(lambda * y.^2);
ftrue = @(t)(1 ./ (1 - lambda * t));

[t,y] = ode45(frhs, [0 5], [1]);
plot(t, y);
hold on;
plot(t, ftrue(t), 'r--');

